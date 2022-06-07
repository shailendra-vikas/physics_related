import numpy
import minuit
import intrinsic_util
import correlation_provider
from numpy.random import multivariate_normal
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages



qso_dict = {""      : '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/qso_final1.npy',
            "low"   : '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/qso_final1_tri_lowlumin.npy',
            "medium": '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/qso_final1_tri_medlumin.npy',
            "high"  : '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/qso_final1_tri_highlumin.npy'}

abs_dict = {""      : '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1.npy',
            "low"   : '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_tri_lowlumin.npy',
            "medium": '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_tri_medlumin.npy',
            "high"  : '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_tri_highlumin.npy'}


def make_selection(qso_data,beta_range):
    '''Return a array which descirbe if the given quasar and given beta we could have any onserver or not.
    '''
    zmax,zmin=intrinsic_util.zab_limits(qso_data['Z'],beta=-.01)
    ZQSO,BETA_RANGE=numpy.meshgrid(qso_data['Z'],beta_range)
    ZMAX,BETA_RANGE=numpy.meshgrid(zmax,beta_range)
    ZMIN,BETA_RANGE=numpy.meshgrid(zmin,beta_range)
    Z=intrinsic_util.redshift_for_beta(ZQSO,BETA_RANGE)
    selection = numpy.where(
                        numpy.logical_and( Z < ZMAX, Z > ZMIN),1.0,0.0)
    selection = numpy.where(
                        numpy.logical_and( Z > intrinsic_util.remove_redshift_range1[0],
                                            Z < intrinsic_util.remove_redshift_range1[1]) ,
                        0.0,selection)
    selection = numpy.where(
                        numpy.logical_and( Z > intrinsic_util.remove_redshift_range2[0],
                                            Z < intrinsic_util.remove_redshift_range2[1]) ,
                        0.0,selection)
    return selection

def prepare_covarinace(covariance,labels):
    cov=numpy.empty((len(labels),len(labels)))
    #cov=[[covariance[lab1,lab2] for lab2 in labels] for lab1 in labels]
    cov=[[covariance.get((lab1,lab2),0.0) for lab2 in labels] for lab1 in labels]
    return numpy.array(cov)

def read_fit_all(info_file):
    labels,values,cov,total_1halo,total_2halo,total_outflow,error_total_1halo,error_total_2halo,error_total_outflow,selections,data_bins,civ_hist,civ_hist_err= \
            numpy.load(info_file)
    return labels,values,cov


def dP_dbeta_parameter_estimation(lum="", usefit=True, info_file=None):
    qso_file, abs_file = qso_dict[lum], abs_dict[lum]

    # read the qso catalog and absorber catalog
    qso_data,abs_data = map(numpy.load,(qso_file,abs_file))
    print 'No of quasar',len(qso_data)
    print 'No of absorber',len(abs_data)

    #prepare the data point
    #data_bins = numpy.concatenate((numpy.linspace(-0.01,.02,30,endpoint=False),numpy.linspace(0.02,.2,30)))
    offset=.05
    data_bins = numpy.logspace(numpy.log10(-.01+offset),numpy.log10(.2+offset))-offset
    data_bin_mids = 0.5 * data_bins[1:] + 0.5 * data_bins[:-1]

    abs_data_beta = (abs_data['VAC_Z'] - abs_data['Z']) * (2 + abs_data['VAC_Z'] + abs_data['Z'])/\
                        ((1+abs_data['VAC_Z'])**2 + (1+abs_data['Z'])**2)
    civ_hist,junk = numpy.histogram(abs_data_beta,bins = data_bins)
    civ_hist_err = numpy.sqrt(civ_hist) / (data_bins[1:] - data_bins[0:-1])
    civ_hist = civ_hist.astype('float') / (data_bins[1:] - data_bins[0:-1])

    #start preparing theory points
    beta_range=numpy.linspace(-.01,.2,1001) # the range for theory points
    selections= make_selection(qso_data,beta_range)

    ZQSO,BETA_RANGE=numpy.meshgrid(qso_data['Z'],beta_range)
    Z=intrinsic_util.redshift_for_beta(ZQSO,BETA_RANGE)
    R=intrinsic_util.fast_redshift_to_comoving_instance.distance(ZQSO)- \
                intrinsic_util.fast_redshift_to_comoving_instance.distance(Z)
    
    #precalculate for the single halo term
    log_one_plus_zqso=numpy.log(1+ZQSO)
    BETA_RANGE_SQ = BETA_RANGE**2
    

    #precalculation for the two halo term
    precalc_2halo = ((1+Z)**2 + (1+ZQSO)**2)**2 * intrinsic_util.cosmo.E_inv(Z) *intrinsic_util.cosmo.Dh/\
                            (4*intrinsic_util.cosmo.h * (1+ZQSO)**2)
    correlation_estimates = correlation_provider.get_corr(R,ZQSO)
    geometric_err = numpy.sqrt(correlation_estimates[1] * correlation_estimates[2])
    geometric_err = precalc_2halo * numpy.where(R < 0.0,0.0, geometric_err) * selections
    precalc_2halo = precalc_2halo * numpy.where(R < 0.0,0.0, 1. + correlation_estimates[0]) * selections
    log_one_plus_z = numpy.log(1+Z)

    weights=intrinsic_util.get_gaussian_window(0.003,beta_range[1]-beta_range[0])

    if usefit:
        # get values, error, covariance from fit
        print 'Fitting for the best fir parameter' 

        def minimization(n0_sigma0,epsilon,sigma_beta,f,beta_01,beta_02,beta_03):
            postcalc_1halo = numpy.log(f) - 0.5 * BETA_RANGE_SQ/sigma_beta**2 + epsilon * log_one_plus_zqso -.5 * numpy.log(2*numpy.pi)
            postcalc_2halo = n0_sigma0* numpy.exp((1+epsilon) * log_one_plus_z)
            total = numpy.exp(postcalc_1halo) * selections + postcalc_2halo * precalc_2halo +\
                        selections*numpy.interp(beta_range,(0.0,0.02/3.,0.04/3.,0.02,0.03),(0.0,beta_01,beta_02,beta_03,0.0)).reshape((len(beta_range),1))
            corr_err= numpy.sum(postcalc_2halo * geometric_err,axis=1)
            theory_points= numpy.sum(total,axis=1)
            theory_points= numpy.convolve(theory_points,weights,mode='same')        
            theory_points_to_compare=numpy.interp(data_bin_mids,beta_range,theory_points)
            error_theory_point=numpy.interp(data_bin_mids,beta_range,corr_err)
            total_err_sq= civ_hist_err**2 + error_theory_point**2
            return numpy.sum(numpy.where(total_err_sq > 0.0,(civ_hist-theory_points_to_compare)**2/total_err_sq,0.0))

        if info_file is not None:
            default_labels,default_value,default_cov =read_fit_all(info_file)
            default_values=dict(zip(default_labels,default_value))
        else:
            default_values = {'n0_sigma0' : 0.003, 'epsilon' : -2.0,  'sigma_beta' : 0.001, 'f' : 2000., 'beta_01' : 2., 'beta_02' : 1.3, 'beta_03' : 0.4}

        if lum == "":
            minimizer=minuit.Minuit(minimization, **default_values)
        else:
            #Only beta1,beta2,beta3 are variable
            minimizer=minuit.Minuit(minimization, fix_n0_sigma0= True, fix_epsilon=True, fix_sigma_beta=True, fix_f =True, **default_values)
            
        minimizer.maxcalls=10000
        minimizer.strategy=0
        #minimizer.printMode=1
        print "####", minimizer.fixed
        minimizer.migrad()
        minimizer.hesse()

        print "####",minimizer.values , minimizer.covariance
    
        labels=['n0_sigma0','epsilon','sigma_beta','f','beta_01','beta_02','beta_03']
        values=numpy.array([minimizer.values[lab] for lab in labels])
        cov=prepare_covarinace(minimizer.covariance,labels)
        print "####", cov

        if lum != "":
            print "Correcting the values and cov:"
            for i in range(0,4):
                values[i]=default_value[i]
            for i in range(0,4):
                for j in range(0,4):
                    cov[i,j] = default_cov[i,j]
    else:
        # get values, error, covariance from the file
        labels,values,cov =read_fit_all(info_file)

    print "Values:", values, cov
    print "Starting to sample from generating the function"
    param_values=multivariate_normal(values,cov,1000)
    total_1halo_all,total_2halo_all,total_outflow_all = [],[],[]
    for i in range(1000):
        param_dict=dict(zip(labels,param_values[i]))
        current_epsilon,current_f,current_sigma_beta,current_n0_sigma0,current_beta_02=\
            param_dict['epsilon'],param_dict['f'],param_dict['sigma_beta'],param_dict['n0_sigma0'],param_dict['beta_02']
        current_beta_01,current_beta_03=param_dict['beta_01'],param_dict['beta_03']
        postcalc_1halo = numpy.log(current_f) - 0.5 * BETA_RANGE_SQ/current_sigma_beta**2 + current_epsilon * log_one_plus_zqso -.5 * numpy.log(2*numpy.pi)
        total_1halo=numpy.sum(numpy.exp(postcalc_1halo) * selections ,axis=1)
        total_2halo=numpy.sum(current_n0_sigma0* numpy.exp((1+current_epsilon) * log_one_plus_z) * precalc_2halo, axis=1)
        total_outflow=numpy.sum(selections*numpy.interp(beta_range,(0.0,0.02/3.,0.04/3.,0.02,0.03),(0.0,current_beta_01,current_beta_02,current_beta_03,0.0)).reshape((len(beta_range),1)),axis=1)
        total_1halo=numpy.convolve(total_1halo,weights,mode='same')
        total_2halo=numpy.convolve(total_2halo,weights,mode='same')
        total_outflow=numpy.convolve(total_outflow,weights,mode='same')
        total_1halo_all.append(total_1halo)
        total_2halo_all.append(total_2halo)
        total_outflow_all.append(total_outflow)

    total_1halo_all,total_2halo_all,total_outflow_all =map(numpy.array,(total_1halo_all,total_2halo_all,total_outflow_all))
    error_total_1halo=numpy.array([numpy.std(total_1halo_all[:,i][numpy.logical_not(numpy.isnan(total_1halo_all[:,i]))]) for i in range(total_1halo_all.shape[1])])
    error_total_2halo=numpy.array([numpy.std(total_2halo_all[:,i][numpy.logical_not(numpy.isnan(total_2halo_all[:,i]))]) for i in range(total_2halo_all.shape[1])])
    error_total_outflow=numpy.array([numpy.std(total_outflow_all[:,i][numpy.logical_not(numpy.isnan(total_outflow_all[:,i]))]) for i in range(total_outflow_all.shape[1])])

    value_dict = dict(zip(labels,values))
    best_epsilon, best_f, best_sigma_beta, best_n0_sigma0= value_dict['epsilon'], value_dict['f'], value_dict['sigma_beta'], value_dict['n0_sigma0']
    best_beta_01, best_beta_03, best_beta_02 = value_dict['beta_01'], value_dict['beta_03'], value_dict['beta_02']
    postcalc_1halo = numpy.log(best_f) - 0.5 * BETA_RANGE_SQ/best_sigma_beta**2 + best_epsilon * log_one_plus_zqso -.5 * numpy.log(2*numpy.pi)
    total_1halo=numpy.sum(numpy.exp(postcalc_1halo) * selections ,axis=1)
    total_2halo=numpy.sum(best_n0_sigma0* numpy.exp((1+best_epsilon) * log_one_plus_z) * precalc_2halo, axis=1)
    total_outflow=numpy.sum(selections*numpy.interp(beta_range,(0.0,0.02/3.,0.04/3.,0.02,0.03),(0.0,best_beta_01,best_beta_02,best_beta_03,0.0)).reshape((len(beta_range),1)),axis=1)
    total_1halo=numpy.convolve(total_1halo,weights,mode='same')
    total_2halo=numpy.convolve(total_2halo,weights,mode='same')
    total_outflow=numpy.convolve(total_outflow,weights,mode='same')
    theory_points= total_1halo + total_2halo + total_outflow
    corr_err= numpy.sum(best_n0_sigma0* numpy.exp((1+best_epsilon) * log_one_plus_z) * geometric_err,axis=1)
    no_of_quasars_before_cov=numpy.sum(selections,axis=1)
    no_of_quasars=numpy.convolve(no_of_quasars_before_cov,weights,mode='same')
    theory_points_to_compare=numpy.interp(data_bin_mids,beta_range,theory_points)

    print 'Area under one halo:',numpy.sum(numpy.where(no_of_quasars>0,total_1halo/no_of_quasars,0.0)),' error:', numpy.sqrt(numpy.sum(numpy.where(no_of_quasars>0,error_total_1halo/no_of_quasars,0.0))**2)
    print 'Area under two halo:',numpy.sum(numpy.where(no_of_quasars>0,total_2halo/no_of_quasars,0.0)),' error:', numpy.sqrt(numpy.sum(numpy.where(no_of_quasars>0,error_total_2halo/no_of_quasars,0.0))**2)
    print 'Area under outflow:',numpy.sum(numpy.where(no_of_quasars>0,total_outflow/no_of_quasars,0.0)),' error:', numpy.sqrt(numpy.sum(numpy.where(no_of_quasars>0,error_total_outflow/no_of_quasars,0.0))**2)
    print ' Multiply by ',beta_range[1]-beta_range[0]


    numpy.save('fit_info%s.npy'%lum, (labels,values,cov,total_1halo,total_2halo,total_outflow,error_total_1halo,error_total_2halo,error_total_outflow,selections,data_bins,civ_hist,civ_hist_err))

    
def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option( "--lum", dest="lum", default="", help="luminosity option")
    parser.add_option( "--nofit", action="store_false", default=True, dest="usefit")
    parser.add_option( "--fitinfo", dest="fitinfo", default=None, help="fitinfo")
    parser.add_option( "--test", action="store_true",dest="test", default=False, help="test")
    
    (options, args) = parser.parse_args()
    if not options.usefit and options.fitinfo is None:
        print "Invalid parameter. When usefit is false the fitinfo must be valid."
        return
    if options.lum not in ["","low","high","medium"]:
        print "Invalid value of luminosity parameter. The valid values are 'low', 'high' and 'medium' when provided"

    fitinfo = options.fitinfo 
    if fitinfo is not None:
        fitinfo = "/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/" + fitinfo

    if options.test:
        print "Calling dP_dbeta_parameter_estimation with parameters", qso_dict[options.lum], abs_dict[options.lum], options.usefit , fitinfo
        return
    dP_dbeta_parameter_estimation(lum=options.lum, usefit=options.usefit, info_file=fitinfo)        
     

if __name__=='__main__':
    main()
