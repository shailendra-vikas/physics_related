import numpy
from scipy.integrate import trapz
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
one_page_output=False

def  get_numbers(component,condition,err_component,beta_range):
    count=trapz(component[condition],beta_range[condition])
    count_err= numpy.sqrt(numpy.sum(err_component[condition]**2))*(beta_range[1]-beta_range[0])
    return (count,count_err)
   
def get_distribution_two_hist(abs_data,column,bins,count_dict):
    # get the 2halo term 
    two_halo_cond = abs_data['BETA'] > 0.05
    N_2halo=numpy.sum(two_halo_cond)
    property_dist_2halo=numpy.histogram(abs_data[column][two_halo_cond],bins=bins)[0]
    property_dist_2halo_norm,property_dist_2halo_norm_err = property_dist_2halo.astype('float')/N_2halo,numpy.sqrt(property_dist_2halo)/N_2halo

    one_plus_outflow_cond = numpy.logical_and(abs_data['BETA']> -.01,abs_data['BETA'] <= 0.04)
    N_one_plus_outflow=numpy.sum(one_plus_outflow_cond)
    property_data=abs_data[column][one_plus_outflow_cond]
    property_dist_one_plus_outflow = numpy.histogram(property_data,bins=bins)[0]
    property_dist_one_plus_outflow_norm = property_dist_one_plus_outflow.astype('float') *\
                 (count_dict['beta_n01_04_one_plus_outflow'] + count_dict['beta_n01_04_2halo'])/(count_dict['beta_n01_04_one_plus_outflow']*N_one_plus_outflow) \
                - (count_dict['beta_n01_04_2halo']/count_dict['beta_n01_04_one_plus_outflow'])* property_dist_2halo_norm

    property_dist_one_plus_outflow_norm_err = numpy.sqrt(property_dist_one_plus_outflow.astype('float') * \
                ((count_dict['beta_n01_04_one_plus_outflow'] + count_dict['beta_n01_04_2halo'])/(count_dict['beta_n01_04_one_plus_outflow']*N_one_plus_outflow))**2 +\
                ((count_dict['beta_n01_04_2halo']/count_dict['beta_n01_04_one_plus_outflow'])* property_dist_2halo_norm_err)**2)

    return property_dist_2halo_norm,property_dist_2halo_norm_err,property_dist_one_plus_outflow_norm,property_dist_one_plus_outflow_norm_err


def get_distribution_three_hist(abs_data,column, bins, count_dict, random=False):
    # get the 2halo term 
    two_halo_cond = abs_data['BETA'] > 0.05
    N_2halo=numpy.sum(two_halo_cond)
    property_dist_2halo=numpy.histogram(abs_data[column][two_halo_cond],bins=bins)[0]
    property_dist_2halo_norm,property_dist_2halo_norm_err = property_dist_2halo.astype('float')/N_2halo,numpy.sqrt(property_dist_2halo)/N_2halo

    cond_n01_04 = numpy.logical_and(abs_data['BETA']>-.01, abs_data['BETA'] <= 0.04)
    beta_data = abs_data['BETA'][cond_n01_04]
    property_data = abs_data[column][cond_n01_04]
    
    if random:
        numpy.random.shuffle(property_data)
    
    cond_01_04= numpy.logical_and(beta_data > 0.01,beta_data <= 0.04)
    N_01_04 = numpy.sum(cond_01_04)
    property_01_04_hist = numpy.histogram(property_data[cond_01_04],bins=bins)[0]
    property_dist_outflow_norm  = property_01_04_hist.astype('float')* \
                                  (count_dict['beta_01_04_2halo']+count_dict['beta_01_04_outflow'])/(count_dict['beta_01_04_outflow']*N_01_04) \
                                  - (count_dict['beta_01_04_2halo']/count_dict['beta_01_04_outflow']) * property_dist_2halo_norm

    property_dist_outflow_norm_err = numpy.sqrt(property_01_04_hist.astype('float') *\
                                    ((count_dict['beta_01_04_2halo']+count_dict['beta_01_04_outflow'])/(count_dict['beta_01_04_outflow']*N_01_04))**2 +\
                                    ((count_dict['beta_01_04_2halo']/count_dict['beta_01_04_outflow'])* property_dist_2halo_norm_err)**2)

    cond_n01_01 = numpy.logical_and(beta_data >-0.01, beta_data <= 0.01)
    N_n01_01 = numpy.sum(cond_n01_01)
    property_n01_01_hist = numpy.histogram(property_data[cond_n01_01],bins=bins)[0]
    property_dist_onehalo_norm  = property_n01_01_hist.astype('float')* \
                                (count_dict['beta_n01_01_2halo']+count_dict['beta_n01_01_outflow']+count_dict['beta_n01_01_1halo']) / (count_dict['beta_n01_01_1halo']*N_n01_01)\
                                -(count_dict['beta_n01_01_2halo']/count_dict['beta_n01_01_1halo']) * property_dist_2halo_norm \
                                -(count_dict['beta_n01_01_outflow']/count_dict['beta_n01_01_1halo']) * property_dist_outflow_norm
    property_dist_onehalo_norm_err = numpy.sqrt(property_n01_01_hist.astype('float') *\
                                    ((count_dict['beta_n01_01_2halo']+count_dict['beta_n01_01_outflow']+count_dict['beta_n01_01_1halo']) / (count_dict['beta_n01_01_1halo']*N_n01_01))**2 \
                                    + ((count_dict['beta_n01_01_2halo']/count_dict['beta_n01_01_1halo']) * property_dist_2halo_norm_err)**2 \
                                    + ((count_dict['beta_n01_01_outflow']/count_dict['beta_n01_01_1halo']) * property_dist_outflow_norm_err)**2)

    return property_dist_2halo_norm,property_dist_2halo_norm_err,property_dist_outflow_norm,property_dist_outflow_norm_err,property_dist_onehalo_norm,property_dist_onehalo_norm_err

def kstest_func(hist1,hist2):
    cum_sum_hist1 = numpy.cumsum(hist1)
    cum_sum_hist1= cum_sum_hist1.astype(float)/cum_sum_hist1[-1]
    cum_sum_hist2 = numpy.cumsum(hist2)
    cum_sum_hist2= cum_sum_hist2.astype(float)/cum_sum_hist2[-1]
    arg=  numpy.argmax(numpy.abs(cum_sum_hist1-cum_sum_hist2))
    return cum_sum_hist1[arg]-cum_sum_hist2[arg]

def Cramer_von_mises(hist1,hist2):
    cum_sum_hist1 = numpy.cumsum(hist1)
    cum_sum_hist1= cum_sum_hist1.astype(float)/cum_sum_hist1[-1]
    cum_sum_hist2 = numpy.cumsum(hist2)
    cum_sum_hist2= cum_sum_hist2.astype(float)/cum_sum_hist2[-1]
    return numpy.sum((cum_sum_hist1-cum_sum_hist2)**2) 
   
kstest = kstest_func 

def cumulate(hist1):
    cum_sum_hist1 = numpy.cumsum(hist1)
    cum_sum_hist1= cum_sum_hist1.astype(float)/cum_sum_hist1[-1]
    return cum_sum_hist1
    
def fwhm_ew_analysis(abs_file,info_file):
    # read the qso catalog and absorber catalog
    #abs_data = map(numpy.load,(abs_file,))
    abs_data = numpy.load(abs_file)

    beta_range=numpy.linspace(-.01,.2,1001)
    labels,values,cov,total_1halo,total_2halo,total_outflow,error_total_1halo,error_total_2halo,error_total_outflow,selections,data_bins,civ_hist,civ_hist_err=\
        numpy.load(info_file)

    abs_data_beta = (abs_data['VAC_Z'] - abs_data['Z']) * (2 + abs_data['VAC_Z'] + abs_data['Z'])/\
                            ((1+abs_data['VAC_Z'])**2 + (1+abs_data['Z'])**2)
    abs_data['BETA']=abs_data_beta
    cond_beta_2halo=abs_data_beta>.05
    cond_beta_2halo_range=beta_range >.05
    N_2halo=numpy.sum(cond_beta_2halo)
    
    beta_outflow_min,beta_outflow_max=.01,.04
    cond_outflow=numpy.logical_and(abs_data_beta> 0.01,abs_data_beta <= 0.04)
    cond_outflow_range=numpy.logical_and(beta_range> 0.01,beta_range <= 0.04)
    N_outflow=numpy.sum(cond_outflow)

    beta_onehalo_min,beta_onehalo_max=-0.01,0.01
    cond_onehalo=numpy.logical_and(abs_data_beta > -0.01, abs_data_beta <=0.01)
    cond_onehalo_range=numpy.logical_and(beta_range>-0.01,beta_range<=0.01)
    N_onehalo=numpy.sum(cond_onehalo)

    beta_one_plus_outflow_min,beta_one_plus_outflow_max= -.01,.04
    cond_one_plus_outflow=numpy.logical_and(abs_data_beta > -0.01, abs_data_beta <=0.04)
    cond_one_plus_outflow_range=numpy.logical_and(beta_range>-0.01,beta_range<=0.04)
    N_one_plus_outflow = numpy.sum(cond_one_plus_outflow_range)

    beta_05_2halo=trapz(total_2halo[cond_beta_2halo_range],beta_range[cond_beta_2halo_range])

    count_dict,count_err_dict={},{}
    count_dict['beta_01_04_2halo'],count_err_dict['beta_01_04_2halo'] =get_numbers(total_2halo,cond_outflow_range,error_total_2halo,beta_range)
    count_dict['beta_01_04_1halo'],count_err_dict['beta_01_04_1halo']=get_numbers(total_1halo,cond_outflow_range,error_total_1halo,beta_range)
    count_dict['beta_01_04_outflow'],count_err_dict['beta_01_04_outflow']=get_numbers(total_outflow,cond_outflow_range,error_total_outflow,beta_range)

    count_dict['beta_n01_01_2halo'],count_err_dict['beta_n01_01_2halo']=get_numbers(total_2halo,cond_onehalo_range,error_total_2halo,beta_range)
    count_dict['beta_n01_01_1halo'],count_err_dict['beta_n01_01_1halo']=get_numbers(total_1halo,cond_onehalo_range,error_total_1halo,beta_range)
    count_dict['beta_n01_01_outflow'],count_err_dict['beta_n01_01_outflow']=get_numbers(total_outflow,cond_onehalo_range,error_total_outflow,beta_range)

    count_dict['beta_n01_04_2halo'],count_err_dict['beta_n01_04_2halo']=get_numbers(total_2halo,cond_one_plus_outflow_range,error_total_2halo,beta_range)
    count_dict['beta_n01_04_one_plus_outflow'],count_err_dict['beta_n01_04_one_plus_outflow']=get_numbers(total_1halo+total_outflow,cond_one_plus_outflow_range,numpy.sqrt(error_total_1halo**2+error_total_outflow**2),beta_range)

    fwhm_bins=numpy.linspace(0,300,31)
    fwhm_bin_mid=0.5 * fwhm_bins[1:] + 0.5 * fwhm_bins[:-1]

    #=======================#
    twohalo_fwhm_hist=numpy.histogram(abs_data['R_FWHM'][cond_beta_2halo],bins=fwhm_bins)[0]
    twohalo_fwhm_hist_err=numpy.sqrt(twohalo_fwhm_hist)
    outflow_fwhm_hist=numpy.histogram(abs_data['R_FWHM'][cond_outflow],bins=fwhm_bins)[0]
    outflow_fwhm_hist_err=numpy.sqrt(outflow_fwhm_hist)
    onehalo_fwhm_hist=numpy.histogram(abs_data['R_FWHM'][cond_onehalo],bins=fwhm_bins)[0]
    onehalo_fwhm_hist_err=numpy.sqrt(onehalo_fwhm_hist)
    one_plus_outflow_hist=numpy.histogram(abs_data['R_FWHM'][cond_one_plus_outflow],bins=fwhm_bins)[0]
    one_plus_outflow_hist_err= numpy.sqrt(one_plus_outflow_hist)

    
    if one_page_output:
        pp = PdfPages('fwhm_ew.pdf')

    sum_onehalo,sum_twohalo,sum_outflow,sum_one_plus_outflow=float(numpy.sum(onehalo_fwhm_hist)),float(numpy.sum(twohalo_fwhm_hist)),float(numpy.sum(outflow_fwhm_hist)),float(numpy.sum(one_plus_outflow_hist))
    pyplot.figure()
    pyplot.errorbar(fwhm_bin_mid,twohalo_fwhm_hist/sum_twohalo,yerr=twohalo_fwhm_hist_err/sum_twohalo,fmt='o-',color='b',label=r'$\beta > 0.05$')
    pyplot.errorbar(fwhm_bin_mid,outflow_fwhm_hist/sum_outflow,yerr=outflow_fwhm_hist_err/sum_outflow,fmt='o-',color='g',label=r'$0.01 < \beta < 0.04$')
    pyplot.errorbar(fwhm_bin_mid,onehalo_fwhm_hist/sum_onehalo,yerr=onehalo_fwhm_hist_err/sum_onehalo,fmt='o-',color='r',label=r'$-0.01 < \beta < 0.01$')
    pyplot.ylabel('Fraction in bin of 10 km/s')
    pyplot.xlabel('FWHM km/s')
    pyplot.legend()
    if one_page_output:
        pp.savefig()
    else:
        pyplot.savefig('normed_fwhm_dist_3range.pdf',fmt='pdf')

    pyplot.figure()
    pyplot.errorbar(fwhm_bin_mid,twohalo_fwhm_hist/sum_twohalo,yerr=twohalo_fwhm_hist_err/sum_twohalo,fmt='o-',color='b',label=r'$\beta > 0.05$')
    pyplot.errorbar(fwhm_bin_mid,one_plus_outflow_hist/sum_one_plus_outflow,yerr=one_plus_outflow_hist_err/sum_one_plus_outflow,fmt='o-',color='g',label=r'$-0.01 < \beta < 0.04$')
    pyplot.ylabel('Fraction in bin of 10 km/s')
    pyplot.xlabel('FWHM km/s')
    pyplot.legend()
    if one_page_output:
        pp.savefig()
    else:
        pyplot.savefig('normed_fwhm_dist_2range.pdf',fmt='pdf')
    
    #=======================#
    print "For FWHM property"

    fwhm_dist_2halo_norm,fwhm_dist_2halo_norm_err,fwhm_dist_one_plus_outflow_norm,fwhm_dist_one_plus_outflow_norm_err = \
                            get_distribution_two_hist(abs_data,'R_FWHM',fwhm_bins,count_dict)
    pyplot.figure()
    pyplot.errorbar(fwhm_bin_mid,fwhm_dist_2halo_norm,yerr=fwhm_dist_2halo_norm_err,fmt='o-',color='b',label='Two halo')
    pyplot.errorbar(fwhm_bin_mid,fwhm_dist_one_plus_outflow_norm,yerr=fwhm_dist_one_plus_outflow_norm_err,fmt='o-',color='g',label='Outflow+ One halo')
    ymin,ymax= pyplot.ylim()
    pyplot.ylim((0,ymax))
    pyplot.ylabel('Fraction in bin of 10 km/s')
    pyplot.xlabel('FWHM km/s')
    pyplot.legend()
    if one_page_output:
        pp.savefig()
    else:
        pyplot.savefig('normed_fwhm_dist_2comp.pdf',fmt='pdf')
    print 'For two component the distance between the distribution is ',kstest(fwhm_dist_2halo_norm,fwhm_dist_one_plus_outflow_norm)

    fig1,fig2=pyplot.figure(),pyplot.figure()
    ax1, ax2 = fig1.add_subplot(111), fig2.add_subplot(111)
    fwhm_dist_2halo_data_norm,fwhm_dist_2halo_data_norm_err,fwhm_dist_outflow_data_norm,fwhm_dist_outflow_data_norm_err,fwhm_dist_1halo_data_norm,\
            fwhm_dist_1halo_data_norm_err = get_distribution_three_hist(abs_data,'R_FWHM',fwhm_bins,count_dict)
    ks_dist_distrib=[]
    for i in range(1000):
        fwhm_dist_2halo_norm,fwhm_dist_2halo_norm_err,fwhm_dist_outflow_norm,fwhm_dist_outflow_norm_err,fwhm_dist_1halo_norm,fwhm_dist_1halo_norm_err = \
            get_distribution_three_hist(abs_data,'R_FWHM',fwhm_bins,count_dict,random=True)
        ks_dist_distrib.append(kstest(fwhm_dist_outflow_norm,fwhm_dist_1halo_norm))
        ax1.plot(fwhm_bin_mid,fwhm_dist_outflow_norm,'-',color=(.8,.8,.8))
        ax2.plot(fwhm_bin_mid,fwhm_dist_1halo_norm,'-',color=(.8,.8,.8))
        if i==500:
            fwhm_dist_2halo_rand_norm,fwhm_dist_2halo_rand_norm_err,fwhm_dist_outflow_rand_norm,fwhm_dist_outflow_rand_norm_err,fwhm_dist_1halo_rand_norm,\
                fwhm_dist_1halo_rand_norm_err = fwhm_dist_2halo_norm,fwhm_dist_2halo_norm_err,fwhm_dist_outflow_norm,\
                fwhm_dist_outflow_norm_err,fwhm_dist_1halo_norm, fwhm_dist_1halo_norm_err 

    
    ax1.plot(fwhm_bin_mid,fwhm_dist_outflow_data_norm,'-',color='b',lw=2)
    ax2.plot(fwhm_bin_mid,fwhm_dist_1halo_data_norm,'-',color='b',lw=2)
    ax1.set_ylabel('Fraction in bin of 10 km/s')
    ax1.set_xlabel('FWHM km/s')
    ax2.set_ylabel('Fraction in bin of 10 km/s')
    ax2.set_xlabel('FWHM km/s')
    fig1.savefig('outflow_random_FWHM.eps',fmt='eps')
    fig2.savefig('onehalo_random_FWHM.eps',fmt='eps')
    
    
    
    
    pyplot.figure()
    pyplot.errorbar(fwhm_bin_mid,fwhm_dist_2halo_data_norm,yerr=fwhm_dist_2halo_data_norm_err,fmt='o-',color='b',label='Two halo')
    pyplot.errorbar(fwhm_bin_mid,fwhm_dist_outflow_data_norm,yerr=fwhm_dist_outflow_data_norm_err,fmt='o-',color='g',label='Outflow')
    pyplot.errorbar(fwhm_bin_mid,fwhm_dist_1halo_data_norm,yerr=fwhm_dist_1halo_data_norm_err,fmt='o-',color='r',label='One halo')
    ymin,ymax= pyplot.ylim()
    pyplot.ylim((0,ymax))
    pyplot.ylabel('Fraction in bin of 10 km/s')
    pyplot.xlabel('FWHM km/s')
    pyplot.legend()
    if one_page_output:
        pp.savefig()
    else:
        pyplot.savefig('normed_fwhm_dist_3comp.pdf',fmt='pdf')
    ks_dist_data=kstest(fwhm_dist_outflow_data_norm,fwhm_dist_1halo_data_norm)
    print "For three component distance between FWHM distribution",ks_dist_data
   
    pyplot.figure()
    pyplot.plot(fwhm_bin_mid,cumulate(fwhm_dist_outflow_data_norm),'-',color='b',label='Outflow Data') 
    pyplot.plot(fwhm_bin_mid,cumulate(fwhm_dist_1halo_data_norm),'--',color='b',label='One halo Data')
    pyplot.plot(fwhm_bin_mid,cumulate(fwhm_dist_outflow_rand_norm),'-',color='r',label='Outflow Rand') 
    pyplot.plot(fwhm_bin_mid,cumulate(fwhm_dist_1halo_rand_norm),'--',color='r',label='One halo Rand')
    print "Range of the distance", min(ks_dist_distrib), max(ks_dist_distrib)  
    ymin,ymax= pyplot.ylim()
    pyplot.ylim((0,ymax))
    pyplot.ylabel('Fraction in bin of 10 km/s')
    pyplot.xlabel('FWHM km/s')
    pyplot.legend(loc=0) 
    pyplot.savefig('fwhm_example.pdf',fmt='pdf')

    pyplot.figure()
    
    pyplot.hist(ks_dist_distrib,bins=numpy.linspace(min(ks_dist_distrib),max(ks_dist_distrib),41))
    pyplot.xlabel('ks distance')
    pyplot.ylabel('# in bin of .01')
    pyplot.axvline(ks_dist_data,ymin=0,ymax=0.7,lw=2,color='r',ls='--')
    pyplot.text(ks_dist_data,60,'distance for data %1.2f'%ks_dist_data,fontsize=20,color='r')
    pyplot.savefig('ks_dist_fwhm.pdf',fmt='pdf')

    print "For the EW property"
    ew_bins=numpy.linspace(.5,2.,31)
    ew_bin_mid=0.5 * ew_bins[1:] + 0.5 * ew_bins[:-1]

    ew_dist_2halo_norm,ew_dist_2halo_norm_err,ew_dist_one_plus_outflow_norm,ew_dist_one_plus_outflow_norm_err = \
                            get_distribution_two_hist(abs_data,'R_EW',ew_bins,count_dict)
    pyplot.figure()
    pyplot.errorbar(ew_bin_mid,ew_dist_2halo_norm,yerr=ew_dist_2halo_norm_err,fmt='o-',color='b',label='Two halo')
    pyplot.errorbar(ew_bin_mid,ew_dist_one_plus_outflow_norm,yerr=ew_dist_one_plus_outflow_norm_err,fmt='o-',color='g',label='Outflow + One halo')
    ymin,ymax= pyplot.ylim()
    pyplot.ylim((0,ymax))
    pyplot.ylabel(r'Fraction in bin of 0.05 $\AA$')
    pyplot.xlabel(r'EW $\AA$')
    pyplot.legend()
    if one_page_output:
        pp.savefig()
    else:
        pyplot.savefig('normed_ew_dist_2comp.pdf',fmt='pdf')

    print 'Distance between two components',kstest(fwhm_dist_2halo_norm,fwhm_dist_one_plus_outflow_norm)

    ew_dist_2halo_data_norm,ew_dist_2halo_data_norm_err,ew_dist_outflow_data_norm,ew_dist_outflow_data_norm_err,ew_dist_1halo_data_norm,ew_dist_1halo_data_norm_err = \
                            get_distribution_three_hist(abs_data,'R_EW',ew_bins,count_dict)
    pyplot.figure()
    pyplot.errorbar(ew_bin_mid,ew_dist_2halo_data_norm,yerr=ew_dist_2halo_data_norm_err,fmt='o-',color='b',label='Two halo')
    pyplot.errorbar(ew_bin_mid,ew_dist_outflow_data_norm,yerr=ew_dist_outflow_data_norm_err,fmt='o-',color='g',label='Outflow')
    pyplot.errorbar(ew_bin_mid,ew_dist_1halo_data_norm,yerr=ew_dist_1halo_data_norm_err,fmt='o-',color='r',label='one halo')
    ymin,ymax= pyplot.ylim()
    pyplot.ylim((0,ymax))
    pyplot.ylabel(r'Fraction in bin of 0.05 $\AA$')
    pyplot.xlabel(r'EW $\AA$')
    pyplot.legend()
    if one_page_output:
        pp.savefig()
    else:
        pyplot.savefig('normed_ew_dist_3comp.pdf',fmt='pdf')

    ks_dist_data=kstest(ew_dist_outflow_data_norm,ew_dist_1halo_data_norm)
    print "The distance for 3 component", ks_dist_data

    ks_dist_distrib=[]
    for i in range(1000):
        ew_dist_2halo_norm,ew_dist_2halo_norm_err,ew_dist_outflow_norm,ew_dist_outflow_norm_err,ew_dist_1halo_norm,ew_dist_1halo_norm_err = \
            get_distribution_three_hist(abs_data,'R_EW',ew_bins,count_dict,random=True)
        if i==500:
            ew_dist_2halo_rand_norm,ew_dist_2halo_rand_norm_err,ew_dist_outflow_rand_norm,ew_dist_outflow_rand_norm_err,ew_dist_1halo_rand_norm,\
                ew_dist_1halo_rand_norm_err = ew_dist_2halo_norm,ew_dist_2halo_norm_err,ew_dist_outflow_norm,\
                ew_dist_outflow_norm_err,ew_dist_1halo_norm, ew_dist_1halo_norm_err 

        ks_dist_distrib.append(kstest(ew_dist_outflow_norm,ew_dist_1halo_norm))

    pyplot.figure()
    pyplot.plot(ew_bin_mid,cumulate(ew_dist_outflow_data_norm),'-',color='b',label='Outflow Data') 
    pyplot.plot(ew_bin_mid,cumulate(ew_dist_1halo_data_norm),'--',color='b',label='One halo Data')
    pyplot.plot(ew_bin_mid,cumulate(ew_dist_outflow_rand_norm),'-',color='r',label='Outflow Rand') 
    pyplot.plot(ew_bin_mid,cumulate(ew_dist_1halo_rand_norm),'--',color='r',label='One halo Rand')
    print "Range of the distance distribution", min(ks_dist_distrib), max(ks_dist_distrib)  
    ymin,ymax= pyplot.ylim()
    pyplot.ylim((0,ymax))
    pyplot.ylabel(r'Fraction in bin of  0.05 $\AA$')
    pyplot.xlabel(r'EW $\AA$')
    pyplot.legend(loc=0) 
    pyplot.savefig('ew_example.pdf',fmt='pdf')

 
    pyplot.figure()
    pyplot.hist(ks_dist_distrib,bins=numpy.linspace( min(ks_dist_distrib),max(ks_dist_distrib),41))
    pyplot.xlabel('ks distance')
    pyplot.ylabel('# in bin of .01')
    pyplot.axvline(ks_dist_data,ymin=0,ymax=0.7,lw=2,color='r',ls='--')
    pyplot.text(ks_dist_data-.45,60,'distance for data %1.2f'%ks_dist_data,fontsize=20,color='r')
    pyplot.savefig('ks_dist_ew.pdf',fmt='pdf')
    
    if one_page_output:
        pp.close()


if __name__=='__main__':
    fwhm_ew_analysis('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_all.npy')
    

    
    

