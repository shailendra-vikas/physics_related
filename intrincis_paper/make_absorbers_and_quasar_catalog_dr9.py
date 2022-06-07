import numpy
import pyfits
import intrinsic_util
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot
from scipy.stats.mstats import mquantiles

"This will take the raw catalog of vac5 and all CIV absorbers and make relavent catalogs and plots about data"

def make_catalogs(vac5_file,absorber_file,spAll_file,out_dir):
    dtype={'names':('pk','thingid','plate-mjd-fiber','RA','DEC','type','Z','umag','gmag','rmag','imag','zmag','BAL','DLA','AI','BI','dr7'),
                'formats':(numpy.int64,numpy.uint64,'S20',numpy.float64,numpy.float64,'S10',numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.int64,numpy.int64,numpy.int64,numpy.int64,'S20')}
    vac5_data=numpy.loadtxt(vac5_file,skiprows=1,dtype=dtype)

    #Remove BAL and other things
    vac5_data=vac5_data[vac5_data['type']=='QSO']
    vac5_dict=dict(zip(vac5_data['plate-mjd-fiber'],vac5_data))

    # get info from spall
    spall=pyfits.open(spAll_file)[1].data
    pmf_spall=['%s-%s-%s'%(record['PLATE'],record['MJD'],record['FIBERID']) for record in spall]
    spall_dict=dict(zip(pmf_spall,spall))
    sn_median_r_vac=numpy.array([spall_dict[pmf_i]['SN_MEDIAN'][2] for pmf_i in vac5_data['plate-mjd-fiber']])

    #save the vac catalog
    newdtype={'names':('pk','thingid','plate-mjd-fiber','RA','DEC','type','Z','umag','gmag','rmag','imag','zmag','BAL','DLA','AI','BI','dr7','SN_MEDIAN_R'),
              'formats':(numpy.int64,numpy.uint64,'S20',numpy.float64,numpy.float64,'S10',numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.int64,numpy.int64,numpy.int64,numpy.int64,'S20',numpy.float64)}
    new_vac=numpy.empty(len(vac5_data),dtype=newdtype)
    for i in range(17):
        new_vac[newdtype['names'][i]]=vac5_data[newdtype['names'][i]]
    new_vac['SN_MEDIAN_R']=sn_median_r_vac
    numpy.save('%s/qso_cat.npy'%out_dir,new_vac)

    #Get the absorbers
    absorbers=numpy.load(absorber_file)
    pmf_absorber=['%s-%s-%s'%(record['PLATE'],record['MJD'],record['FIBER']) for record in absorbers]
    print 'Size of absorbers',len(absorbers)
    cond=numpy.array([key in vac5_dict for key in pmf_absorber])
    absorbers=absorbers[cond] #remove the absorbers which are not there in vac catalog
    pmf_absorber=['%s-%s-%s'%(record['PLATE'],record['MJD'],record['FIBER']) for record in absorbers]
    print 'Size of absorbers from vac',len(absorbers)
    sn_median_r=numpy.array([spall_dict[pmf_i]['SN_MEDIAN'][2] for pmf_i in pmf_absorber])

    #Save absorbers with extra info
    newdtype={'names':('RA', 'DEC', 'Z', 'ZQSO', 'R_EW', 'R_EWERR', 'IMAG', 'SN', 'R_FWHM', 'TARGET', 'BETA', 'PLATE', 'MJD', 'FIBER', 'CHUNK','SN_MEDIAN_R','VAC_Z','RMAG') , 
            'formats':(numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.float64,numpy.int64,numpy.float64,numpy.int32,numpy.int32,numpy.int32,'S10',numpy.float64,numpy.float64,numpy.float64)}
    new_absorbers=numpy.empty(len(absorbers),dtype=newdtype)
    for i in range(15):
        new_absorbers[newdtype['names'][i]]=absorbers[newdtype['names'][i]]
    new_absorbers['SN_MEDIAN_R']=sn_median_r
    new_absorbers['VAC_Z']=numpy.array([vac5_dict[pmf_i]['Z'] for pmf_i in pmf_absorber])
    new_absorbers['RMAG']=numpy.array([vac5_dict[pmf_i]['rmag'] for pmf_i in pmf_absorber])
    numpy.save('%s/absorbers_info.npy'%out_dir,new_absorbers)

    #create final cut
    beta_min=-.01
    zab_high,zab_low=intrinsic_util.zab_limits(new_absorbers['VAC_Z'],beta=beta_min)
    cond=numpy.logical_and(new_absorbers['Z']<zab_high,new_absorbers['Z']>zab_low)
    cond=numpy.logical_and(cond, numpy.logical_or(new_absorbers['Z']<intrinsic_util.remove_redshift_range1[0],new_absorbers['Z']>intrinsic_util.remove_redshift_range1[1]))
    cond=numpy.logical_and(cond, numpy.logical_or(new_absorbers['Z']<intrinsic_util.remove_redshift_range2[0],new_absorbers['Z']>intrinsic_util.remove_redshift_range2[1]))
    numpy.save('%s/abs_final.npy'%out_dir,new_absorbers[cond])

    # one with snr>7
    quasar_cond=new_vac['SN_MEDIAN_R']>6.0
    abs_cond=new_absorbers['SN_MEDIAN_R']>6.0
    abs_cond=numpy.logical_and(abs_cond,new_absorbers['R_EW']>.5)
    abs_cond=numpy.logical_and(abs_cond,new_absorbers['R_EW']<3.)
    abs_cond=numpy.logical_and(abs_cond,new_absorbers['R_FWHM'] <600.)
    abs_cond=numpy.logical_and(abs_cond,cond)
    numpy.save('%s/qso_final1.npy'%out_dir,new_vac[quasar_cond])
    numpy.save('%s/abs_final1.npy'%out_dir,new_absorbers[abs_cond])

    #      furthe devide the sample in low luminosity and high luminosity
    absolute_mag_abs=new_absorbers['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(new_absorbers['VAC_Z'])*(1+new_absorbers['VAC_Z']))
    absolute_mag_qso=new_vac['rmag']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(new_vac['Z'])*(1+new_vac['Z']))
    median_mag_abs=numpy.median(absolute_mag_abs[abs_cond])
    numpy.save('%s/qso_final1_lowlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso > median_mag_abs)])
    numpy.save('%s/abs_final1_lowlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs > median_mag_abs)])

    numpy.save('%s/qso_final1_highlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso <= median_mag_abs)])
    numpy.save('%s/abs_final1_highlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs <= median_mag_abs)])

    tritile_points = mquantiles(absolute_mag_abs[abs_cond],prob=(1./3,2./3))
    print "tritile_points",tritile_points
    numpy.save('%s/qso_final1_tri_lowlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso > tritile_points[1])])
    numpy.save('%s/abs_final1_tri_lowlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs > tritile_points[1])])
   
    numpy.save('%s/qso_final1_tri_highlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso <= tritile_points[0])]) 
    numpy.save('%s/abs_final1_tri_highlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs <= tritile_points[0])])

    numpy.save('%s/qso_final1_tri_medlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,\
                                                numpy.logical_and(absolute_mag_qso > tritile_points[0],absolute_mag_qso <= tritile_points[1]))])
    numpy.save('%s/abs_final1_tri_medlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,\
                                                numpy.logical_and(absolute_mag_abs > tritile_points[0],absolute_mag_abs <= tritile_points[1]))])



    # one with snr>10
    quasar_cond=new_vac['SN_MEDIAN_R']>7.0
    abs_cond=new_absorbers['SN_MEDIAN_R']>7.0
    abs_cond=numpy.logical_and(abs_cond,new_absorbers['R_EW']>.45)
    abs_cond=numpy.logical_and(abs_cond,new_absorbers['R_EW']<3.)
    abs_cond=numpy.logical_and(abs_cond,new_absorbers['R_FWHM'] <600.)
    abs_cond=numpy.logical_and(abs_cond,cond)
    numpy.save('%s/qso_final2.npy'%out_dir,new_vac[quasar_cond])
    numpy.save('%s/abs_final2.npy'%out_dir,new_absorbers[abs_cond])

    #      furthe devide the sample in low luminosity and high luminosity
    absolute_mag_abs=new_absorbers['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(new_absorbers['VAC_Z'])*(1+new_absorbers['VAC_Z']))
    absolute_mag_qso=new_vac['rmag']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(new_vac['Z'])*(1+new_vac['Z']))
    median_mag_abs=numpy.median(absolute_mag_abs[abs_cond])
    numpy.save('%s/qso_final2_lowlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso > median_mag_abs)])
    numpy.save('%s/abs_final2_lowlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs > median_mag_abs)])

    numpy.save('%s/qso_final2_highlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso <= median_mag_abs)])
    numpy.save('%s/abs_final2_highlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs <= median_mag_abs)])

    tritile_points = mquantiles(absolute_mag_abs[abs_cond],prob=(1./3,2./3))
    print "tritile_points",tritile_points
    numpy.save('%s/qso_final2_tri_lowlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso > tritile_points[1])])
    numpy.save('%s/abs_final2_tri_lowlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs > tritile_points[1])])

    numpy.save('%s/qso_final2_tri_highlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,absolute_mag_qso <= tritile_points[0])])
    numpy.save('%s/abs_final2_tri_highlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,absolute_mag_abs <= tritile_points[0])])

    numpy.save('%s/qso_final2_tri_medlumin.npy'%out_dir,new_vac[numpy.logical_and(quasar_cond,\
                                                numpy.logical_and(absolute_mag_qso > tritile_points[0],absolute_mag_qso <= tritile_points[1]))])
    numpy.save('%s/abs_final2_tri_medlumin.npy'%out_dir,new_absorbers[numpy.logical_and(abs_cond,\
                                                numpy.logical_and(absolute_mag_abs > tritile_points[0],absolute_mag_abs <= tritile_points[1]))])




def make_cut_table(abs_file):
    abs_data=numpy.load(abs_file)
    print 'No cut &',numpy.size(abs_data)
    cond =abs_data['R_EW']< 3.
    print '\$EW~<~3\$~\\AA &', numpy.sum(cond)
    cond =numpy.logical_and(cond , abs_data['R_EW']> .5)
    print '\$EW~>~0.5\$~\\AA &', numpy.sum(cond)
    cond =numpy.logical_and(cond , abs_data['SN_MEDIAN_R'] >6)
    print 'SN\$_r~>~6\$ &',numpy.sum(cond)
    cond =numpy.logical_and(cond , abs_data['R_FWHM'] < 600)
    print 'FWHM\$~<~\$ 600 km s\$^{-1}\$ &',numpy.sum(cond)
    zab_high,zab_low=intrinsic_util.zab_limits(abs_data['VAC_Z'],beta=-.01)
    temp_cond=numpy.logical_and(abs_data['Z']<zab_high,abs_data['Z']>zab_low)
    cond =numpy.logical_and(cond ,temp_cond)
    print 'Redshift range cut',numpy.sum(cond)
    cond=numpy.logical_and(cond, numpy.logical_or(abs_data['Z']<intrinsic_util.remove_redshift_range1[0],abs_data['Z']>intrinsic_util.remove_redshift_range1[1]))
    print '\$\lambda_o~<~\$5550 \\AA\\ ,\$\\lambda_o~>~\$5620 \\AA\\ &', numpy.sum(cond)
    cond=numpy.logical_and(cond, numpy.logical_or(abs_data['Z']<intrinsic_util.remove_redshift_range2[0],abs_data['Z']>intrinsic_util.remove_redshift_range2[1]))
    print '\$\\lambda_o~<~\$5850 \\AA\\ ,\$\\lambda_o~>~$5950 \\AA\\  &', numpy.sum(cond)


def make_vaious_plot(directory,absorber_file,qso_file):
    pp = PdfPages('intrinsic_absorbers.pdf')
    qso_data=numpy.load('%s/%s'%(directory,qso_file))
    absorbers_data=numpy.load('%s/%s'%(directory,absorber_file))
    
    #quasar histogram
    pyplot.figure()
    pyplot.hist(qso_data['Z'],bins=20,histtype='step',label='%i'%len(qso_data))
    pyplot.xlabel('z')
    pyplot.ylabel('#')
    pyplot.legend()
    pp.savefig()

    #absorber redshift histogram
    pyplot.figure()
    pyplot.hist(absorbers_data['Z'],bins=20,histtype='step',label='%i'%len(absorbers_data))
    pyplot.xlabel('z')
    pyplot.ylabel('#')
    pyplot.legend()
    pp.savefig()

    #absorber SN_r vs Wr
    hist_snr_wr,xedges,yedges=numpy.histogram2d(absorbers_data['R_EW'],absorbers_data['SN_MEDIAN_R'],bins=(numpy.linspace(0,5,31),numpy.linspace(0,30,31)))
    pyplot.figure()
    pyplot.imshow(numpy.transpose(hist_snr_wr),interpolation='nearest',extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],aspect='auto')
    pyplot.xlabel(r'$W_r\; \AA$')
    pyplot.ylabel(r'SN_r')
    pyplot.colorbar()
    pp.savefig()

    #absorber SN_r vs Wr
    tempcond= absorbers_data['SN_MEDIAN_R'] >6
    data_ew = absorbers_data['R_EW'][tempcond]
    data_sn =absorbers_data['SN_MEDIAN_R'][tempcond]
    ew_bin=numpy.arange(0.,3.1,.1)
    ew_mid= .5*ew_bin[1:]+.5*ew_bin[:-1]
    sort_indices=numpy.argsort(data_ew)
    indices = numpy.searchsorted(data_ew[sort_indices],ew_bin)
    median_bin =numpy.array( map(numpy.median,numpy.array_split(data_sn[sort_indices] ,indices[1:])))
    mean_bin =numpy.array( map(numpy.mean,numpy.array_split(data_sn[sort_indices] ,indices[1:])))
    pyplot.figure()
    pyplot.semilogy(absorbers_data['R_EW'],absorbers_data['SN_MEDIAN_R'],'.',markersize=2)
    pyplot.semilogy(ew_mid,median_bin[:-1],'-k')
    pyplot.semilogy(ew_mid,mean_bin[:-1],'--r')
    pyplot.xlim((0,3))
    pyplot.ylim((1,100))
    pp.savefig()
    

    #absorber significance vs Wr
    hist_sn_wr,xedges,yedges=numpy.histogram2d(absorbers_data['R_EW'],absorbers_data['SN'],bins=(numpy.linspace(0,5,31),numpy.linspace(0,10,31)))
    pyplot.figure()
    pyplot.imshow(numpy.transpose(hist_sn_wr),interpolation='nearest',extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],aspect='auto')
    pyplot.xlabel(r'$W_r\; \AA$')
    pyplot.ylabel(r'Significance')
    pyplot.colorbar()
    pp.savefig()

    #absorber SN_r vs fwhm 
    hist_snr_fwhm,xedges,yedges=numpy.histogram2d(absorbers_data['R_FWHM'],absorbers_data['SN_MEDIAN_R'],bins=(numpy.linspace(0,600,31),numpy.linspace(0,30,31)))
    pyplot.figure()
    pyplot.imshow(numpy.transpose(hist_snr_fwhm),interpolation='nearest',extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],aspect='auto')
    pyplot.xlabel(r'FWHM')
    pyplot.ylabel(r'SN_r')
    pyplot.colorbar()
    pp.savefig()

    for cut in [5.,10.]:
        #absorber SN_r vs Wr for Sn_r > cut
        cond=absorbers_data['SN_MEDIAN_R']>cut
        cond=numpy.logical_and(cond,absorbers_data['R_EW'] >0.3)
        cond=numpy.logical_and(cond,absorbers_data['R_EW'] <3.)
        cond_qso=qso_data['SN_MEDIAN_R']>cut

        pyplot.figure()
        bins=numpy.linspace(1.5,4,41)
        pyplot.hist(qso_data['Z'][cond_qso],bins=bins,histtype='step',label='qso(%i)'%numpy.sum(cond_qso))
        pyplot.hist(absorbers_data['Z'][cond],bins=bins,histtype='step',label='absorber(%i)'%numpy.sum(cond))
        pyplot.title('SN_MEDIAN_R> %f'%cut)
        pyplot.xlabel('z')
        pyplot.ylabel('#')
        pyplot.legend()
        pp.savefig()
    
        hist_snr_wr,xedges,yedges=numpy.histogram2d(absorbers_data['R_EW'][cond],absorbers_data['SN_MEDIAN_R'][cond],bins=(numpy.linspace(.3,3,31),numpy.linspace(cut,20,31)))
        pyplot.figure()
        pyplot.imshow(numpy.transpose(hist_snr_wr),interpolation='nearest',extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],aspect='auto')
        pyplot.title('SN_MEDIAN_R> %f'%cut)
        pyplot.xlabel(r'$W_r\; \AA$')
        pyplot.ylabel(r'SN_r')
        pyplot.colorbar()
        pp.savefig()

        pyplot.figure()
        pyplot.hist(absorbers_data['R_EW'][cond],bins=50,label='absorber(%i)'%numpy.sum(cond),histtype='step')
        pyplot.title('SN_MEDIAN_R> %f'%cut)
        pyplot.xlabel(r'$W_r\; \AA$')
        pyplot.ylabel('#')
        pyplot.legend()
        pp.savefig()

        pyplot.figure()
        bins=numpy.linspace(2,15,51)
        pyplot.hist(absorbers_data['SN'][cond],bins=bins,label='absorber(%i)'%numpy.sum(cond),normed=True,histtype='step')
        pyplot.hist(absorbers_data['SN'],bins=bins,label='absorber(%i)'%len(cond),normed=True,histtype='step')
        pyplot.title('SN_MEDIAN_R> %f'%cut)
        pyplot.xlabel(r'SN')
        pyplot.ylabel('#')
        pyplot.legend()
        pp.savefig()

    pyplot.figure()
    bins=numpy.linspace(1.5,4.,51)
    for s in ['1','2','3']:
        abs_data=numpy.load('%s/abs_final%s.npy'%(directory,s))
        pyplot.hist(abs_data['Z'],bins=bins,histtype='step',label='dataset%s(%i)'%(s,len(abs_data)))
    pyplot.ylabel('#')
    pyplot.xlabel('z')
    pyplot.legend()
    pp.savefig()

    pyplot.figure()
    bins=numpy.linspace(1.5,4.,51)
    for s in ['1','2']:
        qso_data=numpy.load('%s/qso_final%s.npy'%(directory,s))
        pyplot.hist(qso_data['Z'],bins=bins,histtype='step',label='dataset%s(%i)'%(s,len(qso_data)))
    pyplot.ylabel('#')
    pyplot.xlabel('z')
    pyplot.legend()
    pp.savefig()

    pyplot.figure()
    bins=numpy.linspace(.3,3.0,28)
    bin_mids=(bins[1:]+bins[:-1])/2.

    abs_data=numpy.load('%s/abs_final.npy'%directory)
    abs_cond1=abs_data['SN_MEDIAN_R']>10.0
    #abs_cond1=numpy.logical_and(abs_cond1,abs_data['SN']>5.)
    abs_cond1=numpy.logical_and(abs_cond1,abs_data['R_EW']>.3)
    abs_cond1=numpy.logical_and(abs_cond1,abs_data['R_EW']<3.)
    abs_cond1=numpy.logical_and(abs_cond1,abs_data['R_FWHM'] <600.)

    hist1,junk = numpy.histogram(abs_data['R_EW'][abs_cond1],bins=bins)
    pyplot.errorbar(bin_mids,hist1,yerr=numpy.sqrt(hist1),label='10')

    for sn_r_i in [5.,6.,7.,8.,9.]:
        abs_cond2=abs_data['SN_MEDIAN_R']>sn_r_i
        #abs_cond2=numpy.logical_and(abs_cond2,abs_data['SN']>5.)
        abs_cond2=numpy.logical_and(abs_cond2,abs_data['R_EW']>.3)
        abs_cond2=numpy.logical_and(abs_cond2,abs_data['R_EW']<3.)
        abs_cond2=numpy.logical_and(abs_cond2,abs_data['R_FWHM'] <600.)

        hist2,junk = numpy.histogram(abs_data['R_EW'][abs_cond2],bins=bins)
        ratio=numpy.float(numpy.sum(hist1[bin_mids>1.]))/numpy.sum(hist2[bin_mids>1.])
        pyplot.errorbar(bin_mids+.01,ratio*hist2,yerr=numpy.sqrt(hist2)*ratio,label='%i(%i)'%(int(sn_r_i),numpy.sum(abs_cond2)))
    pyplot.ylabel('#')
    pyplot.xlabel('$W_r\; \AA$')
    pyplot.legend()
    pp.savefig()

    pyplot.figure()
    bins=numpy.linspace(.3,3.0,28)
    bin_mids=(bins[1:]+bins[:-1])/2.

    abs_data=numpy.load('%s/abs_final.npy'%directory)
    abs_cond1=abs_data['SN_MEDIAN_R']>10.0
    abs_cond1=numpy.logical_and(abs_cond1,abs_data['R_EW']>.3)
    abs_cond1=numpy.logical_and(abs_cond1,abs_data['R_EW']<3.)
    abs_cond1=numpy.logical_and(abs_cond1,abs_data['R_FWHM'] <600.)

    hist1,junk = numpy.histogram(abs_data['R_EW'][abs_cond1],bins=bins)
    pyplot.errorbar(bin_mids,hist1,yerr=numpy.sqrt(hist1),label='10(%i)'%numpy.sum(abs_cond1))

    for sn_i in [5.,6.,7.,8.,9.]:
        abs_cond2=abs_data['SN']>sn_i
        abs_cond2=numpy.logical_and(abs_cond2,abs_data['R_EW']>.3)
        abs_cond2=numpy.logical_and(abs_cond2,abs_data['R_EW']<3.)
        abs_cond2=numpy.logical_and(abs_cond2,abs_data['R_FWHM'] <600.)

        hist2,junk = numpy.histogram(abs_data['R_EW'][abs_cond2],bins=bins)
        ratio=numpy.float(numpy.sum(hist1[bin_mids>1.]))/numpy.sum(hist2[bin_mids>1.])
        pyplot.errorbar(bin_mids+.01,ratio*hist2,yerr=numpy.sqrt(hist2)*ratio,label='%i(%i)'%(int(sn_i),numpy.sum(abs_cond2)))
    pyplot.ylabel('#')
    pyplot.xlabel('$W_r\; \AA$')
    pyplot.legend()
    pp.savefig()


    pyplot.figure()
    abs_data=numpy.load('%s/abs_final2.npy'%directory)
    abs_data_beta = (abs_data['VAC_Z'] - abs_data['Z']) * (2 + abs_data['VAC_Z'] + abs_data['Z'])/\
                        ((1+abs_data['VAC_Z'])**2 + (1+abs_data['Z'])**2)
    pyplot.plot(abs_data_beta,abs_data['R_FWHM'],'.')
    pyplot.xlabel(r'$\beta$')
    pyplot.ylabel('FWHM')
    pp.savefig()
 

    pp.close()

if __name__=='__main__':
    if False:
        make_catalogs('/n/panlfs/vikas/boss_data/v_5_4_45/vac_catalog/vac5/fpg_vac5.txt',                            #vac5 file
                    '/n/panlfs/vikas/boss_data/v_5_4_45/absorbers/preliminary_catalogs/boss_all_no_cuts.npy',    #All CIV absorbers
                    '/n/panlfs/vikas/boss_data/groups/boss/spectro/redux/spAll-v5_4_45.fits',                    # spall file
                    '/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs')                                # output dir
    if False:
        make_vaious_plot('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs',
                        'absorbers_info.npy',
                        'qso_cat.npy')
    if True:
        make_cut_table('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/absorbers_info.npy')
    
