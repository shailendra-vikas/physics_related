import numpy
import intrinsic_util
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter
from scipy.integrate import trapz
size_val=25
format='eps'

def probablity_vs_r():
    rest_frame_r=numpy.linspace(-100,800,10001)
    zqso=2.4
    epsilon=-1.45
    sigma_beta=0.003
    dp_dr_1h_values=intrinsic_util.dp_dr_1h(rest_frame_r,zqso,epsilon=epsilon,sigma_beta=sigma_beta)
    dp_dr_2h_values=intrinsic_util.dp_dr_2h(rest_frame_r,zqso,epsilon=epsilon)
    dp_dr_total_values=dp_dr_2h_values+5000*dp_dr_1h_values
    sigma_z_in_r=2.0# error 1 mpc
    weight=intrinsic_util.get_gaussian_window(sigma_z_in_r,rest_frame_r[1]-rest_frame_r[0])
    pyplot.figure()
    pyplot.plot(rest_frame_r,5000*dp_dr_1h_values,'-b',label='1h')
    pyplot.plot(rest_frame_r,dp_dr_2h_values,'-g',label='2h')
    pyplot.plot(rest_frame_r,dp_dr_2h_values+5000*dp_dr_1h_values,'-k',label='total')
    pyplot.plot(rest_frame_r,numpy.convolve(dp_dr_total_values,weight,mode='same'),'-r',label=r'convolved with error=2 h$^{-1}$ Mpc ')
    pyplot.xlabel(r'r (h$^{-1}$ Mpc)')
    pyplot.ylabel(r'$\propto dP/dr$')
    pyplot.xlim((-50,200))
    pyplot.ylim((0,50))
    pyplot.legend()
    pyplot.savefig('dp_dr.%s'%format,fmt=format)

def various_plots(file_name):
    labels,values,cov,total_1halo,total_2halo,total_outflow,error_total_1halo,error_total_2halo,error_total_outflow,selections,data_bins,civ_hist,civ_hist_err =\
        numpy.load(file_name)
    param_errors=numpy.sqrt([cov[i,i] for i  in range(len(labels))])
    data_bin_mids=0.5*data_bins[1:]+ 0.5*data_bins[:-1]
    beta_range=numpy.linspace(-.01,.2,1001)




    pyplot.figure()
    pyplot.errorbar(data_bin_mids,civ_hist,civ_hist_err,fmt='o',color='b',label='data')
    pyplot.fill_between(beta_range,total_1halo+total_2halo+total_outflow+error_total_1halo+error_total_2halo+error_total_outflow,\
                        total_1halo+total_2halo+total_outflow-error_total_1halo-error_total_2halo-error_total_outflow,color=(.5,1.,.5))
    pyplot.plot(beta_range,total_1halo+total_2halo+total_outflow,color='k',label='Theory')
    pyplot.xlabel(r'$\beta = v/c$')
    pyplot.ylabel(r'$dN/d\beta$')
    pyplot.legend()
    pyplot.savefig('therory_vs_data.%s'%format,fmt=format)

    pyplot.figure()
    pyplot.fill_between(beta_range, total_1halo + error_total_1halo, total_1halo - error_total_1halo, color=(1.,.5,0.5))
    pyplot.fill_between(beta_range, total_2halo + error_total_2halo, total_2halo - error_total_2halo, color=(0.5,1.,0.5))
    pyplot.fill_between(beta_range, total_outflow + error_total_outflow, total_outflow - error_total_outflow, color=(0.5,0.5,1.))
    pyplot.plot(beta_range, total_1halo, color='k', ls='solid', label='One halo')
    pyplot.plot(beta_range, total_2halo, color='k', ls='dashed', label='Two halo')
    pyplot.plot(beta_range, total_outflow, color='k', ls='dotted', label='Outflow')
    pyplot.xlabel(r'$\beta = v/c$')
    pyplot.ylabel(r'$dN/d\beta$')
    ylim_min,ylim_max=pyplot.ylim()
    pyplot.ylim((0,ylim_max))
    pyplot.legend()
    pyplot.savefig('components.%s'%format,fmt=format)

    pyplot.figure()
    ax1 = pyplot.subplot(1,1,1)
    ax1.fill_between(beta_range, total_1halo + error_total_1halo, total_1halo - error_total_1halo, color=(1.,.5,0.5))
    ax1.fill_between(beta_range, total_2halo + error_total_2halo, total_2halo - error_total_2halo, color=(0.5,1.,0.5))
    ax1.fill_between(beta_range, total_outflow + error_total_outflow, total_outflow - error_total_outflow, color=(0.5,0.5,1.))
    ax1.plot(beta_range, total_1halo, color='k', ls='solid', label='One halo')
    ax1.plot(beta_range, total_2halo, color='k', ls='dashed', label='Two halo')
    ax1.plot(beta_range, total_outflow, color='k', ls='dotted', label='Outflow')
    ax1.set_xlabel(r'$\beta = v/c$')
    ax1.set_ylabel(r'$dN/d\beta$')
    ylim_min,ylim_max=ax1.get_ylim()
    ax1.set_ylim((0,ylim_max))
    ax1.set_xlim((-.01,.2))
    ax2=pyplot.axes((.3,.3,.6,.6))
    ax2.errorbar(data_bin_mids,civ_hist,civ_hist_err,fmt='o',color='b',label='data')
    ax2.fill_between(beta_range,total_1halo+total_2halo+total_outflow+error_total_1halo+error_total_2halo+error_total_outflow,\
                        total_1halo+total_2halo+total_outflow-error_total_1halo-error_total_2halo-error_total_outflow,color=(.5,1.,.5))
    ax2.plot(beta_range,total_1halo+total_2halo+total_outflow,color='k',label='Theory')
    pyplot.setp(ax2, xlim=(-.01,.2), xticks=[], yticks=[])
    ax1.legend()
    pyplot.savefig('two_in_one.%s'%format,fmt=format,transparent=True)



    latex_labels={'n0_sigma0':r'$n_0 \sigma_0$', 'epsilon':r'$\epsilon$', 'sigma_beta':r'$\sigma_\beta$', 'f':r'$f$', 'beta_01':r'$n_{2k}$', 'beta_02':r'$n_{4k}$', 'beta_03':r'$n_{6k}$'}
    TEMP_ERR1,TEMP_ERR2= numpy.meshgrid(param_errors,param_errors)
    normed_cov= cov/(TEMP_ERR1*TEMP_ERR2)
    pyplot.figure()
    pyplot.imshow(normed_cov,interpolation='nearest')
    pyplot.xticks(numpy.arange(7),[latex_labels[lab] for lab in labels ],fontsize=20)
    pyplot.yticks(numpy.arange(7),[latex_labels[lab] for lab in labels ],fontsize=20)
    pyplot.colorbar(shrink=0.78)
    pyplot.savefig('cov.%s'%format,fmt=format)

    pyplot.figure()
    no_of_quasars=numpy.sum(selections,axis=1)
    pyplot.plot(beta_range,no_of_quasars/1000.,color='b',lw=2)
    pyplot.xlabel(r'$\beta = v/c$')
    pyplot.ylabel(r'$N_\mathrm{qso}$ (in thousands)')
    pyplot.savefig('no_quasars.%s'%format,fmt=format)

    no_of_quasars_at_data = numpy.interp(data_bin_mids,beta_range,no_of_quasars)
    pyplot.figure()
    tempcond=no_of_quasars_at_data>.5
    pyplot.errorbar(data_bin_mids[tempcond],civ_hist[tempcond]/no_of_quasars_at_data[tempcond], civ_hist_err[tempcond]/no_of_quasars_at_data[tempcond],fmt='o-',color='b')
    pyplot.ylabel(r'$dN/d\beta$')
    pyplot.xlabel(r'$\beta = v/c$')
    pyplot.savefig('only_data.%s'%format,fmt=format)
    
    
def make_luminosity_2bin_plots(abs_file_low,abs_file_high,info_low,info_high):
    abs_low, abs_high = map(numpy.load, (abs_file_low, abs_file_high))
    abs_low_mag = abs_low['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(abs_low['VAC_Z'])*(1+abs_low['VAC_Z']))
    abs_high_mag = abs_high['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(abs_high['VAC_Z'])*(1+abs_high['VAC_Z']))
    #absolute_mag_qso=new_vac['rmag']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(new_vac['Z'])*(1+new_vac['Z']))
    median_mag=numpy.max(abs_high_mag)
    print median_mag
    bins=numpy.linspace(-30,-24,40)
    pyplot.figure()
    pyplot.hist(numpy.append(abs_low_mag,abs_high_mag),bins=bins,histtype='step',lw=2)
    pyplot.axvline(median_mag,ymin=0,ymax=.95,ls='--',color='r',lw=2)
    pyplot.text(-29.5,600,'Median=%1.2f'%median_mag,fontsize=20,color='r')
    pyplot.xlabel('Absolute Magnitude')
    pyplot.ylabel('# per 0.2 mag')
    pyplot.savefig('abs_mag.%s'%format,fmt=format)

    beta_range=numpy.linspace(-.01,.2,1001)
    labels_low,values_low,cov_low,total_1halo_low,total_2halo_low,total_outflow_low,error_total_1halo_low,error_total_2halo_low,error_total_outflow_low,\
            selections_low,data_bins_low,civ_hist_low,civ_hist_err_low = numpy.load(info_low)
    labels_high,values_high,cov_high,total_1halo_high,total_2halo_high,total_outflow_high,error_total_1halo_high,error_total_2halo_high,error_total_outflow_high,\
            selections_high,data_bins_high,civ_hist_high,civ_hist_err_high = numpy.load(info_high)
    
    print numpy.sum(beta_range*total_outflow_high)/numpy.sum(total_outflow_high)
    print numpy.sum(beta_range*total_outflow_low)/numpy.sum(total_outflow_low)
    pyplot.figure()
    pyplot.fill_between(beta_range,total_outflow_high+error_total_outflow_high,total_outflow_high-error_total_outflow_high,color=(.7,.7,1.0))
    pyplot.plot(beta_range,total_outflow_high,color='b',label=r'$M_r < %1.2f$'%median_mag)
    pyplot.fill_between(beta_range,total_outflow_low+error_total_outflow_low,total_outflow_low-error_total_outflow_low,color=(1.0,.7,.7))
    pyplot.plot(beta_range,total_outflow_low,color='r',label=r'$M_r > %1.2f$'%median_mag)
    pyplot.xlabel(r'$\beta = v/c$')
    pyplot.ylabel(r'$dN/d\beta$')
    pyplot.xlim((-.01,.05))
    pyplot.legend()
    pyplot.savefig('outflow.%s'%format,fmt=format)

def make_luminosity_3bin_plots(abs_file_low,abs_file_med,abs_file_high,info_low,info_med,info_high, info_all):
    abs_low, abs_med, abs_high = map(numpy.load, (abs_file_low, abs_file_med, abs_file_high))
    abs_low_mag = abs_low['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(abs_low['VAC_Z'])*(1+abs_low['VAC_Z']))
    abs_med_mag = abs_med['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(abs_med['VAC_Z'])*(1+abs_med['VAC_Z']))
    abs_high_mag = abs_high['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(abs_high['VAC_Z'])*(1+abs_high['VAC_Z']))
    first_tritile= numpy.max(abs_med_mag)
    second_tritile = numpy.max(abs_high_mag)

    bins=numpy.linspace(-30,-24,40)
    pyplot.figure()
    pyplot.hist(numpy.concatenate((abs_low_mag,abs_med_mag,abs_high_mag)),bins=bins,histtype='step',lw=2)
    pyplot.axvline(second_tritile,ymin=0,ymax=.95,ls='--',color='r',lw=2)
    pyplot.axvline(first_tritile,ymin=0,ymax=.95,ls='--',color='r',lw=2)
    pyplot.text(-29.5,600,'$M_r=%1.2f$'%first_tritile,fontsize=20,color='r')
    pyplot.text(-26,600,'$M_r=%1.2f$'%second_tritile,fontsize=20,color='r')
    pyplot.xlabel('Absolute Magnitude')
    pyplot.ylabel('# per 0.2 mag')
    pyplot.savefig('abs_mag_3bin.%s'%format,fmt=format)

    beta_range=numpy.linspace(-.01,.2,1001)
    labels_low,values_low,cov_low,total_1halo_low,total_2halo_low,total_outflow_low,error_total_1halo_low,error_total_2halo_low,error_total_outflow_low,\
            selections_low,data_bins_low,civ_hist_low,civ_hist_err_low = numpy.load(info_low)
    labels_med,values_med,cov_med,total_1halo_med,total_2halo_med,total_outflow_med,error_total_1halo_med,error_total_2halo_med,error_total_outflow_med,\
            selections_med,data_bins_med,civ_hist_med,civ_hist_err_med = numpy.load(info_med)
    labels_high,values_high,cov_high,total_1halo_high,total_2halo_high,total_outflow_high,error_total_1halo_high,error_total_2halo_high,error_total_outflow_high,\
            selections_high,data_bins_high,civ_hist_high,civ_hist_err_high = numpy.load(info_high)
    labels_all,values_all,cov_all,total_1halo_all,total_2halo_all,total_outflow_all,error_total_1halo_all,error_total_2halo_all,error_total_outflow_all,\
            selections_all,data_bins_all,civ_hist_all,civ_hist_err_all = numpy.load(info_all)

    pyplot.figure()
    pyplot.fill_between(beta_range,total_outflow_high+error_total_outflow_high,total_outflow_high-error_total_outflow_high,color=(.7,.7,1.0))
    pyplot.plot(beta_range,total_outflow_high,color='b',label=r'$M_r < %1.2f$'%first_tritile)
    pyplot.fill_between(beta_range,total_outflow_med+error_total_outflow_med,total_outflow_med-error_total_outflow_med,color=(.7,1.0,.7))
    pyplot.plot(beta_range,total_outflow_med,color='g',label=r'$%1.2f < M_r < %1.2f$'%(first_tritile,second_tritile))
    pyplot.fill_between(beta_range,total_outflow_low+error_total_outflow_low,total_outflow_low-error_total_outflow_low,color=(1.0,.7,.7))
    pyplot.plot(beta_range,total_outflow_low,color='r',label=r'$M_r > %1.2f$'%second_tritile)
    pyplot.xlabel(r'$\beta = v/c$')
    pyplot.ylabel(r'$dN/d\beta$')
    pyplot.xlim((-.01,.05))
    pyplot.legend()
    pyplot.savefig('outflow_3bin.%s'%format,fmt=format)


    count_high = map(lambda x:trapz(x,beta_range),(total_1halo_high,total_2halo_high,total_outflow_high))
    count_med = map(lambda x:trapz(x,beta_range),(total_1halo_med,total_2halo_med,total_outflow_med))
    count_low = map(lambda x:trapz(x,beta_range),(total_1halo_low,total_2halo_low,total_outflow_low))
    count_all = map(lambda x:trapz(x,beta_range),(total_1halo_all,total_2halo_all,total_outflow_all))

    count_high_err = map(lambda x_err: numpy.sqrt(numpy.sum(x_err**2))*(beta_range[1]-beta_range[0]),(error_total_1halo_high,error_total_2halo_high,error_total_outflow_high))
    count_med_err = map(lambda x_err: numpy.sqrt(numpy.sum(x_err**2))*(beta_range[1]-beta_range[0]),(error_total_1halo_med,error_total_2halo_med,error_total_outflow_med))
    count_low_err = map(lambda x_err: numpy.sqrt(numpy.sum(x_err**2))*(beta_range[1]-beta_range[0]),(error_total_1halo_low,error_total_2halo_low,error_total_outflow_low))
    count_all_err = map(lambda x_err: numpy.sqrt(numpy.sum(x_err**2))*(beta_range[1]-beta_range[0]),(error_total_1halo_all,error_total_2halo_all,error_total_outflow_all))

    count_total_high = count_high[0] + count_high[1] + count_high[2] 
    count_total_med = count_med[0] + count_med[1] + count_med[2] 
    count_total_low = count_low[0] + count_low[1] + count_low[2] 
    count_total_all = count_all[0] + count_all[1] + count_all[2] 

    print 'Individual component as % of total ' 
    print 'High brightness: one_halo:', 100.*count_high[0]/count_total_high,' +/- ', 100.*count_high_err[0]/count_total_high
    print 'High brightness: two_halo:', 100.*count_high[1]/count_total_high,' +/- ', 100.*count_high_err[1]/count_total_high
    print 'High brightness: outflow:', 100.*count_high[2]/count_total_high,' +/- ', 100.*count_high_err[2]/count_total_high

    print 'Med brightness: one_halo:', 100.*count_med[0]/count_total_med,' +/- ', 100.*count_med_err[0]/count_total_med
    print 'Med brightness: two_halo:', 100.*count_med[1]/count_total_med,' +/- ', 100.*count_med_err[1]/count_total_med
    print 'Med brightness: outflow:', 100.*count_med[2]/count_total_med,' +/- ', 100.*count_med_err[2]/count_total_med

    print 'Low brightness: one_halo:', 100.*count_low[0]/count_total_low,' +/- ', 100.*count_low_err[0]/count_total_low
    print 'Low brightness: two_halo:', 100.*count_low[1]/count_total_low,' +/- ', 100.*count_low_err[1]/count_total_low
    print 'Low brightness: outflow:', 100.*count_low[2]/count_total_low,' +/- ', 100.*count_low_err[2]/count_total_low

    print 'all brightness: one_halo:', 100.*count_all[0]/count_total_all,' +/- ', 100.*count_all_err[0]/count_total_all
    print 'all brightness: two_halo:', 100.*count_all[1]/count_total_all,' +/- ', 100.*count_all_err[1]/count_total_all
    print 'all brightness: outflow:', 100.*count_all[2]/count_total_all,' +/- ', 100.*count_all_err[2]/count_total_all

    peak_high_index = numpy.argmax(total_outflow_high)
    peak_med_index = numpy.argmax(total_outflow_med)
    peak_low_index = numpy.argmax(total_outflow_low)
    peak_all_index = numpy.argmax(total_outflow_all)

    print 'Peak fraction and location'
    print 'High brightness: Peak location:', beta_range[peak_high_index], ' fraction: ', total_outflow_high[peak_high_index]/(total_1halo_high[peak_high_index]+total_2halo_high[peak_high_index]+total_outflow_high[peak_high_index]), ' err:', error_total_outflow_high[peak_high_index]/(total_1halo_high[peak_high_index]+total_2halo_high[peak_high_index]+total_outflow_high[peak_high_index])
    print 'Med brightness: Peak location:', beta_range[peak_med_index], ' fraction: ', total_outflow_med[peak_med_index]/(total_1halo_med[peak_med_index]+total_2halo_med[peak_med_index]+total_outflow_med[peak_med_index]), ' err:', error_total_outflow_med[peak_med_index]/(total_1halo_med[peak_med_index]+total_2halo_med[peak_med_index]+total_outflow_med[peak_med_index])
    print 'Low brightness: Peak location:', beta_range[peak_low_index], ' fraction: ', total_outflow_low[peak_low_index]/(total_1halo_low[peak_low_index]+total_2halo_low[peak_low_index]+total_outflow_low[peak_low_index]), ' err:', error_total_outflow_low[peak_low_index]/(total_1halo_low[peak_low_index]+total_2halo_low[peak_low_index]+total_outflow_low[peak_low_index])
    print 'all brightness: Peak location:', beta_range[peak_all_index], ' fraction: ', total_outflow_all[peak_all_index]/(total_1halo_all[peak_all_index]+total_2halo_all[peak_all_index]+total_outflow_all[peak_all_index]), ' err:', error_total_outflow_all[peak_all_index]/(total_1halo_all[peak_all_index]+total_2halo_all[peak_all_index]+total_outflow_all[peak_all_index])


    beta_high = trapz(total_outflow_high*beta_range,beta_range)/count_high[2]
    beta_high_err = numpy.sqrt(numpy.sum((error_total_outflow_high*beta_range)**2))*(beta_range[1]-beta_range[0])/count_high[2]
    beta_med = trapz(total_outflow_med*beta_range,beta_range)/count_med[2]
    beta_med_err = numpy.sqrt(numpy.sum((error_total_outflow_med*beta_range)**2))*(beta_range[1]-beta_range[0])/count_med[2]
    beta_low = trapz(total_outflow_low*beta_range,beta_range)/count_low[2]
    beta_low_err = numpy.sqrt(numpy.sum((error_total_outflow_low*beta_range)**2))*(beta_range[1]-beta_range[0])/count_low[2]
    beta_all = trapz(total_outflow_all*beta_range,beta_range)/count_all[2]
    beta_all_err = numpy.sqrt(numpy.sum((error_total_outflow_all*beta_range)**2))*(beta_range[1]-beta_range[0])/count_all[2]
    print 'Mean valocity'
    print 'High luminosity : ',beta_high ,' +/- ',beta_high_err
    print 'Med luminosity : ',beta_med,' +/- ', beta_med_err
    print 'Low luminosity:' , beta_low ,' +/- ', beta_low_err
    print 'All luminosity:', beta_all ,' +/- ', beta_all_err
    







    
def make_absorbers_plots(abs_file):
    abs_data = numpy.load(abs_file)
    ew_bins=numpy.linspace(.3,3.,55)
    sn_r_bins =numpy.linspace(7,40,67)
    wr_i=numpy.histogram2d(abs_data['R_EW'],abs_data['SN_MEDIAN_R'],bins=(ew_bins,sn_r_bins))[0]
    left, width = 0.1*.8, 0.65*.8
    bottom, hight = 0.1, 0.65
    up_hight=.2
    side_width=.2*.8

    rect_scatter = [left, bottom, width, hight]
    rect_histx = [left, bottom+hight, width, up_hight]
    rect_histy = [left+width, bottom, side_width, hight]
    rect_bar = [left+width+side_width+side_width/2,bottom,side_width*.1,hight+up_hight]
    nullfmt   = NullFormatter()

    fig=pyplot.figure(figsize=(10,8))
    axScatter = pyplot.axes(rect_scatter)
    axHistx = pyplot.axes(rect_histx)
    axHisty = pyplot.axes(rect_histy)
    axcolorbar=pyplot.axes(rect_bar)

    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    image=axScatter.imshow(numpy.transpose(wr_i),origin='bottom',interpolation='nearest',extent=(.3,3,4,40),aspect='auto')
    fig.colorbar(image,cax=axcolorbar)
    #axScatter.scatter(civ['R_EW'],civ['IMAG'],s=1)
    axScatter.set_xlabel(r'$W_r (\AA)$',size=size_val)
    axScatter.set_ylabel(r'Median SN$_r$',size=size_val)

    axScatter.set_xlim( (.3, 3) )
    axScatter.set_ylim( (7, 40) )

    axHistx.hist(abs_data['R_EW'], bins=ew_bins ,histtype='step',color='black')
    #axHistx.hist(civ['R_EW'][bonus], bins=100 ,histtype='step',color='red',label='bonus')
    axHistx.legend(frameon=False)
    axHisty.hist(abs_data['SN_MEDIAN_R'], bins=sn_r_bins, orientation='horizontal',histtype='step',color='black')
    #axHisty.hist(civ['IMAG'][bonus], bins=100, orientation='horizontal',histtype='step',color='red',label='bonus')
    #axHisty.legend(frameon=False,loc=4)
    axHisty.set_xticks([300,800])

    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    pyplot.savefig('wr_r.%s'%format,fmt=format)
    
def make_completeness_plot(abs_file):
    absorbers_data=numpy.load(abs_file)

    absorbers_data=absorbers_data[absorbers_data['R_FWHM'] <600.]    
    #abs_cond=abs_data['SN_MEDIAN_R']>6.0
    #abs_cond=numpy.logical_and(abs_cond,abs_data['R_EW']>.5)
    #abs_cond=numpy.logical_and(abs_cond,abs_data['R_EW']<3.)

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
    pyplot.axvline(0.5,color='g',lw=2)
    pyplot.axhline(6,color='g',lw=2)
    pyplot.xlabel(r'EW ($\AA$)')
    pyplot.ylabel(r'SN$_r$')
    pyplot.xlim((.3,3))
    pyplot.ylim((1,100))
    pyplot.savefig('completeness.%s'%format,fmt=format)

    tempcond=numpy.logical_and(absorbers_data['R_EW']>.5 , absorbers_data['R_EW']<3.)
    tempcond=numpy.logical_and(tempcond,absorbers_data['SN_MEDIAN_R']>6.)
    abs_mag=absorbers_data['RMAG']-25-5*numpy.log10(intrinsic_util.fast_redshift_to_comoving_instance.distance(absorbers_data['VAC_Z'])*(1+absorbers_data['VAC_Z']))
    abs_data_beta = (absorbers_data['VAC_Z'] - absorbers_data['Z']) * (2 + absorbers_data['VAC_Z'] + absorbers_data['Z'])/\
                            ((1+absorbers_data['VAC_Z'])**2 + (1+absorbers_data['Z'])**2)

    pyplot.figure()
    pyplot.plot(abs_mag[tempcond],abs_data_beta[tempcond],'.',markersize=2)
    pyplot.xlabel('Absolute Magnitude')
    pyplot.ylabel(r'$\beta$')
    pyplot.xlim((-30,-24))
    pyplot.savefig('beta_vs_mag.%s'%format,fmt=format)
    
    


def main():
    #probablity_vs_r()
    #various_plots('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_all.npy')
    #make_luminosity_2bin_plots('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_lowlumin.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_highlumin.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_low.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_high.npy')
    #make_luminosity_3bin_plots('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_tri_lowlumin.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_tri_highlumin.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1_tri_medlumin.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_tri_low.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_tri_med.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_tri_high.npy','/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/fit_info_final1_all.npy')
    #make_absorbers_plots('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final1.npy')
    #make_completeness_plot('/n/panlfs/vikas/boss_data/v_5_4_45/intrinsic_data/catalogs/abs_final.npy')

if __name__=='__main__':
    main()
