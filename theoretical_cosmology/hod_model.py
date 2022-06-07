import NumberDensity
import numpy
from scipy import integrate
import Cosmology
import NFW

def NFW_func(param,mass,m_star,kk):
    nfw=NFW.NFWProfile(param,mass,m_star)
    return nfw.u(kk)

def large_scale_bias(param,Mmin,massfraction='ShethTorman',occupation=None):
    numberdensity=NumberDensity.NumberDensity(param,massfraction=massfraction)
    if occupation is None:
        number_occupation=lambda m:1
    else:
        number_occupation=lambda m:m/Mmin
    def integrand_number(logm):
        return (10**logm) * numberdensity(10**logm)*number_occupation(10**logm)
    N_object=integrate.quad(integrand_number,numpy.log10(Mmin),20)[0]*numpy.log(10)

    def integrand_bias(logm):
        return (10**logm) * numberdensity(10**logm)*number_occupation(10**logm)*numberdensity.halo_bias(10**logm)
    bias_object=integrate.quad(integrand_bias,numpy.log10(Mmin),20)[0]*numpy.log(10)

    return (bias_object/N_object,N_object)


def intergraltest(param,massfraction='ShethTorman'):
    numberdensity=NumberDensity.NumberDensity(param,massfraction=massfraction)
    def integrand_bias(logm):
        return (10**logm)**2 * numberdensity(10**logm)*numberdensity.halo_bias(10**logm)
    bias_object=integrate.quad(integrand_bias,4,20)[0]*numpy.log(10)
    print 'Expected 1 got:',bias_object/(Cosmology.CRITICAL_DENSITY_Msun_Mpc*(param.omg_b+param.omg_c)*(1+param.z)**3)
    return bias_object/(Cosmology.CRITICAL_DENSITY_Msun_Mpc*(param.omg_b+param.omg_c)*(1+param.z)**3)

def get_minmass_for_bias(param,bias,bias_err,massfraction='ShethTorman',occupation=None):
    Mmin=numpy.logspace(10,13,20)
    bias_value,number_value=[],[]
    for Mmin_i in Mmin:
        bias_and_number=large_scale_bias(param,Mmin_i,massfraction='ShethTorman',occupation=None)
        bias_value.append(bias_and_number[0])
    bias_value=numpy.array(bias_value)
    result_Mmin=numpy.interp(bias,bias_value,Mmin)
    result_Mmin_upper=numpy.interp(bias+bias_err,bias_value,Mmin)
    result_Mmin_lower=numpy.interp(bias-bias_err,bias_value,Mmin)
    return result_Mmin,result_Mmin_upper,result_Mmin_lower

def small_scale_bias(param,kk,Mmin,massfraction='ShethTorman',occupation=None):
    numberdensity=NumberDensity.NumberDensity(param,massfraction=massfraction)
    if occupation is None:
        number_occupation=lambda m:1
    else:
        number_occupation=lambda m:m/Mmin

    def integrand_number(logm):
        return (10**logm) * numberdensity(10**logm)*number_occupation(10**logm)
    N_object=integrate.quad(integrand_number,numpy.log10(Mmin),20)[0]*numpy.log(10)

    k_value=None
    def integrand_bias(logm):
        return (10**logm)*numberdensity(10**logm)*number_occupation(10**logm)*numberdensity.halo_bias(10**logm)* NFW_func(param,10**logm,numberdensity.sigmaFunc.getM_star(),k_value)

    bias_object=[]
    for kk_i in kk:
        k_value=kk_i
        bias_object.append(integrate.quad(integrand_bias,numpy.log10(Mmin),20)[0]*numpy.log(10))
    return numpy.array(bias_object/N_object)

def qso_bias_corr():
    param=Cosmology.CosmologyParam(omg_b=0.04,omg_c=0.22,omg_l=0.7,omg_k=0.0,sigma8=0.8,Tcmb=2.728,h=0.72,z=2.4)
    bias=3.70732
    bias_err=.090377
    kk=numpy.logspace(-4,2,151)
    qso_min_mass,qso_min_mass_upper,qso_min_mass_lower=get_minmass_for_bias(param,bias,bias_err)
    print qso_min_mass,qso_min_mass_upper,qso_min_mass_lower
    #qso_min_mass=2.60566e12
    bias_kk=small_scale_bias(param,kk,qso_min_mass)

    numberdensity=NumberDensity.NumberDensity(param,massfraction='ShethTorman')
    delta_qso_kk=numberdensity.sigmaFunc.deltaSq(kk)*bias_kk
    
    rr=numpy.logspace(-2,3,31)
    KK,RR=numpy.meshgrid(kk,rr)
    KKRR=KK*RR
    sinc_KK=numpy.sin(KKRR)/KKRR
    corr=integrate.trapz(delta_qso_kk*sinc_KK,x=kk,axis=1)*numpy.log(10)
    
    from matplotlib import pyplot
    pyplot.figure()
    pyplot.loglog(kk,bias_kk)
    pyplot.savefig('qso_power.png',fmt='png')
    pyplot.figure()
    pyplot.loglog(rr,corr)
    pyplot.savefig('qso_corr.png',fmt='png')


def civ_bias_corr():
    param=Cosmology.CosmologyParam(omg_b=0.04,omg_c=0.22,omg_l=0.7,omg_k=0.0,sigma8=0.8,Tcmb=2.728,h=0.72,z=2.3)
    bias=2.38
    bias_err=.62
    kk=numpy.logspace(-4,2,51)
    civ_min_mass,civ_min_mass_upper,civ_min_mass_lower=get_minmass_for_bias(param,bias,bias_err,occupation=1)
    print civ_min_mass,civ_min_mass_upper,civ_min_mass_lower
    bias_kk=small_scale_bias(param,kk,civ_min_mass,occupation=1)

    numberdensity=NumberDensity.NumberDensity(param,massfraction='ShethTorman')
    delta_civ_kk=numberdensity.sigmaFunc.deltaSq(kk)*bias_kk

    rr=numpy.logspace(-2,3,31)
    KK,RR=numpy.meshgrid(kk,rr)
    KKRR=KK*RR
    sinc_KK=numpy.sin(KKRR)/KKRR
    corr=integrate.trapz(delta_civ_kk*sinc_KK,x=kk,axis=1)*numpy.log(10)

    from matplotlib import pyplot
    pyplot.loglog(rr,corr)
    pyplot.savefig('corr.png',fmt='png')

def civ_quasar_cross_correlation():
    param=Cosmology.CosmologyParam(omg_b=0.04,omg_c=0.22,omg_l=0.7,omg_k=0.0,sigma8=0.8,Tcmb=2.728,h=0.72,z=2.4)
    civ_bias,quasar_bias=2.3769,3.7073

    civ_minmass,qso_minmass = 53342669458.1,2.32453368181e+12
    civ_minmass_upper,qso_minmass_upper = 344668013021.0, 2.54668076862e+12
    civ_minmass_lower,qso_minmass_lower = 981882900.677, 2.102386595e+12
    kk=numpy.logspace(-4,3,151)
    numberdensity=NumberDensity.NumberDensity(param,massfraction='ShethTorman')

    civ_bias_kk=small_scale_bias(param,kk,civ_minmass,occupation=1)
    qso_bias_kk=small_scale_bias(param,kk,qso_minmass)
    delta_civ_qso_kk=numberdensity.sigmaFunc.deltaSq(kk)*civ_bias_kk*qso_bias_kk

    civ_bias_kk_upper=small_scale_bias(param,kk,civ_minmass_upper,occupation=1)
    qso_bias_kk_upper=small_scale_bias(param,kk,qso_minmass_upper)
    delta_civ_qso_kk_upper=numberdensity.sigmaFunc.deltaSq(kk)*civ_bias_kk_upper*qso_bias_kk_upper

    civ_bias_kk_lower=small_scale_bias(param,kk,civ_minmass_lower,occupation=1)
    qso_bias_kk_lower=small_scale_bias(param,kk,qso_minmass_lower)
    delta_civ_qso_kk_lower=numberdensity.sigmaFunc.deltaSq(kk)*civ_bias_kk_lower*qso_bias_kk_lower

    #temp=numpy.load('data/power_spectrum_cross.npy')
    #delta_civ_qso_kk=temp[1]
    numpy.save('data/power_spectrum_cross',(kk,delta_civ_qso_kk))
    numpy.save('data/power_spectrum_cross_upper',(kk,delta_civ_qso_kk_upper))
    numpy.save('data/power_spectrum_cross_lower',(kk,delta_civ_qso_kk_lower))

    rr=numpy.logspace(-3,3,101)
    delta_kk=delta_civ_qso_kk
    def integrand(x,r):
        return numpy.interp(10**x,kk,delta_kk)*numpy.sin(10**x *r)/(10**x *r)
    corr=numpy.array([integrate.romberg(integrand,numpy.log10(kk[0]),numpy.log10(kk[-1]),args=(rr_i,),divmax=20) for rr_i in rr])*numpy.log(10)
    numpy.save('data/corr_cross',(rr,corr))
    delta_kk=delta_civ_qso_kk_upper
    corr_upper=numpy.array([integrate.romberg(integrand,numpy.log10(kk[0]),numpy.log10(kk[-1]),args=(rr_i,),divmax=20) for rr_i in rr])*numpy.log(10)
    numpy.save('data/corr_cross_upper',(rr,corr_upper))
    delta_kk=delta_civ_qso_kk_lower
    corr_lower=numpy.array([integrate.romberg(integrand,numpy.log10(kk[0]),numpy.log10(kk[-1]),args=(rr_i,),divmax=20) for rr_i in rr])*numpy.log(10)
    numpy.save('data/corr_cross_lower',(rr,corr_lower))

    
if __name__=='__main__':
    param=Cosmology.CosmologyParam(omg_b=0.04,omg_c=0.22,omg_l=0.7,omg_k=0.0,sigma8=0.8,Tcmb=2.728,h=0.72,z=2.3)
    #print get_minmass_for_bias(param,4.0,occupation=1)
    #intergraltest(param,massfraction='PressSchecter')
    #civ_bias_corr()
    #qso_bias_corr()
    civ_quasar_cross_correlation()
