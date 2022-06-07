import NumberDensity
import Cosmology
import numpy
plot_it=False
def bias(nu,delta_c):
    overdensity=200
    y=numpy.log10(overdensity)
    A=1.0+0.24*y*numpy.exp(-(4./y)**4)
    a=0.44*y-0.88
    B=0.183
    b=1.5
    C=0.019+0.107*y+0.19*numpy.exp(-(4./y)**4)
    c=2.4
    return 1-A*nu**a/(nu**a+delta_c**a)+B*nu**b+C*nu**c

def mass_for_bias(target_bias):
    param=Cosmology.CosmologyParam(omg_b=0.045,omg_c=0.215,omg_l=0.74,omg_k=0.,sigma8=0.80,Tcmb=2.728,h=0.72,n=1.,z=2.3)
    numberDensity=NumberDensity.NumberDensity(param,massfraction='ShethTorman')
    R_min,R_max=.005,4
    current_bias=0.0
    while numpy.abs(current_bias/target_bias -1) > .00001:
        current_bias=numberDensity.halobias(numberDensity.del_sc,numberDensity.sigmaFunc((R_min+R_max)/2.))
        #current_bias=bias(numberDensity.del_sc/numberDensity.sigmaFunc((R_min+R_max)/2.),numberDensity.del_sc)
        if current_bias <target_bias:
            R_min=(R_min+R_max)/2.
        else:
            R_max=(R_min+R_max)/2.
        #print R_min,R_max,current_bias

    mass_range=numpy.logspace(11,13,101)
    bias_values=numberDensity.halo_bias(mass_range)
    if plot_it:
        from matplotlib import pyplot
        pyplot.plot(8.3E11,2.38,'sb',label=r'$M_\mathrm{h,CIV} \simeq 8\times 10^{11} h^{-1} M_\odot$')
        pyplot.plot(4.74E12,3.71,'ok',label=r'$M_\mathrm{h,QSO}=4.7\times 10^{12} h^{-1} M_\odot$')
        pyplot.axhspan(2.38-.62,2.38+.62,xmin=(numpy.log10(1.841250e+11)-11)/2,xmax=(numpy.log10(2.175568e+12)-11)/2,color='#aaaaff')
        pyplot.axhspan(3.71-.09,3.71+.09,xmin=(numpy.log10(4.353676e+12)-11)/2,xmax=(numpy.log10(5.157296e+12)-11)/2,color='#aaaaaa')
        pyplot.legend()
        pyplot.semilogx(mass_range,bias_values,color='k')
        #pyplot.text(1e9,5,r'$M_\mathrm{halo}=%1.2E h^{-1}$ M$_\odot$'%(Cosmology.RadiusToMass(param,(R_min+R_max)/2.)),fontsize=20)
        #pyplot.text(2e8,2,r'$M_\mathrm{h,CIV}=1.14\times 10^{11} h^{-1} M_\odot$',fontsize=20,color='r')
        #pyplot.axhspan(2.38-.62,2.38+.62,facecolor=(.7,1.,.7))
        #pyplot.axhline(2.38,xmin=.4,xmax=.7,color='red',ls='--',lw=2)
        #pyplot.axvline(1.14e11,ymin=0.05,ymax=.3,color='red',ls='--',lw=2)
        #pyplot.text(5e9,5,r'$M_\mathrm{h,QSO}=1.88\times 10^{12} h^{-1} M_\odot$',fontsize=20,color='b')
        #pyplot.axhline(3.71,xmin=.5,xmax=.8,color='blue',ls='--',lw=2)
        #pyplot.axvline(1.884314e+12,ymin=.1,ymax=.4,color='blue',ls='--',lw=2)
        pyplot.xlabel('Mass $h^{-1} M_\odot$')
        pyplot.ylabel('Bias')
        pyplot.savefig('bias_mass.eps',fmt='eps')


    return Cosmology.RadiusToMass(param,(R_min+R_max)/2.) 

if __name__=='__main__':
    print '%e'%mass_for_bias(2.38) 
    print '%e'%mass_for_bias(2.38-.62) 
    print '%e'%mass_for_bias(2.38+.62) 
    #print '%e'%mass_for_bias(3.71+.09) 
    print '%e'%mass_for_bias(3.71) 
    #print mass_for_bias(2.37691019+.81692762) 
    #print mass_for_bias(2.37691019-.61692762) 
    #mass_for_bias(3.70731768) 


