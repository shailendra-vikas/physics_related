#!/usr/bin/env python

import numpy
from scipy import integrate
from scipy import special
import Cosmology 
import TransferFunc

class WindowFunction:
    """ Window function 
        Unit of R is h^-1 Mpc
    """ 
    class realTopHat:
        def __init__(self):
            pass

        def __call__(self,x):
            return numpy.where(x < 0.001,1.-x**2/10.,3*(numpy.sin(x)-x*numpy.cos(x))/x**3)
            
        def derivative(self,x):
            return numpy.where(x <.01,3.*(x/15.-(39./(42.*120.))*x**3.),3*(x**2*numpy.sin(x)+3*x*numpy.cos(x)-3*numpy.sin(x))/x**4)

    def __init__(self,windowtype=None):
        if windowtype is None or windowtype == 'realTopHat':
            self.window = self.realTopHat()

    def __call__(self,x):
        return self.window(x)
        
    def derivative(self,x): # Derivative w.r.t R
        return self.window.derivative(x)

class Delta:
    """ "Big trianglesq" function = k^3+n*|T(k)|^2 
        Assume the input k is in unit of h Mpc^-1
        Normalized to have value 1 at k=1
    """
    def __init__(self,cosmoparam):
        self.cosmoparam=cosmoparam
        self.tf = TransferFunc.TransferFunction(cosmoparam)
        n_tilda=self.cosmoparam.n-1
        omega_o=self.cosmoparam.omg_c+self.cosmoparam.omg_b
        delta_h=1.94E-5 * omega_o**(-.785-0.05*numpy.log(omega_o)) * numpy.exp( -.95*n_tilda -0.169* n_tilda**2)
        self.norm = delta_h**2 * 3000.**(3.+self.cosmoparam.n)
        

    def __call__(self,k_hMpc):
        return self.norm*k_hMpc**(3.+self.cosmoparam.n)*self.tf(k_hMpc)**2 

class BigTriangleSq(Delta):
    def __init__(self,cosmoparam):
        Delta.__init__(self,cosmoparam)
        zero_cosmoparam=Cosmology.CosmologyParam(omg_b=cosmoparam.omg_b,omg_c=cosmoparam.omg_c,omg_l=cosmoparam.omg_l,omg_k=cosmoparam.omg_k,sigma8=cosmoparam.sigma8,Tcmb=cosmoparam.Tcmb,h=cosmoparam.h,n=cosmoparam.n,z=0.0)
        temp_delta=Delta(zero_cosmoparam)
        self.wf=WindowFunction()
        integral=integrate.quad(lambda logx: temp_delta((10**logx)/8.)*self.wf(10**logx)**2,-6.,7.)
        self.newnorm=cosmoparam.sigma8**2/(numpy.log(10)*integral[0])

    def __call__(self,k_hMpc):
        return self.newnorm*Delta.__call__(self,k_hMpc)


class Sigma:
    """
     Calculate sigma for given cosmological parameter ans mass.
        R is in the unit of h^{-1}Mpc
        sigma(m) = Integral( d(lnk)*big_trangle^2 * Window(kr)^2)
    """ 
    def __init__(self,cosmoparam):
        self.cosmoparam=cosmoparam
        self.deltaSq=BigTriangleSq(cosmoparam)
        self.M_star=None
        self.wf=WindowFunction()

    def __call__(self,R):
        if numpy.isscalar(R):
            output=integrate.quad(lambda logx:self.deltaSq(10**logx/R)*self.wf(10**logx)**2,-6.,5.)
            return numpy.sqrt(numpy.log(10)*output[0])
        output=numpy.array([integrate.quad(lambda logx:self.deltaSq(10**logx/R_i)*self.wf(10**logx)**2,-6.,5.) for R_i in R ])
        return numpy.sqrt(numpy.log(10)*output[:,0])
        
    def sigma_into_d_sigma_by_d_R(self,R):
        """
            calulcate sigma*d(sigma)/dR =Integrate(2*deltaSq(k)*WindowFunction(k,R)*d(WindowFunction)/dR,log_k)
        """
        if numpy.isscalar(R):
            output=integrate.quad(lambda x:self.deltaSq(x/R)*self.wf(x)*self.wf.derivative(x),0.,30.)
            return output[0]/R
        output=numpy.array([integrate.quad(lambda x:self.deltaSq(x/R_i)*self.wf(x)*self.wf.derivative(x),0.,30.) for R_i in R])
        return output[:,0]/R
        
    def getM_star(self):
        """
            Try to binary search
        """
        if self.M_star is not None:
            return self.M_star
            
        log_Rmin=-4
        log_Rmax=2
        del_sc=1.686

        while((log_Rmax-log_Rmin) > 0.0001):
            if self(10**((log_Rmin+log_Rmax)/2.)) > del_sc:
                log_Rmin=(log_Rmin+log_Rmax)/2.
            else:
                log_Rmax=(log_Rmin+log_Rmax)/2.
            
        self.M_star= Cosmology.RadiusToMass(self.cosmoparam,10**((log_Rmin+log_Rmax)/2.))
        return self.M_star

class NumberDensity:
    """
        n(m)= (\rho^bar/m^2)*d ln(sigma)/d ln(m) nu*f(nu)
    """
    class NumberDensityFunction:
        """
            nu*f(nu) = A(p)*(1+(q*nu)^-p)(q*nu/2pi)^.5*exp(-q*nu/2)
            for pressschecter p=0 q=1
        """
        def __init__(self,p=0.3,q=.75):
            self.A_p=1./(1+special.gamma(.5-p)/(numpy.sqrt(numpy.pi)*2**p))
            self.p,self.q=p,q

        def __call__(self,nu):
            qnu=self.q*nu
            return self.A_p*(1.+1./qnu**self.p)*numpy.sqrt(qnu/(2.*numpy.pi))*numpy.exp(-qnu/2.)

    class HaloBias:
        """
            b(p,q)= 1+ (qnu-1)/del_sc + (2p/del_sc)/(1+ qnu**p)
        """
        def __init__(self,p=0.3,q=.75):
            self.p,self.q=p,q

        def __call__(self,del_sc,sigma_val):
            nu=(del_sc/sigma_val)**2
            qnu=self.q*nu
            return 1.+(qnu-1)/del_sc + 2*self.p/(del_sc*(1+qnu**self.p))
            
    def __init__(self,param,massfraction='PressSchecter'):
        self.param=param
        self.massfraction= self.NumberDensityFunction(p=0.3,q=0.75) if massfraction=='ShethTorman' else self.NumberDensityFunction(p=0,q=1.)
        self.halobias= self.HaloBias(p=0.3,q=0.75) if massfraction=='ShethTorman' else self.HaloBias(p=0,q=1.)
        self.sigmaFunc=Sigma(param)
        self.del_sc=1.686
        
    def __call__(self,mass):
        """
            Assumed unit is in h^-1 M_sun
        """
        R=Cosmology.MassToRadius(self.param,mass)
        sigma=self.sigmaFunc(R)
        nu=(self.del_sc/sigma)**2
        sigma_into_d_sigma_by_d_R=numpy.abs(self.sigmaFunc.sigma_into_d_sigma_by_d_R(R))
        number_density_factors=2/(mass*4*numpy.pi*R**2*sigma**2)
        return number_density_factors*sigma_into_d_sigma_by_d_R * self.massfraction(nu)

    def halo_bias(self,mass):
        R=Cosmology.MassToRadius(self.param,mass)
        sigma=self.sigmaFunc(R)
        return self.halobias(self.del_sc,sigma)

        
if __name__== '__main__':
    param=Cosmology.CosmologyState(omg_b=0.045,omg_c=0.255,omg_l=0.7,omg_k=0.,sigma8=0.93,Tcmb=2.728,h=0.7,n=1.,z=0.0)
    param1=Cosmology.CosmologyState(omg_b=0.045,omg_c=0.255,omg_l=0.7,omg_k=0.,sigma8=0.93,Tcmb=2.728,h=0.7,n=1.,z=2.0)
    R=numpy.logspace(-1,1,11)
    sigma=Sigma(param)
    sigma1=Sigma(param1)
    print sigma.getM_star()
    print sigma1.getM_star()
    #print sigma(R)
    #print sigma.sigma_into_d_sigma_by_d_R(R)
    mass=numpy.logspace(10,16,15)
    nd=NumberDensity(param)
    dn_dm=nd(mass)[0]
    print dn_dm.shape
    rho_matter=Cosmology.CRITICAL_DENSITY_Msun_Mpc*.3
    print mass.shape
    print sigma(R)

    from matplotlib import pyplot
    pyplot.figure()
    pyplot.loglog(mass,mass**2*dn_dm/rho_matter)
    pyplot.savefig('nd.png',fmt='png')

