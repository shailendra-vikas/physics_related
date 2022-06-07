#!/usr/bin/env python
import Cosmology
import NumberDensity
import NFW
import numpy
from scipy import integrate


class SingleHalo:
    """                        dn(m)  (  m  )2           2
        P1(k)= Integrate( -------*(-----)  * |u(k|m)|,dm)
                           dm     ( rho )
    """
    def __init__(self,param):
        self.nd=NumberDensity.NumberDensity(param)
        self.sigma=NumberDensity.Sigma(param)
    
    def __call__(self,k):
        output=integrate.quad(lambda log_m_over_m_star:self.integrand(10**log_m_over_m_star*self.sigma.getM_star(),k)[0],-10,5)
        error=integrate.quad(lambda log_m_over_m_star:self.integrand(10**log_m_over_m_star*self.sigma.getM_star(),k)[1],-10,5)
        return (numpy.log(10)*output[0],numpy.log(10)*error[0])
    
    def integrand(self,m,k):
        nfw= NFW.NFW(m,self.sigma.getM_star())
        factor=(m/Cosmology.CRITICAL_DENSITY_Msun_Mpc)**2*nfw.u(k)**2
        return (factor*self.nd(m)[0],factor*self.nd(m)[1])
        
if __name__ == '__main__':
    param = Cosmology.CosmologyParam(omg_b=0.045,omg_c=0.255,omg_l=0.7,omg_k=0.,sigma8=0.93,Tcmb=2.728,h=0.7,n=1.)
    sh=SingleHalo(param)
    print 'k=',0.1,' => ',sh(0.1)
    print 'k=',1.,' => ',sh(1.)
    print 'k=',10.,' => ',sh(10.)
        
    
    
