#!/usr/bin/env python
import numpy
import Cosmology
from scipy.special import sici

class  NFWProfile:
    """
        This class store all the methods for NFW profile
    """
    def __init__(self,param,mass,m_star):
        self.mass=mass
        self.c_bar=(9./(1+param.z))*(mass/m_star)**-0.13
        self.r_vir=Cosmology.MassToRadius(param,mass)
        self.r_s=self.r_vir/self.c_bar
        self.factor= numpy.log(1.+self.c_bar)-self.c_bar/(1.+self.c_bar)
        
    def __call__(self,r):
        """
            Input:r/Rvir
            Output: Unit in avg density
        """
        self.rhos= self.mass/(4*numpy.pi*self.factor*self.r_s**3)
        r_over_rs=r/self.r_s
        return self.rhos/(r_over_rs *(1+r_over_rs)**2)
    
    def u(self,k):
        """
            Input is k*Rvir
            u(k|m)= 4pi*rho_s*r_s^3{                                    sin(ckr_s)
                --------------{sin(kr_s)[Si((1+c)kr_s)-Si(kr_s)]-  ----------- + cos(kr_s)[Ci((1+c)kr_s)-Ci(kr_s)]
                      m       {                                     (1+c)kr_s
        """
        krs=k*self.r_s
        (si_k,ci_k)=sici(krs)
        (si_kc,ci_kc)=sici((1.+self.c_bar)*krs)
        return (numpy.sin(krs)*(si_kc-si_k)-numpy.sin(self.c_bar*krs)/((1.+self.c_bar)*krs)+numpy.cos(krs)*(ci_kc-ci_k))/self.factor
    
if __name__ == '__main__':
    param=Cosmology.CosmologyState(omg_b=0.04,omg_c=0.22,omg_l=0.74,omg_k=0.0,sigma8=0.9,Tcmb=2.7,h=0.72,z=2.3)
    r= 10**numpy.arange(-2.,0.,.1)
    nfw= NFWProfile(param,1E16,1E16)
    k=numpy.logspace(-2,2,11)
    print nfw(r)
    print nfw.u(k)
    
