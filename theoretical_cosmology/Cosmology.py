#!/usr/bin/env python
import math

VOLUMN_FACTOR=4.*math.pi/3.   #4*pi/3
SOLAR_MASS=1.989E33 #Unit grams
CRITICAL_DENSITY=1.8791E-29 #Unit h^2 g/cm^3
MPC_CM=3.0856E24
CRITICAL_DENSITY_Msun_Mpc=CRITICAL_DENSITY*MPC_CM**3/SOLAR_MASS #unit h^2 M_sun/Mpc^3
OVERDENSITY=200.


class CosmologyParam:
    """This class contains the parameters for cosmology.
    Should do the sanity check for parameter"""
    
    def __init__(self,omg_b=0.05,omg_c=0.25,omg_l=0.7,omg_k=0.0,sigma8=0.8,h=0.6,Tcmb=2.7,n=1.,z=0.0):
        self.omg_b = omg_b
        self.omg_c = omg_c
        self.omg_l = omg_l
        self.omg_k = omg_k
        self.sigma8 = sigma8
        self.Tcmb =Tcmb
        self.h = h
        self.n = n  # P(k) propotional k^n
        self.overdensity=200.
        self.z=z
    
    def write(self):
        print 'omg_b','=>',self.omg_b
        print 'omg_c','=>',self.omg_c
        print 'omg_l','=>',self.omg_l
        print 'omg_k','=>',self.omg_k
        print 'sigma8','=>',self.sigma8
        print 'Tcmb','=>',self.Tcmb
        print 'h','=>',self.h
        print 'n','=>',self.n
        print 'Overdensity','=>',self.overdensity
        print 'z','=>',self.z

def E_sq(param):
    one_plus_z=(1+param.z)
    return one_plus_z**3*(param.omg_b+param.omg_c)+one_plus_z**2*param.omg_k+param.omg_l

def MassToRadius(param,mass):
    """
        Converts Mass to Radius. 
        Input unit h^{-1} M_sun
        Output unit R (h^{-1} Mpc)
    """
    return (mass/(VOLUMN_FACTOR*CRITICAL_DENSITY_Msun_Mpc*(param.omg_c+param.omg_b)))**(1./3.)

def RadiusToMass(param,R):
    """
        Converts Radius to Mass. 
        Input unit R (h^{-1} Mpc)
        Output unit h^{-1} M_sun
    """
    return VOLUMN_FACTOR*CRITICAL_DENSITY_Msun_Mpc*(param.omg_c+param.omg_b)*R**3

if __name__ == "__main__":
    testparam1 =  CosmologyParam(omg_b=0.02,omg_c=0.2,omg_l=0.7,omg_k=0.1,sigma8=0.8,Tcmb=2.7,h=0.6)
    testparam1.write()
    teststate1 = CosmologyState(testparam1,z=10.0)
    print hasattr(teststate1,'z')
    teststate1.write()
    teststate2 = CosmologyState(omg_b=0.01,omg_c=0.125,omg_l=0.17,omg_k=1.0,sigma8=1.9,Tcmb=2.7,h=0.6,z=5.0)
    teststate2.write()


