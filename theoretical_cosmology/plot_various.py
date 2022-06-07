import NFW
import numpy
import Cosmology
import SingleHalo
import NumberDensity
from matplotlib import pyplot

def plot_triangle(params):
    k=numpy.logspace(-3,1,100)

    pyplot.figure()
    for i in range(len(params)):
        bigT=NumberDensity.BigTriangleSq(params[i])
        pyplot.loglog(k,numpy.sqrt(bigT(k)),'-')
    pyplot.xlabel('k')
    pyplot.ylabel(r'$\Delta$')
    pyplot.savefig('triag.png',fmt='png')

def plot_sigma(params):
    R=numpy.logspace(-3,1.2,50)
    pyplot.figure()
    for i in range(len(params)):
        sigma=NumberDensity.Sigma(params[i])
        pyplot.loglog(R,sigma(R),'-')
    pyplot.xlabel('R')
    pyplot.ylabel(r'$\sigma$')
    pyplot.savefig('sigma.png',fmt='png')

    pyplot.figure()
    for i in range(len(params)):
        sigma=NumberDensity.Sigma(params[i])
        M=Cosmology.RadiusToMass(params[i],R)
        pyplot.loglog(M,sigma(R),'-')
    pyplot.xlabel('M')
    pyplot.ylabel(r'$\sigma$')
    pyplot.savefig('sigma2.png',fmt='png')

def plot_bias_nu(params):
    M=numpy.logspace(5,18,25)
    pyplot.figure()
    for i in range(len(params)):
        #numberdensity=NumberDensity.NumberDensity(params[i],massfraction='ShethTorman')
        numberdensity=NumberDensity.NumberDensity(params[i])
        R=Cosmology.MassToRadius(params[i],M)
        nu=(numberdensity.del_sc/numberdensity.sigmaFunc(R))
        pyplot.loglog(nu,numberdensity.halo_bias(M))
    pyplot.xlabel(r'$\nu$')
    pyplot.ylabel(r'b')
    #pyplot.xlim((.5,5))
    #pyplot.ylim((.5,10))
    pyplot.savefig('bias.png',fmt='png')

def plot_bias_mass(params):
    M=numpy.logspace(5,18,25)
    pyplot.figure()
    for i in range(len(params)):
        #numberdensity=NumberDensity.NumberDensity(params[i],massfraction='ShethTorman')
        numberdensity=NumberDensity.NumberDensity(params[i])
        pyplot.loglog(M,numberdensity.halo_bias(M))
    pyplot.xlabel('Mass')
    pyplot.ylabel(r'b')
    #pyplot.xlim((.5,5))
    #pyplot.ylim((.5,10))
    pyplot.savefig('bias_mass.png',fmt='png')
    


def plot_number_density(param):
    nd=NumberDensity.NumberDensity(param)
    mass=numpy.logspace(12,18,15)
    numbers=numpy.array([nd(mass_i)[0] for mass_i in mass])
    print mass
    print numbers
    pyplot.figure()
    pyplot.loglog(mass,mass*numbers,'-')
    pyplot.savefig('number_density.png',fmt='png')


def plot_nfw(param):
    m=1E15
    sigma=NumberDensity.Sigma(param)
    nfw= NFW.NFW(m,sigma.getM_star())
    k=numpy.logspace(-2,2,100)
    values= numpy.array([nfw.u(k_i) for k_i in k])
    pyplot.figure()
    pyplot.loglog(k,values,'-')
    pyplot.savefig('nfw.png',fmt='png')

    
def plot_singlehalo_power(param):
    sh=SingleHalo.SingleHalo(param)
    k=numpy.logspace(0,2,3)
    print sh(1.0)       
    #pyplot.figure()
    #pyplot.loglog(k,numpy.array([k_i**3*sh(k_i)[0]/(2*numpy.pi**2) for k_i in k]),'-')
    #pyplot.savefig('single.png',fmt='png')


if __name__=='__main__':
    param=Cosmology.CosmologyParam(omg_b=0.04,omg_c=0.22,omg_l=0.74,omg_k=0.0,sigma8=0.9,Tcmb=2.7,h=0.72,z=0.0)
    param2=Cosmology.CosmologyParam(omg_b=0.04,omg_c=0.22,omg_l=0.74,omg_k=0.0,sigma8=0.9,Tcmb=2.7,h=0.72,z=2.0)
    params=[param,param2]
    plot_triangle(params)
    plot_sigma(params)
    plot_bias_nu(params)
    plot_bias_mass(params)
    #plot_number_density(param)
    #plot_nfw(param)
    #plot_singlehalo_power(param)
