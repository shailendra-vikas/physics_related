import scipy.integrate as integrate
import numpy 
"""
    The module for calculating different cosmological distance. Based on the paper Hogg(2000) http://arxiv.org/abs/astroph/9905116

    The basic input parameter is Omega_lambda(Ol), Omega_matter(Om), hubble parameter (h) and equation of state of dark energy(w)

"""
class Cosmology:
    """
        The class containing cosmological distance methods.
    """
    def __init__(self,Om=0.3,Ol=0.7,w=-1.0,h=None):
        self.Om,self.Ol,self.w,self.Ok = Om,Ol,w,1.0-Om-Ol
        self.h  = h is None and 1.0 or h
        self.Dh      = 2997.9      # h^-1 Mpc
        self.tH      = 9.778       #h^-1 Gyr
    
    def E_inv(self,z):  # E(z)=H(z)/H0
        return 1.0/numpy.sqrt(self.Om*(1.+z)**3 + self.Ok*(1.+z)**2 + self.Ol*(1.+z)**(3*(1.+self.w)))

    def T(self,z): 
        return self.E_inv(z)/(1.+z)

    """ Lookback time  -> Difference between the age of the Universe now and the age at z"""
    def Tl(self, z):
        z=numpy.asarray(z)
        integral=lambda z:integrate.romberg(self.T, 0, z)
        vec_integral=numpy.vectorize(integral)
        return self.tH*vec_integral(z)/self.h #Gyr

    """ Line of sight comoving distance -> Remains constant with epoch if objects are in the Hubble flow"""
    def Dc(self, z):
        z=numpy.asarray(z)
        integral=lambda z:integrate.romberg(self.E_inv, 0, z)
        vec_integral=numpy.vectorize(integral) # To enable using numpy broadcast
        return self.Dh*vec_integral(z)/self.h #Mpc

    """ Transverse comoving distance -> At same redshift but separated by angle dtheta; Dm * dtheta is transverse comoving distance"""
    def Dm(self, z):
        dc=self.Dc(z)
        x=dc*numpy.sqrt(numpy.abs(self.Ok))*self.h/self.Dh
        #Used taylor expansion for small x for better accuracy
        if self.Ok >=0.0:
            return numpy.where(x<1E-5,
                dc*(1.0+numpy.sign(self.Ok)*x**2/6.+x**4/120.),
                dc*numpy.sinh(x)/x) #Mpc
                
        return numpy.where(x<1E-5,
                dc*(1.0+numpy.sign(self.Ok)*x**2/6.+x**4/120.),
                dc*numpy.sin(x)/x) #Mpc

    """ Angular diameter distance ->  Ratio of an objects physical transvserse size to its angular size in radians"""
    def Da(self, z):
        z=numpy.asarray(z)
        return self.Dm(z)/(1.+z) #Mpc

    """ Luminosity distance ->  Relationship between bolometric flux and bolometric luminosity"""
    def Dl(self, z):
        z=numpy.asarray(z)
        return (1.+z)*self.Dm(z) #Mpc

    """ Distance modulus ->  Recall that Dl is in mpc"""
    def DistMod(self, z):
        z=numpy.asarray(z)
        return 5.*numpy.log10(self.Dl(z))+25 #Mag

    """ Volume -> The volume of the sphere centered on us observed"""
    def volume(self,z):
        z=numpy.asarray(z)
        return self.E_inv(z)*self.Dh*self.Dm(z)**2 #Mpc^3

if __name__ == '__main__':
    c = Cosmology(Om=0.3,Ol=0.6,w=-1.0)
    z1=0.5
    print 'For z = %.2f:' % (z1)
    print 'Lookback time                  %.2f Gyr' %c.Tl(z1)
    print 'Comoving L.O.S. Distance (w)   %.2f Mpc' %c.Dc(z1)
    print 'Angular diameter distance      %.2f Mpc' %c.Da(z1)
    print 'Luminosity distance            %.2f Mpc' %c.Dl(z1)
    print 'Distance modulus               %.2f mag' %c.DistMod(z1)
    print 'Volume                         %.2f Mpc^3'%c.volume(z1)
