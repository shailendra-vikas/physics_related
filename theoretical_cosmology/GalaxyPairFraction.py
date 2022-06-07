#!/usr/bin/env python
import scipy.integrate as integrator
import math
import matplotlib as mpl
mpl.use('ps')
import matplotlib.pyplot as mplpyp
import numpy
mpl.rc('text', usetex=True)


class halo_desity_profile:
	""" The density profile of the a Halo. 
		Assumed to be sperical symmetric 
		\rho(x)=\rho_avg * del *g(c)/3x(1+cx)^2
		where x=r/R_vir and c is concentration parameter.
		c=(9/1+z)*(m/m_*(z))^-0.13
		Note: 
			Profile is not normalized
			Mass is provided in unit of M*
		"""
	mass_over_mass_star =0.	
	def __init__(self,mass_ratio):
		self.mass_over_mass_star=mass_ratio
		
	def __call__(self,r):	#r is in unit of R_vir
		c=9.*self.mass_over_mass_star**-0.13
		if r > 1.:
			return 0.
		else:
			r= r > 1E-7 and r or 1E-7	
			return 1./(r*(1.+c*r)**2)

class baryon_profile_in_halo:
	""" \rho_b(r)= \rho_m(r) * r**del_gamma 
	    Mass is provided in the unit of M*
	"""
	del_gamma=0.
	halo_profile=None
	def __init__(self,mass_ratio,del_gamma):
		self.del_gamma = del_gamma
		self.halo_profile = halo_desity_profile(mass_ratio)
	
	def __call__(self,r):
		r= r > 1E-7 and r or 1E-7	
		return self.halo_profile(r)*r**self.del_gamma
		
class galaxy_cord_to_halo_cord_transformation:
	""" Take the cordinate of galaxy and convert it to halo cordinate """		
	def __init__(self):
		pass
	
	def __call__(self,gal_r,gal_theta,gal_halo_dist):
		return gal_halo_dist+gal_r*math.cos(gal_theta)
		
class fraction_galaxy_pair:
	""" The fraction of galaxy pair with distance greater than r """
	mass_over_mass_star =0.
	def __init__(self,mass_ratio):
		self.mass_over_mass_star=mass_ratio
		
	def __call__(self,gal_gal_dist,density_profile):
		dist_from_halo = galaxy_cord_to_halo_cord_transformation()
		integral = integrator.dblquad(\
			lambda gal_theta,gal_halo_dist: density_profile(dist_from_halo(gal_gal_dist,gal_theta,gal_halo_dist)),\
			0.,math.pi,lambda gal_halo_dist: 0.,lambda gal_halo_dist: 1.)
		factor = 2.*math.pi*gal_gal_dist**2
		return (factor*integral[0],factor*integral[1])
	
def plot_frac_cum_fraction(steps,density_profile,plots,ls,clr,lbl):
	frac = fraction_galaxy_pair(1.)
	gal_gal_dists = numpy.arange(0.,2.,2./float(steps))
	values = numpy.zeros(steps)
	errors = numpy.zeros(steps)
	cum_values= numpy.zeros(steps)
	cum_val=0
	for i in range(steps):
		values[i],errors[i] = frac(gal_gal_dists[i],density_profile)
		cum_val += values[i]
		cum_values[i]=cum_val
	maxval=numpy.amax(values)	
	plot1 = plots['Fraction'].add_subplot(111)
	plot1.errorbar(gal_gal_dists[:],values[:]/maxval,yerr=errors[:]/maxval,linestyle=ls,color=clr,label=lbl)
	plot1.set_ylabel(r'$f(r/R_{vir})$',size='large',style='normal',family='monospace')
	plot1.set_xlabel(r'$r/R_{vir}$',size='large',style='normal',family='monospace')
	plot1.grid(color='r',linestyle=':',linewidth=1)
	plot1.legend()
	plot2 = plots['CumFraction'].add_subplot(111)
	plot2.plot(gal_gal_dists[:],cum_values[:]/cum_values[steps-1],linestyle=ls,color=clr,label=lbl)
	plot2.set_ylabel(r'$F(r/R_{vir})$',size='large',style='normal',family='monospace')
	plot2.set_xlabel(r'$r/R_{vir}$',size='large',style='normal',family='monospace')
	plot2.grid(color='r',linestyle=':',linewidth=1)
	plot2.legend()
	return plots
	
def save_plots(plots):
	for key,value in plots.iteritems():
		value.savefig(key+'.ps')
		
	
if __name__=="__main__":
	plots = {'Fraction':mplpyp.figure(1),'CumFraction':mplpyp.figure(2)}
	plot_frac_cum_fraction(40,baryon_profile_in_halo(1.0,-1.),plots,'-','b',r'$\Delta \gamma = -1$')
	plot_frac_cum_fraction(40,baryon_profile_in_halo(1.0,0.),plots,'-','r',r'$\Delta \gamma = 0$')
	plot_frac_cum_fraction(40,baryon_profile_in_halo(1.0,1.),plots,'-','g',r'$\Delta \gamma = 1$')
	save_plots(plots)
	
