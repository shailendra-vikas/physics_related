import numpy
"""
This file is interface to read correlation data provided in the data file 
and interpolate and provide correlation value for any given distance of 
seperation and the redshift.
"""



r_value,corr_value=numpy.load('/n/home10/vikas/myproduct/theorycorr/data/corr_cross.npy')
r_value_upper,corr_value_upper=numpy.load('/n/home10/vikas/myproduct/theorycorr/data/corr_cross_upper.npy')
r_value_lower,corr_value_lower=numpy.load('/n/home10/vikas/myproduct/theorycorr/data/corr_cross_lower.npy')
omg_b,omg_c,omg_l,h=0.04,0.22,0.7,0.72

def SQR(a):
    return a*a

def growth(redshift):
    omega_matter=omg_b+omg_c
    hubble=h
    theta_cmb = 2.728/2.7
    omega_lambda=omg_l
    omhh = omega_matter*SQR(hubble)

    z_equality = 25000.0*omhh/SQR(SQR(theta_cmb))
    omega_curv = 1.0-omega_matter-omega_lambda
    omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+omega_matter*(1.0+redshift))
    omega_lambda_z = omega_lambda/omega_denom
    omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom
    growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/(pow(omega_matter_z,4.0/7.0)\
                    -omega_lambda_z+(1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0))
    growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)-omega_lambda + \
                (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0))

    growth_to_z0 = growth_k0/growth_to_z0
    return growth_to_z0

def get_corr(r,zqso):
    # Unit of r is Mpc but unit of r_value is h^-1 Mpc 
    if numpy.isscalar(r) or len(r.shape)==1:
        ret_corr=numpy.interp(r,r_value/h,corr_value,right=0.0)
        err_upper=numpy.interp(r,r_value_upper/h,corr_value_upper,right=0.0)-ret_corr
        err_lower=ret_corr-numpy.interp(r,r_value_lower/h,corr_value_lower,right=0.0)
        growth_factor=(growth(zqso)/growth(2.4))**2
        if numpy.isscalar(growth_factor):
    	    return growth_factor*ret_corr,err_upper,err_lower
        GROWTH_FACTOR,RET_CORR=numpy.meshgrid(growth_factor,ret_corr)
        GROWTH_FACTOR,ERR_UPPER=numpy.meshgrid(growth_factor,err_upper)
        GROWTH_FACTOR,ERR_LOWER=numpy.meshgrid(growth_factor,err_lower)
        return GROWTH_FACTOR*RET_CORR,ERR_UPPER,ERR_LOWER
    elif len(r.shape)==2 and len(zqso.shape)==2:
        RET_CORR = numpy.interp(r.flatten(),r_value/h,corr_value,right=0.0).reshape(r.shape)
        ERR_UPPER = numpy.interp(r.flatten(),r_value_upper/h,corr_value_upper,right=0.0).reshape(r.shape) - RET_CORR
        ERR_LOWER = RET_CORR - numpy.interp(r.flatten(),r_value_lower/h,corr_value_lower,right=0.0).reshape(r.shape)
        growth_value = growth(zqso[0,:])
        return growth_value*RET_CORR,ERR_UPPER,ERR_LOWER

