from numpy import sqrt,log
from numpy import power as pow
import Cosmology

def SQR(a):
    return a*a

class TransferFunction:
    def __init__(self,param):
        self.param=param
        omega_matter=param.omg_b+param.omg_c
        omega_baryon=param.omg_b
        omega_lambda=param.omg_l
        hubble=param.h
        redshift=param.z
        omega_hdm=0.0

        self.theta_cmb = 2.728/2.7   
        self.num_degen_hdm = 1  
        if omega_baryon<=0: omega_baryon=1e-5
        if omega_hdm<=0: omega_hdm=1e-5

        omega_curv = 1.0-omega_matter-omega_lambda
        self.omhh = omega_matter*SQR(hubble)
        obhh = omega_baryon*SQR(hubble)
        onhh = omega_hdm*SQR(hubble)
        f_baryon = omega_baryon/omega_matter
        self.f_hdm = omega_hdm/omega_matter
        f_cdm = 1.0-f_baryon-self.f_hdm
        self.f_cb = f_cdm+f_baryon
        f_bnu = f_baryon+self.f_hdm

        # Compute the equality scale. 
        z_equality = 25000.0*self.omhh/SQR(SQR(self.theta_cmb))  
        k_equality = 0.0746*self.omhh/SQR(self.theta_cmb)

        # Compute the drag epoch and sound horizon
        z_drag_b1 = 0.313*pow(self.omhh,-0.419)*(1+0.607*pow(self.omhh,0.674))
        z_drag_b2 = 0.238*pow(self.omhh,0.223)
        z_drag = 1291*pow(self.omhh,0.251)/(1.0+0.659*pow(self.omhh,0.828))*(1.0+z_drag_b1*pow(obhh,z_drag_b2))
        y_drag = z_equality/(1.0+z_drag)

        self.sound_horizon_fit = 44.5*log(9.83/self.omhh)/sqrt(1.0+10.0*pow(obhh,0.75))

        # Set up for the free-streaming & infall growth function 
        p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm))
        self.p_cb = 0.25*(5.0-sqrt(1+24.0*self.f_cb))

        omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+omega_matter*(1.0+redshift))
        omega_lambda_z = omega_lambda/omega_denom
        omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom
        self.growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/(pow(omega_matter_z,4.0/7.0)\
                -omega_lambda_z+(1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0))
        growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)-omega_lambda + \
            (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0))
        self.growth_to_z0 = self.growth_k0/growth_to_z0  
    
        # Compute small-scale suppression 
        alpha_nu = f_cdm/self.f_cb*(5.0-2.*(p_c+self.p_cb))/(5.-4.*self.p_cb)*  pow(1+y_drag,self.p_cb-p_c)*    \
            (1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/(1-0.193*sqrt(self.f_hdm*self.num_degen_hdm)+\
            0.169*self.f_hdm*pow(self.num_degen_hdm,0.2))*(1+(p_c-self.p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*self.p_cb))/(1+y_drag))
        self.alpha_gamma = sqrt(alpha_nu)
        self.beta_c = 1/(1-0.949*f_bnu)
        #Done setting scalar variables 
        self.hhubble = hubble   # Need to pass Hubble constant to TFmdm_onek_hmpc() 

    def __call__(self,kk):
        kk=kk*self.hhubble
        qq = kk/self.omhh*SQR(self.theta_cmb)

        # Compute the scale-dependent growth functions 
        y_freestream = 17.2*self.f_hdm*(1+0.488*pow(self.f_hdm,-7.0/6.0))*SQR(self.num_degen_hdm*qq/self.f_hdm)
        temp1 = pow(self.growth_k0, 1.0-self.p_cb)
        temp2 = pow(self.growth_k0/(1+y_freestream),0.7)
        growth_cb = pow(1.0+temp2, self.p_cb/0.7)*temp1
        growth_cbnu = pow(pow(self.f_cb,0.7/self.p_cb)+temp2, self.p_cb/0.7)*temp1

        # Compute the master function
        gamma_eff =self.omhh*(self.alpha_gamma+(1-self.alpha_gamma)/(1+SQR(SQR(kk*self.sound_horizon_fit*0.43))))
        qq_eff = qq*self.omhh/gamma_eff

        tf_sup_L = log(2.71828+1.84*self.beta_c*self.alpha_gamma*qq_eff)
        tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11))
        tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff))

        qq_nu = 3.92*qq*sqrt(self.num_degen_hdm/self.f_hdm)
        max_fs_correction = 1+1.2*pow(self.f_hdm,0.64)*pow(self.num_degen_hdm,0.3+0.6*self.f_hdm)/(pow(qq_nu,-1.6)+pow(qq_nu,0.8))
        tf_master = tf_sup*max_fs_correction

        # Now compute the CDM+HDM+baryon transfer functions 
        tf_cb = tf_master*growth_cb/self.growth_k0
        tf_cbnu = tf_master*growth_cbnu/self.growth_k0
        #print 'growth_cb/self.growth_k0 at redshift:',growth_cb/self.growth_k0
        #print 'growth_cb/growth_to_z0:',(growth_cb/self.growth_k0)*self.growth_to_z0
        return tf_cb*self.growth_to_z0


if __name__ == '__main__':
    param = Cosmology.CosmologyParam(omg_b=0.04,omg_c=0.22,omg_l=0.74,omg_k=0.,sigma8=0.9,Tcmb=2.728,h=0.72,n=1.,z=2.)
    tf = TransferFunction(param)
    print 'k=',0.1,' => ',tf(0.1)
    print 'k=',1.,' => ',tf(1.)
    print 'k=',10.,' => ',tf(10.)



