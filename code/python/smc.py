import numpy, pylab, matplotlib
import random, os

class Variables(object):
    """variables specifying the particulars of the trace and the preferences of the user."""
    def __init__(self, F, dt, 
                 x=None,   #stimulus (zeros vector of length T that has 1s in the frames that have high spike likelihood)
                 name='oopy', #name for plots/figures
                 Nparticles=99, # number of particles
                 Nspikehist=0, # number of spike history terms
                 condsamp=True, #use conditional sampler?
                 true_n = None,   #true spikes, if available
                 smc_iter_max = 3, # max number of iterations
                 est_c = True, # do we estimate tau_c, A, C_0?
                 est_t = True, # do we estimate tau_c?
                 est_n = True, # b,k
                 est_h = True, #w
                 est_F = True, #alpha, beta
                 showGraphs = True #show graphical results as we go
                 ):
        """
        initializes the variables. requires the Fluorescence trace and the frame length (dt), everything else has defaults.
        @param F: the Fluorescence trace itself.
        @param dt: seconds between timesteps
        @param x=None:   stimulus (zeros vector of length T that has 1s in the frames that have high spike likelihood)
        @param name='oopy': name for plots/figures
        @param Nparticles=99: number of particles
        @param Nspikehist=0: # number of spike history terms
        @param condsamp=True: #use conditional sampler?
        @param true_n = None:   #true spikes, if available
        @param smc_iter_max = 3: # max number of iterations
        @param est_c = True: # do we estimate tau_c, A, C_0?
        @param est_t = True: # do we estimate tau_c?
        @param est_n = True: # b,k
        @param est_h = True: #w
        @param est_F = True: #alpha, beta
        @param showGraphs = True: #show graphical results as we go
         
        """
       
        self.F = F
        self.T = len(F)
        self.dt = dt
        self.name=name
        self.Nparticles=Nparticles
        self.Nspikehist = Nspikehist
        self.condsamp=condsamp
        self.true_n = true_n
        self.smc_iter_max = smc_iter_max
        self.est_c = est_c
        self.est_t = est_t
        self.est_n = est_n
        self.est_h = est_h
        self.est_F = est_F
        self.showGraphs =showGraphs
        
        if(x==None): self.x = numpy.ones(len(T))
        else:self.x = x
        

class Parameters(object):
    """parameters of the model"""
    def __init__(self, 
                 V, 
                 tau_c = 1.0,
                 A=50.0,
                 C_0 = 0.,
                 C_init = 0.,
                 sigma_c = 0.1,
                 n=1.,
                 k_d = 200.,
                 k = None,
                 alpha=None,
                 beta=None, 
                 zeta=None,
                 gamma=None,
                 omega=-1.,
                 tau_h = 0.02,
                 sigma_h = 0.,
                 g=None,
                 sig2_h = None,
                 a=None,
                 sig2_c = None
                 ):
        """initializes model parameters.  needs a Variables object but everything else has default values.
            @param V: a Variables object initialized for this data set 
            @param tau_c = 1.: calcium decay time constant (seconds)
            @param A = 50.: change in [Ca++] after a spike (uM)
            @param C_0 = 0.: baseline [Ca++] (uM)
            @param C_init = 0.: initial [Ca++] (uM)
            @param sigma_c = 0.1: standard dev of noise (uM ?? )
            @param n = 1.: hill eq. exponent
            @param k_d = 200: hill coeff.
            @param k = 0.0001: linear filter related to spike rate
            @param alpha=None: scale for F (defaults to mean(F))
            @param beta = None: offset for F (defaults to min(F))
            @param zeta=None: constant variance (defaults to alpha/5.)
            @param gamma=None: scaled variance (defaults to zeta/5.)
            @param omega=-1,: weight?
            @param tau_h = 0.02: time constant?
            @param sigma_h = 0: standard dev of noise
            @param g=None: dt/tau_h, memoized for convenience
            @param sig2_h = None: square(sigma_h)*dt , memoized for convenience
            @param a=None: dt/tau_c, memoized for convenience
            @param sig2_c = None: square(sigma_c)*dt, memoized for convenience
        """
        
        self.lik = -numpy.inf #initial likelihood: uninitialized state is none too likely!
        
        self.V = V
        self.tau_c = tau_c
        self.A = A
        self.C_0 = C_0
        self.C_init = C_init
        self.sigma_c = sigma_c
        self.n=n
        self.k_d = k_d
        if(k==None): self.k = .001
        else: self.k = k
        self.kx = self.k*V.x; #we'll have to fix this to allow k to be non-scalar. later. 
        if(alpha==None): self.alpha = numpy.mean(V.F)
        else: self.alpha = alpha
        if(beta==None): self.beta = min(V.F)
        else: self.beta = beta
        if(zeta==None): self.zeta = self.alpha / 5.0
        else: self.zeta = zeta
        if(gamma==None): self.gamma = self.zeta / 5.0
        else: self.gamma = gamma
        self.omega = omega
        self.tau_h = tau_h
        self.sigma_h = sigma_h
        if(g==None): self.g = V.dt / self.tau_h
        else: self.g = g
        if(sig2_h==None): self.sig2_h = (self.sigma_h**2)*V.dt
        else:self.sig2_h = sig2_h
        if(a==None): self.a = V.dt / self.tau_c
        else: self.a = a
        if(sig2_c == None): self.sig2_c = (self.sigma_c ** 2)*V.dt
        else: self.sig2_c = sig2_c
                 
        
class Memoized(object):
    ''' creates and holds some convenience vectors and matrices that we'll reuse a lot.  
    '''
    def __init(self, vars):
        '''
        sets it up. 
        @param vars: a Variables object.  
        '''
        self.n_sampl = numpy.random.uniform(size = (vars.Nparticles, vars.T))
        self.C_sampl = numpy.random.uniform(vars.Nparticles, vars.T)
        self.oney = numpy.ones((vars.Nparticles, 1))
        self.zeroy = numpy.zeros((vars.Nparticles, 1))

class ObsLik(object):
    '''
    holds the observation likelihood parameters
    '''
    def __init__(self, vars, pars):
        '''
        initializes the likelihoods.
        @param vars: a Variables object 
        @param pars: a Parameters object 
        '''
        self.p_o = numpy.zeros(2,1)
        self.mu_o = numpy.zeros(2,1)
        self.sig2_o = numpy.zeros(1)
        self.p =  numpy.zeros(1)
        self.mu = numpy.zeros(1)
        self.sig2 = numpy.zeros(1)
        self.s = 1
        self.init_like(pars, vars.F(0))
        
    def init_lik(self, pars, initFluo):
        '''
        haven't gotten here yet.
        '''
        pass
        

class States(object):
    """ states of the model"""
    def __init__(self, vars, pars):
        '''
        set up the states. 
        @param vars: instance of a Variables object
        @param pars: instance of a Parameters object  
        '''
        self.p = numpy.zeros((vars.Nparticles, vars.T)) #rate
        n = numpy.zeros((V.Nparticles, V.T))
        self.n = n.astype('bool')     #spike counts
        self.C = P.C_init * numpy.ones((V.Nparticles, V.T)) #calcium -- probably to be rao-blackwellized away tho!
        self.w_f = (1.0 / V.Nparticles) * numpy.ones((V.Nparticles, V.T))) #forward particle weights
        self.w_b = (1.0 / V.Nparticles) * numpy.ones((V.Nparticles, V.T))) #backward particle weights
        
        #note: i think we shouldn't need this, or it shouldn't be a function of T_o, which we've gotten rid of. 
        # but i want it here commented out so that when i try to use S.Neff i remember why it doesn't exist.
        #self.Neff = (1.0 / V.Nparticles) * numpy.ones((1,V.T_o))  

def forward(vars, pars):
    '''
    the model is F_t = f(C) = alpha C^n/(C^n + k_d) + beta + e_t,
    where e_t ~ N[0, gamma*f(C)+zeta]
    
    @param vars: instance of a variables object 
    @param pars: instance of a parameters object 
    
    @return: instance of a States object  (simulation states)
    '''
    A = Memoized(vars)
    

def backward():
    pass