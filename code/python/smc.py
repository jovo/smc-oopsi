import numpy, pylab, matplotlib
import random, os, bisect, re


class Variables(object):
    """variables specifying the particulars of the trace and the preferences of the user."""
    def __init__(self, F, dt, 
                 x=None,   #stimulus (zeros vector of length T that has 1s in the frames that have high spike likelihood)
                 name='oopy', #name for plots/figures
                 Nparticles=999, # number of particles
                 Nspikehist=0, # number of spike history terms
                 condsamp=False, #use conditional sampler?
                 true_n = None,   #true spikes, if available
                 smc_iter_max = 3, # max number of iterations
                 est_c = True, # do we estimate tau_c, A, C_0?
                 est_t = True,  # do we estimate tau_c?
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
        @param Nspikehist=0:  number of spike history terms
        @param condsamp=True: use conditional sampler?
        @param true_n = None:   true spikes, if available
        @param smc_iter_max = 3:  max number of iterations
        @param est_c = True: do we estimate tau_c, A, C_0?
        @param est_t = True: do we estimate tau_c?
        @param est_n = True: b,k
        @param est_h = True: w
        @param est_F = True: alpha, beta
        @param showGraphs = True: show graphical results as we go
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
        
        if(x==None): self.x = numpy.ones(self.T)
        else:self.x = x
        

class Parameters(object):
    """parameters of the model"""
    def __init__(self, 
                 V, 
                 tau_c = 1.0,
                 A=50.0,
                 C_0 = 5e-5,
                 C_init = 5e-5,
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
                 sigma_h = 5e-5,
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
            @param k = None: for now, ansatz for spikes per second. 
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
        self.k=k
        if(k==None): self.k = .0001
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
    def __init__(self, vars, pars):
        '''
        sets it up. 
        @param vars: a Variables object.  
        '''
        self.n_sampl = numpy.random.uniform(size = (vars.Nparticles, vars.T))
        self.C_sampl = numpy.random.uniform(size=(vars.Nparticles, vars.T))
        self.oney = numpy.ones((vars.Nparticles, 1))
        self.zeroy = numpy.zeros((vars.Nparticles, 1))
        self.U_sampl = numpy.random.uniform(size=(vars.Nparticles, vars.T))
        
        
        diffs = 1.0 / vars.Nparticles
        ints = numpy.arange(0,1,diffs)
        shapedInts = numpy.repeat(ints, vars.T).reshape((vars.Nparticles,vars.T)).T #the second .T is a transpose operation.
        self.U_resamp = shapedInts  + diffs*numpy.random.uniform(size = (vars.T,vars.Nparticles))
        
        
#here's the original code,  and ints was to Nparticles+1. 
#so instead, take ints only to Nparticles, 
#and since V.T_o is 1, we're tiling exactly once, so no more repmat
#A.U_resamp  = repmat(ints(1:end-1),V.T_o,1)+diffs*rand(V.T_o,V.Nparticles); % resampling matrix

        self.epsilon_c = numpy.sqrt(pars.sig2_c) * numpy.random.normal(size=(vars.Nparticles, vars.T))

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
        self.p_o = numpy.zeros((2,1))
        self.mu_o = numpy.zeros((2,1))
        self.sig2_o = numpy.zeros(1)
        self.p =  numpy.zeros(1)
        self.mu = numpy.zeros(1)
        self.sig2 = numpy.zeros(1)
        self.s = 1  # s = v.freq = "intermittent observation frequency". it's always going to be 1 for a while.
        self.V = vars
        self.P = pars
        
        self.init_lik()
        
    def init_lik(self):
        '''
        get the mean (mu1) and variance (sig1) for P[C_t | F_t]
        '''
        F = self.V.F #for brevity
        P = self.P   #for brevity
        
        finv = numpy.power( ((P.k_d * (F-P.beta)) / (P.alpha - F + P.beta)),
                            1/P.n)
        
        mu1 = finv #copying josh's variable names
        if( (mu1 > 0).all() and (numpy.imag(mu1)==0).all()):
            self.sig1 = -1 / (-(-P.alpha * mu1 ** P.n * P.n / mu1 / (mu1 ** P.n + P.k_d) + P.alpha * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1) ** 2 / (P.gamma * mu1 ** P.n / (mu1 ** P.n + P.k_d) + P.zeta) + 2 * (F - P.alpha * mu1 ** P.n / (mu1 ** P.n + P.k_d) - P.beta) / (P.gamma * mu1 ** P.n / (mu1 ** P.n + P.k_d) + P.zeta) ** 2 * (-P.alpha * mu1 ** P.n * P.n / mu1 / (mu1 ** P.n + P.k_d) + P.alpha * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1) * (P.gamma * mu1 ** P.n * P.n / mu1 / (mu1 ** P.n + P.k_d) - P.gamma * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1) - (F - P.alpha * mu1 ** P.n / (mu1 ** P.n + P.k_d) - P.beta) / (P.gamma * mu1 ** P.n / (mu1 ** P.n + P.k_d) + P.zeta) * (-P.alpha * mu1 ** P.n * P.n ** 2 / mu1 ** 2 / (mu1 ** P.n + P.k_d) + P.alpha * mu1 ** P.n * P.n / mu1 ** 2 / (mu1 ** P.n + P.k_d) + 3 * P.alpha * (mu1 ** P.n) ** 2 * P.n ** 2 / mu1 ** 2 / (mu1 ** P.n + P.k_d) ** 2 - 2 * P.alpha * (mu1 ** P.n) ** 3 / (mu1 ** P.n + P.k_d) ** 3 * P.n ** 2 / mu1 ** 2 - P.alpha * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1 ** 2) - (F - P.alpha * mu1 ** P.n / (mu1 ** P.n + P.k_d) - P.beta) ** 2 / (P.gamma * mu1 ** P.n / (mu1 ** P.n + P.k_d) + P.zeta) ** 3 * (P.gamma * mu1 ** P.n * P.n / mu1 / (mu1 ** P.n + P.k_d) - P.gamma * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1) ** 2 + 1 / 2 * (F - P.alpha * mu1 ** P.n / (mu1 ** P.n + P.k_d) - P.beta) ** 2 / (P.gamma * mu1 ** P.n / (mu1 ** P.n + P.k_d) + P.zeta) ** 2 * (P.gamma * mu1 ** P.n * P.n ** 2 / mu1 ** 2 / (mu1 ** P.n + P.k_d) - P.gamma * mu1 ** P.n * P.n / mu1 ** 2 / (mu1 ** P.n + P.k_d) - 3 * P.gamma * (mu1 ** P.n) ** 2 * P.n ** 2 / mu1 ** 2 / (mu1 ** P.n + P.k_d) ** 2 + 2 * P.gamma * (mu1 ** P.n) ** 3 / (mu1 ** P.n + P.k_d) ** 3 * P.n ** 2 / mu1 ** 2 + P.gamma * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1 ** 2) - 1 / 2 * (P.gamma * mu1 ** P.n * P.n ** 2 / mu1 ** 2 / (mu1 ** P.n + P.k_d) - P.gamma * mu1 ** P.n * P.n / mu1 ** 2 / (mu1 ** P.n + P.k_d) - 3 * P.gamma * (mu1 ** P.n) ** 2 * P.n ** 2 / mu1 ** 2 / (mu1 ** P.n + P.k_d) ** 2 + 2 * P.gamma * (mu1 ** P.n) ** 3 / (mu1 ** P.n + P.k_d) ** 3 * P.n ** 2 / mu1 ** 2 + P.gamma * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1 ** 2) / (P.gamma * mu1 ** P.n / (mu1 ** P.n + P.k_d) + P.zeta) + 1 / 2 * (P.gamma * mu1 ** P.n * P.n / mu1 / (mu1 ** P.n + P.k_d) - P.gamma * (mu1 ** P.n) ** 2 / (mu1 ** P.n + P.k_d) ** 2 * P.n / mu1) ** 2 / (P.gamma * mu1 ** P.n / (mu1 ** P.n + P.k_d) + P.zeta) ** 2)
            #sig1=-1/(-(-P.alpha*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)+P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)+2*(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2*(-P.alpha*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)+P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)-(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)*(-P.alpha*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)+P.alpha*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)+3*P.alpha*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2-2*P.alpha*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2-P.alpha*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)-(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^3*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2+1/2*(F-P.alpha*mu1^P.n/(mu1^P.n+P.k_d)-P.beta)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2*(P.gamma*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)-P.gamma*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)-3*P.gamma*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2+2*P.gamma*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2+P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)-1/2*(P.gamma*mu1^P.n*P.n^2/mu1^2/(mu1^P.n+P.k_d)-P.gamma*mu1^P.n*P.n/mu1^2/(mu1^P.n+P.k_d)-3*P.gamma*(mu1^P.n)^2*P.n^2/mu1^2/(mu1^P.n+P.k_d)^2+2*P.gamma*(mu1^P.n)^3/(mu1^P.n+P.k_d)^3*P.n^2/mu1^2+P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1^2)/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)+1/2*(P.gamma*mu1^P.n*P.n/mu1/(mu1^P.n+P.k_d)-P.gamma*(mu1^P.n)^2/(mu1^P.n+P.k_d)^2*P.n/mu1)^2/(P.gamma*mu1^P.n/(mu1^P.n+P.k_d)+P.zeta)^2);
        else:
            self.mu1 = 5e-5
            self.sig1 = 5e-5
            
    def update_moments(self, A, states, t):
        '''
        @param A: A Memoized instance 
        @param states: States instance
        @param t: current frame or timestep  
        '''
        S = states #give me convenience or give me death!
        P = self.P
        V = self.V
        
        self.init_lik()
        self.p[0] = 1.0
        #if we had called init_like we'd now be propagating those vals into self.[mu,sig2]
        
        #two blocks for spike histories
        
        if(t < (V.T-1)):
            phat  = 1-numpy.exp(-numpy.exp(P.kx[t+1])*V.dt)  #if k, kx are non-scalar, then the inner term becomes -numpy.exp( p.kx[t+1).T * V.dt
        
        #an intermittent sampling for loop: tt=s:-1:2
        
        #ok but we must need this code!
        
#        for tt=s:-1:2
#    O.p_o(1:2^(s-tt+1),tt-1)    = repmat(O.p_o(1:2^(s-tt),tt),2,1).*[(1-phat(tt))*ones(1,2^(s-tt)) phat(tt)*ones(1,2^(s-tt))]';
#    O.mu_o(1:2^(s-tt+1),tt-1)   = (1-P.a)^(-1)*(repmat(O.mu_o(1:2^(s-tt),tt),2,1)-P.A*A.spikemat(1:2^(s-tt+1),tt-1)-P.a*P.C_0);     %mean of P[O_s | C_k]
#    O.sig2_o(tt-1)              = (1-P.a)^(-2)*(P.sig2_c+O.sig2_o(tt)); % var of P[O_s | C_k]
#
#    for n=0:s-tt+1
#        nind=A.ninds{n+1};
#        O.p(n+1,tt-1)   = sum(O.p_o(nind,tt-1));
#        ps          = (O.p_o(nind,tt-1)/O.p(n+1,tt-1))';
#        O.mu(n+1,tt-1)  = ps*O.mu_o(nind,tt-1);
#        O.sig2(n+1,tt-1)= O.sig2_o(tt-1) + ps*(O.mu_o(nind,tt-1)-repmat(O.mu(n+1,tt-1)',A.lenn(n+1),1)).^2;
#    end
#end
        
        #if s== 2 , another intermittent sampling block 
        
        #the while loop to get rid of NaN : 
        #python/numpy doesn't seem to allow assignment of NaN, so we'll have to watch for div by zero errors as we go.
        
        self.p += numpy.finfo(float).eps #to avoid div by zeros
        
        
         



class States(object):
    """ states of the model"""
    def __init__(self, vars, pars):
        '''
        set up the states. 
        @param vars: instance of a Variables object
        @param pars: instance of a Parameters object  
        '''
        self.V = vars
        self.P = pars
        #convenience:
        V = self.V
        P = self.P
        
        self.p = (P.k/V.dt) + numpy.zeros((vars.Nparticles, vars.T)) #rate
        n = numpy.zeros((V.Nparticles, V.T))
        self.n = n.astype('bool')     #spike counts
        self.C = P.C_init * numpy.ones((V.Nparticles, V.T)) #calcium -- probably to be rao-blackwellized away tho!
        self.w_f = (1.0 / V.Nparticles) * numpy.ones((V.Nparticles, V.T)) #forward particle weights
        self.w_b = (1.0 / V.Nparticles) * numpy.ones((V.Nparticles, V.T)) #backward particle weights
        
        #note: i think we shouldn't need this, or it shouldn't be a function of T_o, which we've gotten rid of. 
        # but i want it here commented out so that when i try to use S.Neff i remember why it doesn't exist.
        self.Neff = (1.0 / V.Nparticles) * numpy.ones((1,V.T+1))  
        
        
        #just set S.p to the appropriate constant of putative firing rate times dt. 
        #S.p = repmat(1-exp(-exp(P.kx)*V.dt)',1,V.Nparticles)';


    def prior_sampler(self, memoized, t):
        '''
        @param memoized: instance of Memoized
        @param t: current frame/timestep.  
        '''
        
        #convenience! :
        A = memoized
        P = self.P
        V = self.V
        F = V.F
        
        #spike histories block
        
        #self.p_new = self.p[:,t]
        
        self.n[:,t] = A.U_sampl[:,t] < self.p[:,t]
        
        #this line is clearly a calcium thing but it'll have to be something for the hill stuff, which is worthwhile since the dye definitely saturates. 
        self.C[:,t]        = (1-P.a)*self.C[:,t-1]+P.A*self.n[:,t]+P.a*P.C_0+A.epsilon_c[:,t]
        
        #then there's an if for intermittent sampling that is now always true
        S_mu = Hill_v1(P,self.C[:,t])
        F_mu = P.alpha*S_mu+P.beta   #E[F_t]
        F_var = P.gamma*S_mu+P.zeta  #V[F_t]
        ln_w = -0.5* numpy.power((F[t] - F_mu),2.) / F_var - numpy.log(F_var)/2.
        
        for i in xrange(len(ln_w)):
            if(numpy.isnan(ln_w[i])):
                print('particle %d ln_w is NaN. \n F[t]=%f \n n=%s \n C=%f \n S_mu=%f \n F_mu=%f \n F_var=%f' %
                (i,F[t], self.n[i,t], self.C[i,t], S_mu[i], F_mu[i], F_var[i]))
        
        
        if(t>0):
            ln_w += numpy.log(self.w_f[:,t-1])
            
        #print(ln_w)
        
        ln_w = ln_w - numpy.max(ln_w)

        
        
        w = numpy.exp(ln_w)

        self.w_f[:,t] = w/numpy.sum(w)
        #print(self.w_f[:,t])
        
        #print('gamma: %f    zeta: %f'%(P.gamma, P.zeta))
        #print(S_mu)
        #print(w)
    
    
def backward(vars, pars):
    '''
    '''
      


class Z(object):
    '''more memoized shit'''
    def __init__(vars, states):
        self.oney = numpy.ones(size=(vars.Nparticles,1))
        self.zeroy = numpy.zeros(size = (vars.Nparticles,vars.Nparticles))
        self.C0 = states.C[:,vars.T-1]
        self.C0mat = numpy.zeros(shape=(len(oney),len(C0)))
        for i in xrange(len(c0)):
            c0mat[:,i] = c0[i]
        
        
    

def forward(vars, pars):
    
    '''
    the model is F_t = f(C) = alpha C^n/(C^n + k_d) + beta + e_t,
    where e_t ~ N[0, gamma*f(C)+zeta]
    
    @param vars: instance of a variables object 
    @param pars: instance of a parameters object 
    
    @return: instance of a States object  (simulation states)
    '''
    A = Memoized(vars, pars) 
    S = States(vars, pars)
    #skipping the spike history stuff, but it would go roughly here-ish, and in some __init__s. 
    #O = ObsLik(vars, pars)
    
    #O.p[0] = 1.
    #O.mu[0] = O.mu_o[0]
    #O.sig2[0] = O.sig2_o[0]
    
    #skipping another V.freq block
    
    #O.update_moments(A,S,0) # the 0 used to be s, which is initialized to V.Freq, which is 1. but is it even really a constant, or does V.freq incremement somewhwere?
    
    #convenience:
    V = vars
    P = pars
    Nresamp = 0
    #here is the particle filter:
    for t in xrange(1,V.T): #are these the right timestep bounds?
        print(t)
        S.prior_sampler(A,t)

        #spikeHist block
        
        #here is stratified respampling:
        Nresamp = t
        S.Neff[0,Nresamp] =  1/numpy.sum(numpy.power(S.w_f[:,t],2))
        #print(S.Neff[0,Nresamp])
        if(S.Neff[0,Nresamp] < V.Nparticles/2.0):
            #so resample:
            print('resampling')
            edges = numpy.insert(S.w_f[:,t].cumsum(),0,0)
            ind = histc_j(A.U_resamp[Nresamp,:], edges) # this does the strat. resamp.
            random.shuffle(ind) #do a permutation of the inds (to avoid potential biases.?)
            
            S.p[:,t]   =  S.p[ind,t]   #:: if V.freq is 1, then t-V.freq+1 = t. right? #t-V.freq+1:t];      #% resample probabilities (necessary?)
            S.n[:,t]   =  S.n[ind,t]
            S.C[:,t]   =  S.C[ind,t]
            S.w_f[:,t] = (1.0/V.Nparticles)*numpy.ones((V.Nparticles)) #% reset weights
            #Nresamp += 1
            
        #skipping a spikehist block
        #O.update_moments(A,S,t);                      #% estimate P[O_s | C_tt] for all t'<tt<s as a gaussian
        #O.s = t #update the time var
    
        
    return S
        
        
    
def Hill_v1(pars,C):
    '''
    % generalized hill model
    '''
    C[C<0]  = 5e-5;
    return numpy.power(C,pars.n) / ( numpy.power(C,pars.n) +pars.k_d)

def histc_j(x, edges):
    '''
    given a vector, x, and a (sorted) list of N+1 edges defining the N bins, 
    returns a map, m, s.t. m[j] tells you which bin x[j] falls into.
    '''
    inds = numpy.zeros(x.size, dtype=numpy.int)
    for i in xrange(len(x)):
        inds[i] = int(max(0, bisect.bisect_left(edges,x[i])-1) )
        #bisect is intended for mainting a sorted list, so it only returns 0
        #if x[i]<=edges[0],  so the call to max is to ensure that we don't put
        #the smallest value into some imaginary bin -1. 
    
    return inds
        
    

def backward(vars, states, pars):
    z = Z(vars, states)
    for t in reversed(xrange(vars.T)):
        z = step_backward(vars, states, pars, z, t)
        states.w_b[:,t-1] = z.w_b
        
        
def step_backward(vars, states, pars, z, t):
    
    # % compute ln P[n_t^i | h_t^i]
    z.n1 = states.n[:,t]
    ln_Pn = numpy.zeros(shape = z.oney.shape)
    spikingInds = z.n1==1
    nsInds = z.n1 != 1
    ln_Pn[spikingInds] = numpy.log(states.p[spikingInds,t])
    ln_Pn[nsInds]  = numpy.log(1.0 - states.p[nsInds,t])
    
    # % compute ln P[C_t^i | C_{t-1}^j, n_t^i]
    z.C0 = S.c[:,t-1]
    z.C1 = S.c[:,t]
    z.C0mat = numpy.zeros(shape=(len(z.oney),len(z.C0)))
    for i in xrange(len(z.C0)):
        z.C0mat[:,i] = z.C0[i]    
    z.C1mat = numpy.zeros(shape=(len(z.oney),len(z.C1)))
    for i in xrange(len(z.C1)):
        z.C1mat[:,i] = z.C1[i]    
    mu = (1-pars.a) * states.C[:,t-1] + pars.a * z.n1 + pars.a*pars.C_0
    mumat = numpy.zeros(shape = (len(z.oney),len(mu)))
    for i in xrange(len(mu)):
        mumat[:,i] = mu[i]
    ln_PC_Cn = -0.5 * (z.C1mat - mumat)**2 / pars.sig2_c
    
    #% compute ln P[h_t^i | h_{t-1}^j, n_{t-1}^i]
    ln_Ph_hn = z.zeroy
    #spike hist block.. but i think we skip it. 
    
    #% compute P[H_t^i | H_{t-1}^j]
    sum_lns = ln_PC_Cn + ln_Ph_hn
    for i in xrange(len(ln_Pn)):
        sum_lns[:,i] += ln_Pn[i]
    mx = numpy.max(sum_lns, 1) #seems that here, 1 is the rows? 
    T = numpy.exp(sum_lns - mx)
    ncols = T.shape[1] #we want each column to sum to one (here, 1 is cols)
    for i in xrange(ncols):
        normer = numpy.sum(T[:,i])
        for j in xrange(T.shape[0]):
            T[i,j] /= normer
    
    #% compute P[H_t^i, H_{t-1}^j | 0]
    
    
      
    
    
def realDataTest():
    par = os.path.pardir
    sep = os.path.sep
    fstring=par + sep + par + sep + 'data'+sep+'fluo.txt'
    #print(fstring)
    fluofile = open(fstring,'r')
     
    line = fluofile.readline()
    #print(line)
    numtext=re.split('[\t]', line)
    numlist=numpy.zeros(len(numtext)-1)  #omitting the last value in the fluo file b/c it's a zero-pad. 
    for i in xrange(len(numtext)-1):
        #print('%d:  %s'%(i, numtext[i]))
        numlist[i] = float(numtext[i])
        
    fluofile.close()
    
    numlist -= numpy.min(numlist)
    numlist /= numpy.max(numlist)
    
    v = Variables(numlist, 0.075)
    p = Parameters(v, A=20, C_0=numpy.mean(numlist[1:10]),C_init=numpy.mean(numlist[1:10]))
    S = forward(v,p)
    
    cbar = numpy.zeros(p.V.T)
    nbar = numpy.zeros(p.V.T)
    
    for t in xrange(p.V.T):
        for i in xrange(p.V.Nparticles):
            weight = S.w_f[i,t]
            cbar[t] += weight * S.C[i,t]
            nbar[t] += weight * S.n[i,t]
            
        
        
    #nbar /= p.V.Nparticles
    pylab.figure()
    pylab.hold(True)
    pylab.plot(nbar, label='expected spikes')
    pylab.title('expected spikes')
    pylab.plot(numlist, label='F')
    pylab.legend()
    
    pylab.figure()
    pylab.plot(cbar)
    pylab.title('expected Ca')
    
    #pylab.plot(numlist)
    pylab.show()
    
if __name__ == "__main__":
    realDataTest()