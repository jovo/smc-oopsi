import smc
import numpy, pylab















def runsim():
    T = 405     # no. of time steps
    dt = .025   # time step size (seconds)
    Nsec = T*dt #sim time in seconds
    tvec = numpy.arange(Nsec, step=dt) 
    N = 200     # no. of particles
    
    spt     = [26 50 84 128 199 247 355] # spike times?
    rate  = len(spt) / Nsec
    ### A = 1 # jump size
    tau = 0.5 # Ca decay const
    lam = 1. / numpy.sqrt(rate)
    sig = 1.
    
    k = numpy.log(-numpy.log(1-rate*dt)/dt)
    tau_c = tau
    A = 5 # jump size
    C_0 = .1
    C_init = C_0
    sigma_c = sig
    n = 1.
    k_d = 200
    alpha = 1
    beta = 0
    gamma = 0.
    zeta = 5e-5
    a = dt/tau_c
    
    nn = numpy.zeros(T)
    nn[spt] = 1
    C = numpy.zeros(T)
    epsilon_c = sigma_c*numpy.sqrt(dt)*numpy.random.normal(T)
    for t in xrange(1,T):
        C[t] = (1-a)*C[t-1] + a*C_0 + A*n[t]+epsilon_c[t]
        
    S = smc.
    
    
    





if __name__ == "__main__":
    runsim()
     


