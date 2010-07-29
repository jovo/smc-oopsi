import smc
import numpy, pylab




def setupSimData():
    T = 405     # no. of time steps
    dt = .025   # time step size (seconds)
    Nsec = T*dt #sim time in seconds
    tvec = numpy.arange(Nsec, step=dt) 
    
    spt     = [26, 50, 84, 128, 199, 247, 355] # spike times?
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
    k_d = 200.0
    alpha = 1
    beta = 0
    gamma = 0.
    zeta = 5e-5
    a = dt/tau_c
    nn = numpy.zeros(T)
    nn[spt] = 1
    C = numpy.zeros((T,1))
    epsilon_c = sigma_c*numpy.sqrt(dt)*numpy.random.normal(size=(T,1))
    #print(epsilon_c)
    #print(C)
    for t in xrange(1,T):
        C[t] = (1-a)*C[t-1] + a*C_0 + A*nn[t]+epsilon_c[t]

    #uhh not the greatest, but this is just the hill equation from smc.hill_v1.
    #not calling it b/c we don't have a Pars object, b/c Pars requires a Vars, 
    #and Vars requires a timeseries, and the timeseries is what we're building right now !! 
    S = numpy.power(C,n) / ( numpy.power(C,n) + k_d)

    eps_t = numpy.random.normal(size=(T,1))
    F = alpha * S + beta + numpy.sqrt(gamma*S+zeta)*eps_t
    F[F<0] = numpy.finfo(float).eps
    
    V = smc.Variables(F,dt, true_n=spt)
    P = smc.Parameters(V, tau_c=tau_c, A=A, C_0 = C_0, C_init=C_0, sigma_c = sigma_c, k_d=k_d,
                        alpha=alpha, beta=beta,gamma=gamma)
    
    return P

    





if __name__ == "__main__":
    p = setupSimData()
     


