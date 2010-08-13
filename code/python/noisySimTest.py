import smc
import numpy, pylab




def setupSimData(spt = [26, 50, 84, 128, 199, 247, 355] ):
    T = 405     # no. of time steps
    dt = .025   # time step size (seconds)
    Nsec = T*dt #sim time in seconds
    tvec = numpy.arange(Nsec, step=dt) 
    
    # spike times   [50,128,247]#
    rate  = len(spt) / Nsec
    ### A = 1 # jump size
    tau = 0.5 # Ca decay const
    lam = 1. / numpy.sqrt(rate)
    sig = 1.
    
    k = numpy.log(-numpy.log(1-rate*dt)/dt)
    tau_c = tau
    A = 10 # jump size
    C_0 = .1
    C_init = C_0
    sigma_c = sig
    n = 1.
    k_d = 200.0
    alpha = 1
    beta = 0
    gamma = 0.001
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

    eps_t = numpy.random.normal(size=(T,1)) / 10.0
    F = alpha * S + beta + numpy.sqrt(gamma*S+zeta)*eps_t
    F[F<0] = numpy.finfo(float).eps
    
    V = smc.Variables(F,dt, true_n=spt, Nparticles=99)
    P = smc.Parameters(V, tau_c=tau_c, A=A, C_0 = C_0, C_init=C_0, sigma_c = sigma_c, k_d=k_d,
                        alpha=alpha, beta=beta,gamma=gamma, zeta=zeta)
    
    return P

    

def forwardTest():
    spikeTimes = spt = [26, 50, 84, 128, 199, 247, 355] 
    P = setupSimData(spt = spikeTimes )
    S = smc.forward(P.V, P)
    
    pylab.figure()
    pylab.plot(P.V.F)
    pylab.title('F')
    

    cbar = numpy.zeros(P.V.T)
    nbar = numpy.zeros(P.V.T)
    
    for t in xrange(P.V.T):
        for i in xrange(P.V.Nparticles):
            weight = S.w_f[i,t]
            cbar[t] += weight * S.C[i,t]
            nbar[t] += weight * S.n[i,t]
    
    #cbar /= P.V.Nparticles
    #nbar /= P.V.Nparticles
    
    pylab.figure()
    pylab.hold(True)
    pylab.title('expected Ca vs arbitrary particles')
    pylab.plot(S.C[3,:], label='particle 3')
    pylab.plot(S.C[10,:], label='particle 10')
    pylab.plot(S.C[20,:], label='particle 20')
    pylab.plot(cbar, label='E[C]')
    pylab.plot(spikeTimes, 6+numpy.ones(len(spikeTimes)), 'k.', label='simulated spike times')

    pylab.legend()
    
    
    pylab.figure()
    pylab.plot(nbar)
    pylab.title('expected spikes')
    

    #pylab.figure()
    pylab.matshow(S.w_f)
    pylab.title('min s.w_f %f , max s.w_f: %f'%(numpy.min(S.w_f), numpy.max(S.w_f)))
    
    
    pylab.figure()
    pylab.hold(True)
    pylab.plot(S.w_f[3,:], label='particle 3')
    pylab.plot(S.w_f[10,:], label='particle 10')
    pylab.plot(S.w_f[20,:], label='particle 20')
    pylab.legend()
    pylab.title('individual particle weights')
    
    pylab.show()
    
    



if __name__ == "__main__":
    forwardTest()


