import smc
import numpy, pylab




def setupSimData(spt = [26, 50, 84, 128, 199, 247, 355], A=10, tau=0.5, beta=0, gamma=0e-5, zeta=5e-5, alpha=1 ):
    T = 405     # no. of time steps
    dt = .025   # time step size (seconds)
    Nsec = T*dt #sim time in seconds
    tvec = numpy.arange(Nsec, step=dt) 
    
    # spike times   [50,128,247]#
    rate  = len(spt) / Nsec
    ### A = 1 # jump size
    lam = 1. / numpy.sqrt(rate)
    sig = 1.
    
    k = numpy.log(-numpy.log(1-rate*dt)/dt)
    tau_c = tau
    C_0 = .1
    C_init = C_0
    sigma_c = sig
    n = 1.
    k_d = 200.0
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
    
    pylab.figure()
    pylab.plot(F)
    pylab.show()
    
    V = smc.Variables(F,dt, true_n=spt, Nparticles=99)
    P = smc.Parameters(V, tau_c=tau_c, A=A, C_0 = C_0, C_init=C_0, sigma_c = sigma_c, k_d=k_d,
                        alpha=alpha, beta=beta,gamma=gamma, zeta=zeta)
    
    return P



def paramWalkHelper(spikeTimes, P):
    S = smc.forward(P.V, P)
    trueSpikes = numpy.zeros(P.V.T,dtype=bool)
    trueSpikes[spikeTimes] = True
    
    posterior = 0.0
    for t in xrange(P.V.T):
        nbar = 0.0
        for i in xrange(P.V.Nparticles):
            nbar += S.w_f[i,t] * S.n[i,t]
        if(trueSpikes[t]):
            posterior += nbar
        else:
            posterior += (1-nbar)
            
    print('posterior: %f     post/V.T: %f'%(posterior, posterior/P.V.T))
    
    return posterior/P.V.T
    

def forwardParamWalk():
    spikeTimes = [26,50,128,199,247,355]
    AVals = numpy.array((4,4))#((4.,2.,1.,0.5, 0.25,0.0125))
    tauVals = numpy.arange(0.3,0.9,0.05)
    betaVals = numpy.arange(0,10,1)
    gammaVals = numpy.array((.01,.1))#((0.001,.001,.001,.001,.001, .001))
    zetaVals = numpy.arange(0.00001,0.500001,0.025)
    alphaVals = numpy.arange(0.2,2,0.25)
    
    posteriors = []
    for i in xrange(len(AVals)):
        A = AVals[i]
        gamma = gammaVals[i]
        P = setupSimData(spt=spikeTimes,A=A,gamma=gamma)
        posteriors.append(paramWalkHelper(spikeTimes, P))
    
    print(posteriors)
    pylab.figure()
    pylab.plot(AVals/gammaVals,posteriors,'r.')
    pylab.title('posterior weight as a fn of A/gamma')
    pylab.savefig('walk_A_gamma.png')
    
    return

    posteriors=[]
    for tau in tauVals:
        P = setupSimData(spt=spikeTimes,tau=tau)
        posteriors.append(paramWalkHelper(spikeTimes, P))
    pylab.figure()
    pylab.plot(tauVals,posteriors,'r.')
    pylab.title('posterior weight as a fn of tau')
    pylab.savefig('walk_tau.png')

    posteriors=[]
    for beta in betaVals:
        P = setupSimData(spt=spikeTimes,beta=beta)
        posteriors.append(paramWalkHelper(spikeTimes, P))
    pylab.figure()
    pylab.plot(betaVals,posteriors,'r.')
    pylab.title('posterior weight as a fn of beta')
    pylab.savefig('walk_beta.png')
    
    posteriors=[]
    for zeta in zetaVals:
        P = setupSimData(spt=spikeTimes,zeta=zeta)
        posteriors.append(paramWalkHelper(spikeTimes, P))
    pylab.figure()
    pylab.plot(zetaVals,posteriors,'r.')
    pylab.title('posterior weight as a fn of zeta')
    pylab.savefig('walk_zeta.png')
    
    posteriors=[]
    for gamma in gammaVals:
        P = setupSimData(spt=spikeTimes,gamma=gamma)
        posteriors.append(paramWalkHelper(spikeTimes, P))
    pylab.figure()
    pylab.plot(gammaVals,posteriors,'r.')
    pylab.title('posterior weight as a fn of gamma')
    pylab.savefig('walk_gamma.png')    
    
    posteriors=[]
    for alpha in alphaVals:
        P = setupSimData(spt=spikeTimes,alpha=alpha)
        posteriors.append(paramWalkHelper(spikeTimes, P))
    pylab.figure()
    pylab.plot(alphaVals,posteriors,'r.')
    pylab.title('posterior weight as a fn of alpha')
    pylab.savefig('walk_alpha.png')
    
    pylab.show()
    


        

    

def forwardTest():
    spikeTimes = spt = [26, 50, 84, 128, 199, 247, 355] 
    P = setupSimData(spt = spikeTimes )
    S = smc.forward(P.V, P)
    
    pylab.figure()
    pylab.plot(P.V.F)
    pylab.title('Simulated Fluorescence')
    

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
    
    #print(nbar.shape)
    pylab.figure()
    pylab.hold(True)
    pylab.plot(nbar, label='expected spikes')
    pylab.plot(spikeTimes, 0.5*numpy.ones(len(spikeTimes)), 'k.', label='simulated spike times')
    pylab.legend()
    pylab.title('spike detection')
    

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
    forwardParamWalk()


