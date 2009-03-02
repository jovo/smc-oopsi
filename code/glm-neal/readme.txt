This contains few matlab routines for HMM-based analysis of spikes.

* runSMC is a sample script showing parameters assignment, generating samples of observation sequences from HMM, infering state sequence from observations and evaluating inferrence accuracy.
* errorsSMC is spike-train comparison routine using bipartite matching of target and inferred spike trains
* prepMC is a sampler for sequence of states/observations from given HMM with either continuous or discrete state.
* sampleSMC is sampler from state sequences given observation sequence in continuous state HMM following Neal's stochastic grids (pulls) algorithm.
* sampleMC is forward-backward procedure for generating state sequence given observation sequence in discrete state HMM.
* filtMC is forward-backward procedure for generating P(x_t|Y) given observation sequence in discrete state HMM.
* spPXX is an example of P(X_1|X_0,t) definition for an HMM for use with sampleSMC/sampleMC/filtMC.
* spPY is an example of P(Y|X,t) definition for an HMM for use with sampleSMC/sampleMC/filtMC.
* spPG is an example of P(x|t) definition for stochastic grids (pulls) in Neal's algorithm.
* _GOOPSI scripts are HMM definition for Flu-[Ca]-spikes system
 
* assignmentoptimal is bipartite matching mex-function (due Markus Buehren, available at matlab central exchange).