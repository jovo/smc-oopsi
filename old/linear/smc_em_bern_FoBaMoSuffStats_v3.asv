function [S M] = smc_em_bern_FoBaMoSuffStats_v3(Sim,R,P)

%% makes code prettier :-)
B.sig2_c    = P.sigma_c^2*Sim.dt;
B.sig2_o    = P.sigma_o^2;
B.a         = 1-Sim.dt/P.tau_c;
B.beta      = P.beta;
B.kx        = P.k'*Sim.x;
if Sim.M>0
    B.sig2_h    = P.sigma_h.^2*Sim.dt;
    B.g         = 1-Sim.dt/P.tau_h;
    B.omega     = P.omega;
end

%% forward step
fprintf('\nTotal Number of Steps=%d\n',Sim.T)
if Sim.pf==0
    fprintf('\nprior forward step:          ')
    S = smc_em_bern_PriorSampl_v2(Sim,R,B);
elseif Sim.pf==1    
    fprintf('\nbackward forward step:       ')
    S = smc_em_bern_BackSampl_v7(Sim,R,B);
elseif Sim.pf==2
    fprintf('\n2^k-mixture forward step:    ')
    S = smc_em_bern_RealBackSampl_v3(Sim,R,B);
elseif Sim.pf==3
    fprintf('\nk-mixture forward step:      ')
    S = smc_em_bern_k_comp_BackSampl_v1(Sim,R,B);
end

%% initialize stuff for backwards step and suff stats
S.w_b   = 1/Sim.N*ones(Sim.N,Sim.T);
ln_Pn   = zeros(Sim.N,1);
oney    = ones(Sim.N,1);
zeroy2  = zeros(Sim.N,Sim.N);
M.Q     = zeros(2,2);
M.L     = zeros(2,1);
M.u     = 0;
for m = 1:Sim.M,
    M.v{m} = 0;
end

%% get suff stats
fprintf('\nbackward step..........')
C0 = S.C(:,Sim.T);
C0mat = C0(:,oney)';
for t=Sim.T-Sim.freq-1:-1:Sim.freq

    % compute ln P[n_t^i | h_t^i]
    n1          = S.n(:,t);                             %for prettiness sake
    ln_Pn(n1==1)= log(S.p(n1==1,t));                    %P[n=1] for those that spiked
    ln_Pn(~n1)  = log(1-S.p(~n1,t));                    %P[n=0] for those that did not

    % compute ln P[C_t^i | C_{t-1}^j, n_t^i]
    C0          = S.C(:,t-1);                           %for prettiness sake
    C1          = S.C(:,t);
    C1mat       = C1(:,oney);                           %recall from previous time step
    C0mat       = C0(:,oney);                           %faster than repamt
    mu          = B.a*S.C(:,t-1)+B.beta*n1;             %mean
    mumat       = mu(:,oney)';                          %faster than repmat
    ln_PC_Cn    = -0.5*(C1mat - mumat).^2/B.sig2_c;     %P[C_t^i | C_{t-1}^j, n_t^i]

    % compute ln P[h_t^i | h_{t-1}^j, n_{t-1}^i]
    ln_Ph_hn    = zeroy2;                               %reset transition prob for h terms
    for m=1:Sim.M                                       %for each h term
        h1      = S.h(:,t,m);                           %faster than repmat
        h1      = h1(:,oney);                           
        h0      = B.g(m)*S.h(:,t-1,m)+S.n(:,t-1);       
        h0      = h0(:,oney)';                          
        ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/B.sig2_h(m);
    end

    % compute P[H_t^i | H_{t-1}^j]
    sum_lns = ln_Pn(:,oney)+ln_PC_Cn + ln_Ph_hn;        %in order to ensure this product doesn't have numerical errors
    mx      = max(sum_lns,[],1);                        %find max in each of row
    mx      = mx(oney,:);                               %make a matrix of maxes
    T0      = exp(sum_lns-mx);                          %exponentiate subtracting maxes (so that in each row, the max entry is exp(0)=1
    Tn      = sum(T0,1);                                %then normalize
    T       = T0.*repmat(1./Tn(:)', Sim.N, 1);          %such that each column sums to 1        

    % compute P[H_t^i, H_{t-1}^j | O]
    PHHn    = (T*S.w_f(:,t-1))';                        %denominator
    PHHn    = PHHn+eps;
    PHHn2   = PHHn(oney,:)';                            %faster than repmat                     
    PHH     = T .* (S.w_b(:,t)*S.w_f(:,t-1)')./PHHn2;   %normalize such that sum(PHH)=1

    S.w_b(:,t-1)= sum(PHH,1);                           %marginalize to get P[H_t^i | O]

    C0dt    = C0*Sim.dt; 
    bmat    = C1mat-C0mat';
    bPHH    = PHH.*bmat;

    M.Q(1,1)= M.Q(1,1) + C0dt'*PHH*C0dt; 
    M.Q(1,2)= M.Q(1,2) - n1'*PHH*C0dt; 
    M.Q(2,2)= M.Q(2,2) + sum(PHH'*(n1.^2));

    M.L(1)  = M.L(1) + sum(bPHH*C0dt); 
    M.L(2)  = M.L(2) - sum(bPHH'*n1);


    if mod(t,100)==0 && t>=9900
        fprintf('\b\b\b\b\b%d',t)
    elseif mod(t,100)==0 && t>=900
        fprintf('\b\b\b\b%d',t)
    elseif mod(t,100)==0
        fprintf('\b\b\b%d',t)
    end
end
M.Q(2,1)= M.Q(1,2);
M.L     = 2*M.L;

%%   get moments
M.nbar  = sum(S.w_b.*S.n,1);
M.nvar  = sum((repmat(M.nbar,Sim.N,1)-S.n).^2)/Sim.N;

M.Cbar  = sum(S.w_b.*S.C,1);
M.Cvar  = sum((repmat(M.Cbar,Sim.N,1)-S.C).^2)/Sim.N;