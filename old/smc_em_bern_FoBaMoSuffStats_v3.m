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

    ln_Pn(S.n(:,t)==1)  = log(S.p(S.n(:,t)==1,t));
    ln_Pn(~S.n(:,t))    = log(1-S.p(~S.n(:,t),t));

    C1          = S.C(:,t);
    C1          = C1(:,oney);
    C0          = B.a*S.C(:,t-1)+B.beta*S.n(:,t);
    C0dt        = C0*Sim.dt; %n01 = -n1;%PHHT  = PHH'; %Xterms = -n1*S.C(:,t-1)'*Sim.dt;  %bmatT = bmat';
    C0          = C0(:,oney)';
    ln_PC_Cn    = -0.5*(C0 - C1).^2/B.sig2_c;

    ln_Ph_hn = zeroy2;
    for m=1:Sim.M
        h1 = S.h(:,t,m);
        h1 = h1(:,oney);
        h0 = B.g(m)*S.h(:,t-1,m)+S.n(:,t-1);
        h0 = h0(:,oney)';
        ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/B.sig2_h(m);
    end

    sum_lns = ln_Pn(:,oney)+ln_PC_Cn + ln_Ph_hn;
    mx      = max(sum_lns,[],1);
    mx      = mx(oney,:);
    T0      = exp(sum_lns-mx);
    Tn      = sum(T0,1);
    T       = T0./Tn(oney,:);

    PHHn    = (T*S.w_f(:,t-1))';
    PHHn    = PHHn+eps;
    PHHn2   = PHHn(oney,:)';
    PHH     =  T .* (S.w_b(:,t)*S.w_f(:,t-1)')./PHHn2;

    S.w_b(:,t-1)= sum(PHH,1);
    
    bmat    = C1-C0';
    bPHH    = PHH.*bmat;

    M.Q(1,1)= M.Q(1,1) + C0dt'*PHH*C0dt; %sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1));
    M.Q(1,2)= M.Q(1,2) - n1'*PHH*C0dt; %sum(PHH(:).*Xterms(:));
    M.Q(2,2)= M.Q(2,2) + sum(PHH'*(n1.^2));%n1'*PHH*n1; %sum(PHH(:).*repmat((-n1).^2,Sim.N,1));

    M.L(1)  = M.L(1) + sum(bPHH*C0dt); %sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)));
    M.L(2)  = M.L(2) - sum(bPHH'*n1); %sum(PHH(:).*(2*repmat(-n1,Sim.N,1).*bmat(:)));


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