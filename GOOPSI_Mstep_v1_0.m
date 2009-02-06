function Enew = GOOPSI_Mstep_v1_0(Sim,S,M,E,F)
% this function finds the mle of the parameters
%
% Input---
% Sim:  simulation parameters
% R:    real data
% S:    simulation results
% M:    moments and sufficient stats
% E:    old parameter estimates
%
% Output is 'Enew', a structure with a field for each parameter, plus some
% additional fields for various likelihoods

Enew = E;   % initialize parameters
lik = [];   % initialize likelihood

if Sim.n_params == true
    % MLE for spike rate parameters: baseline (b), linear filter (k), and spike history weights (omega)
    fprintf('\nestimating spike rate params\n')
    RateParams=E.k;                                                 % vector of parameters to estimate (changes depending on user input of which parameters to estimate)
    sp      = S.n==1;                                               % find (particles,time step) pairs that spike
    nosp    = S.n==0;                                               % don't spike
    x       = repmat(Sim.x,1,Sim.N);                                % generate matrix for gradinent
    zeroy   = zeros(Sim.N,Sim.T);                                   % make matrix of zeros for evaluating lik

    if Sim.h_params == true
        if Sim.M>0                                                  % if spike history terms are present
            RateParams=[RateParams; E.omega];                       % also estimate omega
            for i=1:Sim.M                                           % and modify stimulus matrix for gradient
                x(Sim.StimDim+i,:)=reshape(S.h(:,:,i),1,Sim.N*Sim.T);
            end
        end

        options     = optimset('Display','off','GradObj','off');    % use gradient
        [bko lik_r] = fminunc(@f_bko,RateParams,options);           % find MLE
        Enew.k      = bko(1:end-Sim.M);                             % set new parameter estimes
        if Sim.M>0                                                  % for omega too
            Enew.omega = bko(end-Sim.M+1:end);
        end

    else
        if Sim.M>0                                                  % if spike history terms are present
            for i=1:Sim.M                                           % and modify stimulus matrix for gradient
                x(Sim.StimDim+i,:)=reshape(S.h(:,:,i),1,Sim.N*Sim.T);
            end
        end

        options     = optimset('Display','off','GradObj','off');    % use gradient
        [bk lik_r]  = fminunc(@f_bk,RateParams,options);            % find MLE
        Enew.k      = bk(1:end);                                    % set new parameter estimes
    end
    Enew.lik_r   = -lik_r;
    lik = [lik Enew.lik_r];
end

    function [lik dlik]= f_bko(RateParams)                          % get lik and grad

        xk      = RateParams(1:end-Sim.M)'*Sim.x;                   % filtered stimulus
        hs      = zeroy;                                            % incorporate spike history terms
        for l=1:Sim.M
            hs  = hs+RateParams(end-Sim.M+l)*S.h(:,:,l);
        end
        s       = repmat(xk,Sim.N,1) + hs;

        f_kdt   = exp(s)*Sim.dt;                                    % shorthand
        lik     = -sum(S.w_b(sp).*log(1-exp(-f_kdt(sp))))...        % liklihood
            +sum(S.w_b(nosp).*f_kdt(nosp));

        if nargout > 1                                              % if gradobj=on
            ef      = exp(f_kdt);                                   % shorthand
            dlik      = x(:,sp)*(S.w_b(sp).*f_kdt(sp)./( ef(sp)-1))... %gradient of lik
                -x(:,nosp)*(S.w_b(nosp).*f_kdt(nosp));
        end
    end %function f_bko


    function [lik dlik]= f_bk(RateParams)                           % get lik and grad

        xk      = RateParams'*Sim.x;                                % filtered stimulus
        hs      = zeroy;                                            % incorporate spike history terms
        for l=1:Sim.M
            hs  = hs+E.omega*S.h(:,:,l);
        end
        s       = repmat(xk,Sim.N,1) + hs;

        f_kdt   = exp(s)*Sim.dt;                                    % shorthand
        lik     = -sum(S.w_b(sp).*log(1-exp(-f_kdt(sp))))...        % liklihood
            +sum(S.w_b(nosp).*f_kdt(nosp));

        if nargout > 1                                              % if gradobj=on
            ef      = exp(f_kdt);                                   % shorthand
            dlik      = x(:,sp)*(S.w_b(sp).*f_kdt(sp)./( ef(sp)-1))... % gradient of lik
                -x(:,nosp)*(S.w_b(nosp).*f_kdt(nosp));
        end
    end %function f_bko

%% MLE for calcium parameters
if Sim.C_params == true
    fprintf('estimating calcium parammeters\n')
    [ve_x fval] = quadprog(M.Q, M.L,[],[],[],[],[0 0 0],[inf inf inf],[1/E.tau_c E.A E.C_0/E.tau_c]+eps);
    Enew.tau_c  = 1/ve_x(1);
    Enew.A      = ve_x(2);
    Enew.C_0    = ve_x(3)/ve_x(1);
    Enew.sigma_c= sqrt(-(0.5*ve_x'*M.Q*ve_x + M.L'*ve_x))/(Sim.T*Sim.dt);
    Enew.lik_c  = -(M.K/2+fval)/(2*Enew.sigma_c*sqrt(Sim.dt)) - M.J*log(Enew.sigma_c); %lik for [Ca]
    lik = [lik Enew.lik_c];
end

% % %% MLE for spike history parameters
% % for m=1:Sim.M
% %     Enew.sigma_h(m)= sum(M.v{m})/Sim.T;
% % end
%
% %% MLE for observation parameters
% O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
% Oind    = find(~isnan(O));
% Enew.sigma_o = sqrt(sum(sum(S.w_b(:,Oind).*(repmat(O(Oind),Sim.N,1)-S.C(:,Oind)).^2)/Sim.T));


%% MLE for observation parameters
if Sim.F_params == true
    fprintf('estimating observation parammeters\n')
    options         = optimset('Display','off','GradObj','off'); %use gradient
    ab_0            = [E.alpha E.beta];
%     [ab Enew.lik_o]   = fminunc(@f_ab,ab_0,options);%find MLE
    [Enew.lik_o ab]   = f1_ab(ab_0);%find MLE for {alpha, beta and gamma*}
    Enew.alpha=ab(1);
    Enew.beta=ab(2);
    Enew.gamma=E.gamma*ab(3);
    Enew.zeta=E.zeta*ab(3);
    lik = [lik Enew.lik_o];
end

%     function [lik x] = f_ab(ab_o)
%         Fmean=Hill_v1(E,S.w_b.*S.C);
%         A=-[Fmean(:) ones(Sim.N*Sim.T,1)];
%         H=A'*A;
%         f=A'*repmat(F,Sim.N,1);
%         [x lik] = quadprog(H,f,[],[],[],[],[0 0],[inf inf]);        
%     end %function f_ab_0
    
  function [lik x] = f1_ab(ab_o)
    %THIS EXPLICITLY ASSUMES WEIGHTS w_b ARE SUM=1 NORMALIZED
    pfS=Hill_v1(E,S.C);
    pfV=E.gamma*pfS+E.zeta;
    % minimize quadratic form of E[(F - ab(1)*pfS - ab(2))^2/pfV]
    % taken as weighted average over all particles (Fmean)
    f1_abH=0; f1_abf=0; f1_abc=0; f1_abn=0;
    for i=1:size(pfS,1)
      f1_abH = f1_abH + ...
        [sum(pfS(i,:).^2./pfV(i,:).*S.w_b(i,:)) sum(pfS(i,:)./pfV(i,:).*S.w_b(i,:));
         sum(pfS(i,:)./pfV(i,:).*S.w_b(i,:))    sum(S.w_b(i,:)./pfV(i,:))];
      f1_abf = f1_abf - ...
        [sum(pfS(i,:).*F(:)'./pfV(i,:).*S.w_b(i,:));
         sum(F(:)'./pfV(i,:).*S.w_b(i,:))];
      f1_abc = f1_abc + sum(F(:)'.^2./pfV(i,:).*S.w_b(i,:));
      f1_abn = f1_abn + sum(S.w_b(i,:));
    end
    % % this is an alternative setup for doing this
    % pfS=sum(S.w_b.*Hill_v1(E,S.C),1);
    % pfV=E.gamma*pfS+E.zeta;
    % A=-[pfS(:)./sqrt(pfV(:)) ones(Sim.T,1)./sqrt(pfV(:))];
    % H=A'*A;
    % f=A'*(F(:)./sqrt(pfV(:)));
    % solve as QP given ab(1)>0 and no bounds on ab(2)
    [x lik] = quadprog(f1_abH,f1_abf,[],[],[],[],[0 -inf],[inf inf]);
    lik=(lik+f1_abc/2);               %estimate the variance
    if(isfield(Sim,'G_params') && Sim.G_params==1)%estimate gamma_new/gamma
      x(3)=lik/f1_abn; 
    else%if estimating gamma doesn't seem to be working out...
      x(3)=1;
    end
    lik=-lik-f1_abn*log(x(3))/2-sum(log(pfV(:)).*S.w_b(:))/2;
  end %function f_ab_0
  %%
Enew.lik=sum(lik);
end
