function Enew = smc_em_bern_mstep_v1(Sim,R,S,M,E)

Enew = E;

%% MLE for spike rate params
fprintf('\nestimating spike rate params\n')
RateParams=E.k;
sp      = S.n==1;
nosp    = S.n==0;
x       = repmat(Sim.x,1,Sim.N);

options = optimset('GradObj','off');
[bk,fval,exitflag,output,grad] = fminunc(@f_bk,RateParams,options);
Enew.k  = bk;


    function [f g]= f_bk(RateParams)

        kx      = RateParams'*Sim.x;
        f_kdt   = repmat(exp(kx)*Sim.dt,Sim.N,1);
        f       = -sum(S.w_b(sp).*log(1-exp(-f_kdt(sp))))...
            +sum(S.w_b(nosp).*f_kdt(nosp));

        ef      = exp(f_kdt);
        g       = x(:,sp)*(S.w_b(sp).*f_kdt(sp)./( ef(sp)-1))...
            -x(:,nosp)*(S.w_b(nosp).*f_kdt(nosp));
    end %function f_bk


%% MLE for calcium parameters
% fprintf('estimating calcium parammeters\n')
% ve_x        = quadprog(M.Q, M.L,[],[],[],[],[0 -inf],[],[1/E.tau_c E.beta]);
% Enew.tau_c  = 1/ve_x(1);
% Enew.beta   = ve_x(2);
% 
% Enew.sigma_c= sqrt((ve_x'*M.Q*ve_x + M.L'*ve_x)/Sim.T);
% 
% % %% MLE for spike history parameters
% % for m=1:Sim.M
% %     Enew.sigma_h(m)= sum(M.v{m})/Sim.T;
% % end
% 
% %% MLE for observation parameters
% O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
% Oind    = find(~isnan(O));
% Enew.sigma_o = sqrt(sum(sum(S.w_b(:,Oind).*(repmat(O(Oind),Sim.N,1)-S.C(:,Oind)).^2)/Sim.T));


%% this code maximizes the fully observed case
% %%%%%%%% maximize alpha and beta
% Ps=S.w_b(:,2:end);
% Ps=sqrt(Ps(:));
% ns=S.n(:,2:end);
% ns=ns(:);
% J=[Ps'.*repmat(R.C(2:end)*Sim.dt,1,Sim.N); -Ps'.*ns']*[Ps'.*repmat(R.C(2:end)*Sim.dt,1,Sim.N); -Ps'.*ns']';
% 
% Ps=S.w_b(:,2:end)';
% Ps=Ps(:);
% ns=S.n(:,2:end)';
% ns=ns(:);
% f=[repmat(R.C(2:end)*Sim.dt,1,Sim.N); -ns']*(repmat(R.C(2:end)-R.C(1:end-1),1,Sim.N)'.*Ps);
% 
% [minab,fval,exitflag,output,lambda] = quadprog(J, f,[],[],[],[],[0 -inf],[],[1 1]);
% Enew.tau_c=1/minab(1);
% Enew.beta=minab(2);
% 
% %%%%%%%% maximize sigma
% DelCs=repmat(R.C(2:end)-R.C(1:end-1)*exp(-minab(1)*Sim.dt),1,Sim.N);
% DelCs=DelCs(:);
% Enew.sigma_c=sqrt(sum(Ps.*(DelCs-minab(2)*ns).^2/Sim.dt)/sum(Ps(:)));

end %function smc_em_bern_mstep
