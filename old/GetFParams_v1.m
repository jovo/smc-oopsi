function xhat = GetFParams_v1

[Sim P] = InitializeStuff;

Sim.freq    = 5;
Sim.Nsec    = 1;               %# of sec
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem         = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                    %fix number of steps
end
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec-Sim.dt*rem;%time vector

Sim.M       = 1;
Sim.StimDim = 5;                %stimulus dimensionality
Sim.x       = rand(Sim.StimDim,Sim.T)*3-1;
% for t=1:Sim.freq:Sim.T
%     Sim.x(:,t:t+Sim.freq-1)=repmat(Sim.x(:,t),1,Sim.freq);
% end
Sim.x(1,:)  = ones(1,Sim.T);
P.k         = 1.3*ones(Sim.StimDim,1);%linear kernel
P.k(2:2:end)= P.k(2:2:end)-.5;
P.omega     = -6*P.k(1);
P.tau_h     = .025;                %decay rate for spike history terms
% P.gamma     = 1e-5;             %var gainh
% P.zeta      = 1e-5;            %var offset
P.C_0       = 1;
% P.A         = .1;
% P.tau_c     = .5;

R = smc_em_bern_real_exp_v5(Sim,P);

S       = Hill_v1(P,R.C);
options = optimset('GradObj','off','Hessian','off');
x0      = [P.alpha P.beta P.gamma P.zeta]';
zeroy   = zeros(Sim.T,1);
oney    = ones(Sim.T,1);
[xhat]  = fmincon(@fnlogL2,x0,[],[],[],[],[0 0 0 0],[inf inf inf inf]); %,[],options);

    function logL = fnlogL2(x) %[logL dlogL ddlogL] = fnlogL2(x)        %this function compute log L = log P(O|H)
        varf=[zeroy zeroy S' oney]*x;
        logL = sum((((R.F'-[S' oney zeroy zeroy]*x).^2)./varf+log(varf))/2);
%         if nargout > 1
%             dlogL=-(o-alpha*C^n/(C^n+k_d)-beta)/(gamma*C^n/(C^n+k_d)+zeta)*(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)+1/2*(o-alpha*C^n/(C^n+k_d)-beta)^2/(gamma*C^n/(C^n+k_d)+zeta)^2*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)-1/2*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)/(gamma*C^n/(C^n+k_d)+zeta);
%             if nargout > 2
%                 ddlogL = -(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)^2/(gamma*C^n/(C^n+k_d)+zeta)+2*(o-alpha*C^n/(C^n+k_d)-beta)/(gamma*C^n/(C^n+k_d)+zeta)^2*(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)-(o-alpha*C^n/(C^n+k_d)-beta)/(gamma*C^n/(C^n+k_d)+zeta)*(-alpha*C^n*n^2/C^2/(C^n+k_d)+alpha*C^n*n/C^2/(C^n+k_d)+3*alpha*(C^n)^2*n^2/C^2/(C^n+k_d)^2-2*alpha*(C^n)^3/(C^n+k_d)^3*n^2/C^2-alpha*(C^n)^2/(C^n+k_d)^2*n/C^2)-(o-alpha*C^n/(C^n+k_d)-beta)^2/(gamma*C^n/(C^n+k_d)+zeta)^3*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)^2+1/2*(o-alpha*C^n/(C^n+k_d)-beta)^2/(gamma*C^n/(C^n+k_d)+zeta)^2*(gamma*C^n*n^2/C^2/(C^n+k_d)-gamma*C^n*n/C^2/(C^n+k_d)-3*gamma*(C^n)^2*n^2/C^2/(C^n+k_d)^2+2*gamma*(C^n)^3/(C^n+k_d)^3*n^2/C^2+gamma*(C^n)^2/(C^n+k_d)^2*n/C^2)-1/2*(gamma*C^n*n^2/C^2/(C^n+k_d)-gamma*C^n*n/C^2/(C^n+k_d)-3*gamma*(C^n)^2*n^2/C^2/(C^n+k_d)^2+2*gamma*(C^n)^3/(C^n+k_d)^3*n^2/C^2+gamma*(C^n)^2/(C^n+k_d)^2*n/C^2)/(gamma*C^n/(C^n+k_d)+zeta)+1/2*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)^2/(gamma*C^n/(C^n+k_d)+zeta)^2;
%             end
%         end
    end

syms x alpha beta gamma zeta F S

varf=[0 0 S 1];
x=[alpha beta gamma zeta]';

lnlik=-log(varf)/2-(F-[0 0 S 1]*x).^2/varf;

end