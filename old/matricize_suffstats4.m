clear; clc, 
T       = 5;
Sim.N   = 4;
Sim.M   = 1;
Sim.dt  = 0.005;
S.C     = rand(Sim.N,T);
S.n     = rand(Sim.N,T);
S.h     = rand(Sim.N,T,Sim.M);
PHH     = rand(Sim.N); PHH = PHH/sum(PHH(:));

Nparams = 3;
M.Q     = zeros(Nparams,Nparams);
M.L     = zeros(Nparams,1);
M.u     = 0;
for m = 1:Sim.M,
    M.v{m} = 0;
end

P.tau_h = 1;
P.tau_c = 1;
P.A     = 2;
P.a     = Sim.dt/P.tau_c;
P.C_0   = 0;

t=T;
for i=1:Sim.N
    for j=1:Sim.N
        A{i,j}  = [S.C(j,t-1)*Sim.dt, -S.n(i,t), -Sim.dt];
        b{i,j}  = S.C(i,t)-S.C(j,t-1);
        Q{i,j}  = A{i,j}'*A{i,j};
        L{i,j}  = 2*A{i,j}'*b{i,j};
        u{i,j}  = (S.C(i,t)-(1-Sim.dt/P.tau_c)*S.C(j,t-1)-P.A*S.n(i,t)-P.a*P.C_0).^2;
        M.Q     = M.Q+PHH(i,j)*Q{i,j};
        M.L     = M.L+PHH(i,j)*L{i,j};
        M.u     = M.u+PHH(i,j)*u{i,j};
        for m=1:Sim.M
            M.v{m}=M.v{m}+PHH(i,j)*...
                (S.h(i,t,m)-(1-Sim.dt/P.tau_h)*S.h(j,t-1,m)-S.n(i,t)).^2/Sim.dt;
        end
    end
end

% A2 = [S.C(:,t-1)*Sim.dt, -S.n(:,t), -repmat(Sim.dt,Sim.N,1)];
% b2 = S.C(:,t)-S.C(:,t-1);
% 
% H = A2'*A2;
% f = A2'*b2;
% 
% norm(A2'*PHH'*PHH*A2-M.Q)
% norm(2*A2'*PHH*b2-M.L)

Z.oney  = ones(Sim.N,1);
Z.n1    = S.n(:,t);
Z.C0    = S.C(:,t-1);
Z.C1    = S.C(:,t);
Z.C1mat = Z.C1(:,Z.oney);
Z.C0mat = Z.C0(:,Z.oney);
Z.PHH   =  PHH/sum(PHH(:));
        
C0dt    = Z.C0*Sim.dt;
bmat    = Z.C1mat-Z.C0mat';
bPHH    = Z.PHH.*bmat;

%%
% u
syms A a C_0
umat=(repmat(S.C(:,t)-A*S.n(:,t),1,Sim.N)-repmat((1-a)*S.C(:,t-1),1,Sim.N)'-a*C_0).^2;
A=P.A; a=P.a; C_0=P.C_0;
umat=eval(umat);
diffu =  sum(sum(Z.PHH.*umat))-M.u;
fprintf('diff u = %d\n',diffu)

% v
vmat1=(repmat(S.h(:,t,1)-S.n(:,t),1,Sim.N)...
    -repmat((1-Sim.dt/P.tau_h)*S.h(:,t-1,1),1,Sim.N)').^2/Sim.dt;
v2=sum(PHH(:).*vmat1(:));
fprintf('diff v2 = %d\n', M.v{1}-v2)

bmat = repmat(S.C(:,t),1,Sim.N)-repmat(S.C(:,t-1)',Sim.N,1);
bmatT=bmat'; 
bPHH = PHH.*bmat;
PHHT = PHH';
C01  = (S.C(:,t-1)*Sim.dt);
n01  = -S.n(:,t);

% L1
L1  = sum(bPHH*C0dt); %= sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)));
fprintf('diff L1 = %d\n',M.L(1)-2*L1)

% L2
L2  = - sum(bPHH'*Z.n1); %sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:)));
fprintf('diff L2 = %d\n', M.L(2)-2*L2)

% L3
L3 = - Sim.dt*sum(bPHH(:)); %-2*Sim.dt*sum(bPHH(:));
fprintf('diff L3 = %d\n', M.L(3)-2*L3)

% Q11
Q11 = sum(Z.PHH*(C0dt.^2)); %sum(sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1)));
fprintf('diff Q11 = %d\n',M.Q(1,1)-Q11)

% Q12
Q12 = - Z.n1'*Z.PHH*C0dt; % Xterms  = -S.n(:,t)*S.C(:,t-1)'*Sim.dt; Q12b    = sum(sum(PHH.*Xterms));
fprintf('diff Q12 = %d\n',M.Q(1,2)-Q12)

% Q13
Q13 = sum(sum(-Z.PHH.*Z.C0mat'*Sim.dt^2)); % Xterms2 = repmat(-S.C(:,t-1)',Sim.N,1)*Sim.dt^2; Q13     = sum(sum(PHH.*Xterms2));
fprintf('diff Q13 = %d\n',M.Q(1,3)-Q13)

% Q22
Q22 = sum(Z.PHH'*(Z.n1.^2)); % Q22b = sum(sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1)));
fprintf('diff Q22 = %d\n',M.Q(2,2)-Q22)

% Q23
Q23 = sum(sum(Z.PHH(:).*repmat(Z.n1,Sim.N,1))*Sim.dt); % Q23 = sum(sum(PHH(:).*repmat((S.n(:,t)),Sim.N,1)))*Sim.dt;
fprintf('diff Q23 = %d\n',M.Q(2,3)-Q23)

% Q33
Q33 = sum(Z.PHH(:))*Sim.dt^2; % Q33 = sum(PHH(:))*Sim.dt^2;
fprintf('diff Q33 = %d\n\n',M.Q(3,3)-Q33)

% %%
% clear, %clc
% T       = 5;
% Sim.N   = 100;
% Sim.M   = 1;
% Sim.dt  = 0.005;
% S.C     = rand(Sim.N,T);
% S.n     = rand(Sim.N,T);
% S.h     = rand(Sim.N,T,Sim.M);
% PHH     = rand(Sim.N); PHH = PHH/sum(PHH(:));
% 
% P.tau_h = 1;
% P.tau_c = 1;
% P.A     = 2;
% t       = 5;
% 
% Nparams = 3;
% M.Q     = zeros(Nparams,Nparams);
% M.L     = zeros(Nparams,1);
% M.u     = 0;
% for m = 1:Sim.M,
%     M.v{m} = 0;
% end
% 
% Z.oney  = ones(Sim.N,1);
% Z.n1    = S.n(:,t);                        
% Z.C0    = S.C(:,t-1);                      
% Z.C1    = S.C(:,t);
% Z.C1mat = Z.C1(:,Z.oney);                  
% Z.C0mat = Z.C0(:,Z.oney);
% Z.PHH   =  PHH/sum(PHH(:));
% 
% C0dt    = Z.C0*Sim.dt;
% bmat    = Z.C1mat-Z.C0mat';
% bPHH    = Z.PHH.*bmat;
% 
% 
% PHHT    = PHH';
% Xterms  = -S.n(:,t)*S.C(:,t-1)'*Sim.dt;
% Xterms2 = repmat(-S.C(:,t-1)',Sim.N,1)*Sim.dt^2;
% C01     = (S.C(:,t-1)*Sim.dt);
% n01     = -S.n(:,t);
% C11     = S.C(:,t);
% bPHH    = PHH.*bmat;
% bmatT   = bmat';
% 
% frac_err_Q11 = (sum(Z.PHH*(C0dt.^2)) - sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1)))/sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1));
% fprintf('frac_err_Q11 = %d\n',frac_err_Q11)
% 
% frac_err_Q12 = (-Z.n1'*Z.PHH*C0dt - sum(PHH(:).*Xterms(:)))/sum(PHH(:).*Xterms(:));
% fprintf('frac_err_Q12 = %d\n',frac_err_Q12)
% 
% frac_err_Q13 = (sum(sum(-Z.PHH.*Z.C0mat'*Sim.dt^2)) - sum(PHH(:).*Xterms2(:)))/sum(PHH(:).*Xterms2(:));
% fprintf('frac_err_Q13 = %d\n',frac_err_Q13)
% 
% frac_err_Q22 = (sum(Z.PHH'*(Z.n1.^2)) - sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1)))/sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1));
% fprintf('frac_err_Q22 = %d\n',frac_err_Q22)
% 
% frac_err_Q23 =  (sum(sum(Z.PHH(:).*repmat(Z.n1,Sim.N,1))*Sim.dt) - sum(sum(PHH(:).*repmat(S.n(:,t),Sim.N,1)))*Sim.dt)/(sum(sum(PHH(:).*repmat((S.n(:,t)),Sim.N,1)))*Sim.dt);
% fprintf('frac_err_Q23 = %d\n',frac_err_Q23)
% 
% frac_err_Q33 = (sum(Z.PHH(:))*Sim.dt^2 - sum(PHH(:))*Sim.dt^2)./(sum(PHH(:))*Sim.dt^2);
% fprintf('frac_err_Q33 = %d\n',frac_err_Q33)
% 
% frac_err_L1  = (2*sum(bPHH*C0dt) - sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:))))/sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)));
% fprintf('frac_err_L1 = %d\n',frac_err_L1)
% 
% frac_err_L2  = (-2*sum(bPHH'*Z.n1) - sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:))))/sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:)));
% fprintf('frac_err_L2 = %d\n',frac_err_L2)
% 
% fprintf('frac_err_L3 too easy to check\n')