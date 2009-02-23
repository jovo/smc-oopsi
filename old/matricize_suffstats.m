clear all, clc, 
T       = 5;
Sim.N   = 3;
Sim.M   = 1;
Sim.dt  = 0.005;
S.C     = rand(Sim.N,T);
S.n     = rand(Sim.N,T);
S.h     = rand(Sim.N,T,Sim.M);
PHH     = rand(Sim.N); PHH = PHH/sum(PHH(:));
P.tau_h = 1;

M.Q     = zeros(2,2);
M.L     = zeros(2,1);
M.u     = 0;
for m = 1:Sim.M,
    M.v{m} = 0;
end
syms    Ptau_c Pbeta
P.tau_c = 1;
P.beta  = 2;

t=T;
for i=1:Sim.N
    for j=1:Sim.N
        A{i,j}  = [S.C(j,t-1)*Sim.dt, -S.n(i,t)];
        b{i,j}  = S.C(i,t)-S.C(j,t-1);
        Q{i,j}  = A{i,j}'*A{i,j};
        L{i,j}  = 2*A{i,j}'*b{i,j};
        u{i,j}  = (S.C(i,t)-(1-Sim.dt/P.tau_c)*S.C(j,t-1)-P.beta*S.n(i,t)).^2/Sim.dt;
        M.Q     = M.Q+PHH(i,j)*Q{i,j};
        M.L     = M.L+PHH(i,j)*L{i,j};
        M.u     = M.u+PHH(i,j)*u{i,j};
        for m=1:Sim.M
            M.v{m}=M.v{m}+PHH(i,j)*...
                (S.h(i,t,m)-(1-Sim.dt/P.tau_h)*S.h(j,t-1,m)-S.n(i,t)).^2/Sim.dt;
        end
    end
end
Qmat=[Q{1,1} Q{1,2} Q{1,3}; Q{2,1} Q{2,2} Q{2,3}; Q{3,1} Q{3,2} Q{3,3}];


%% u
umat=(repmat(S.C(:,t)-P.beta*S.n(:,t),1,Sim.N)...
    -repmat((1-Sim.dt/P.tau_c)*S.C(:,t-1),1,Sim.N)').^2/Sim.dt;
diffu =  sum(PHH(:).*umat(:))-M.u;
fprintf('u2 diff = %d\n',diffu)

%% v
vmat1=(repmat(S.h(:,t,1)-S.n(:,t),1,Sim.N)...
    -repmat((1-Sim.dt/P.tau_h)*S.h(:,t-1,1),1,Sim.N)').^2/Sim.dt;
v2=sum(PHH(:).*vmat1(:));
fprintf('diff v2 = %d\n', M.v{1}-v2)

%% L2
bmat = repmat(S.C(:,t),1,Sim.N)-repmat(S.C(:,t-1)',Sim.N,1);
L2b  = sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:)));
fprintf('diff L2 = %d\n', M.L(2)-L2b)

bPHH = PHH.*bmat;
C01         = (S.C(:,t-1)*Sim.dt);
n01         = -S.n(:,t);
fprintf('\ntrue L2 = %g,\n mat L2 = %g\n\n',M.L(2), 2*sum((PHH'.*bmat')*n01))
 
%% L1
bmatT=bmat'; PHHT=PHH';
L1c  = sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)));
fprintf('diff L1 = %d\n',M.L(1)-L1c)


%% Q11
ind11 = [1,3,5,13,15,17,25,27,29];
Q11 = sum(sum(PHH(:).*Qmat(ind11)'));
PHHT = PHH';
Q11b = sum(sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1)));
fprintf('diff Q11 = %d\n',M.Q(1,1)-Q11b)

%% Q12
ind12 = [2,4,6,14,16,18,26,28,30];
Q12 = sum(sum(PHH(:).*Qmat(ind12)'));

Xterms=-S.n(:,t)*S.C(:,t-1)'*Sim.dt;
Q12b = sum(sum(PHH(:).*Xterms(:)));
fprintf('diff Q12 = %d\n',M.Q(1,2)-Q12b)

%% Q22
ind22 = [8,10,12,20,22,24,32,34,36];
Q22 = sum(sum(PHH(:).*Qmat(ind22)'));
Q22b = sum(sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1)));
fprintf('diff Q22 = %d\n',M.Q(2,2)-Q22b)
fprintf('\ntrue Q22 = %g,\n mat Q22 = %g\n\n',M.Q(2,2), n01'*PHH'*n01)
fprintf('\ntrue Q22 = %g,\n mat Q22 = %g\n\n',M.Q(2,2), sum(PHH'*(n01.^2)))


% A = [C01 n01];
% A'*PHH*A
% 
% A = [C01'; n01'];
% A*PHH*A'

% Amat1=[A{1,1} A{1,2} A{1,3}; A{2,1} A{2,2} A{2,3}; A{3,1} A{3,2} A{3,3}];
%
% for i=1:Sim.N
%     for j=1:Sim.N
%         A   = [S.C(j,t-1)*Sim.dt, -S.n(i,t)];
%         b   = S.C(i,t)-S.C(j,t-1);
%         Q   = A'*A;
%         L   = 2*A'*b;
%         M.Q = M.Q+PHH(i,j)*Q;
%         M.L = M.L+PHH(i,j)*L;
%     end
% end

% inda = [1:3,7:9,13:15];
% L1a  = sum(PHH(:).*(2*(Amat1(inda)'.*bmat(:))))
% indb = [1,4,7,2,5,8,3,6,9];
% L1b  = sum(PHH(indb)'.*(2*repmat(S.C(:,t-1)*Sim.dt,3,1).*bmat(indb)'))
% ind  = [4:6,10:12,16:18];
% L2   = sum(PHH(:).*(2*(Amat1(ind)'.*bmat(:))));

% Ptau_c = 1;
% Pbeta  = 2;
% M.u
