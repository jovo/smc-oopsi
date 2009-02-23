clear all, clc
T       = 5;
Sim.N   = 100;
Sim.M   = 1;
Sim.dt  = 0.005;
S.C     = rand(Sim.N,T);
S.n     = rand(Sim.N,T);
S.h     = rand(Sim.N,T,Sim.M);
PHH     = rand(Sim.N); PHH = PHH/sum(PHH(:));
P.tau_h = 1;
P.tau_c = 1;
P.beta  = 2;
t       = 5;

M.Q     = zeros(2,2);
M.L     = zeros(2,1);
M.u     = 0;
for m = 1:Sim.M,
    M.v{m} = 0;
end

PHHT        = PHH';
Xterms      = -S.n(:,t)*S.C(:,t-1)'*Sim.dt;
C01         = (S.C(:,t-1)*Sim.dt);
n01         = -S.n(:,t);
C11         = S.C(:,t);

M.Q(1,1)    = M.Q(1,1) + sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1));
diffQ11 = (C01'*PHH*C01 - sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1)));
fprintf('diffQ11 = %d\n',diffQ11)

M.Q(1,2)    = M.Q(1,2) + sum(PHH(:).*Xterms(:));
diffQ12 = n01'*PHH*C01 - sum(PHH(:).*Xterms(:));
fprintf('diffQ12 %d\n',diffQ12)

M.Q(2,2)    = M.Q(2,2) + sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1));
diffQ22 = n01'*PHH*n01 - sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1));
fprintf('diffQ22 %d\n',diffQ22)

diffQ22 = n01'*PHH*n01 - sum(PHH(:).*repmat((-S.n(:,t))'.^2,1,Sim.N)');
fprintf('diffQ22 %d\n',diffQ22)

bmat        = repmat(S.C(:,t),1,Sim.N)-repmat(S.C(:,t-1)',Sim.N,1);
bmatT       = bmat';
M.L(1)      = M.L(1) + sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)));
bPHH = PHH.*bmat;
diffL1 = 2*sum(bPHH*C01) - sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)));
fprintf('diffL1 %d\n',diffL1)

M.L(2)      = M.L(2) + sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:)));
diffL2 = 2*sum((PHH.*bmat)*n01) - sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:)));
fprintf('diffL2 %d\n',diffL2)
