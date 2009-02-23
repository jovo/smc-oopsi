PHHT        = PHH';
Xterms      = -S.n(:,t)*S.C(:,t-1)'*Sim.dt;
C01         = (S.C(:,t-1)*Sim.dt);
n01         = -S.n(:,t);
C11         = S.C(:,t);

M.Q(1,1)    = M.Q(1,1) + sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1));
C01'*PHH*C01 == sum(PHHT(:).*repmat((S.C(:,t-1)*Sim.dt).^2,Sim.N,1))

M.Q(1,2)    = M.Q(1,2) + sum(PHH(:).*Xterms(:));
-S.n(:,t)'*PHH*S.C(:,t-1)*Sim.dt == sum(PHH(:).*Xterms(:))
n01'*PHH*C01 ==sum(PHH(:).*Xterms(:))

M.Q(2,2)    = M.Q(2,2) + sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1));
n01'*PHH*n01  == sum(PHH(:).*repmat((-S.n(:,t)).^2,Sim.N,1))

bmat        = repmat(S.C(:,t),1,Sim.N)-repmat(S.C(:,t-1)',Sim.N,1);
bmatT       = bmat';
M.L(1)      = M.L(1) + sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)));
bPHH = PHH.*bmat
2*sum(bPHH*C01) == sum(PHHT(:).*(2*repmat(S.C(:,t-1)*Sim.dt,Sim.N,1).*bmatT(:)))

M.L(2)      = M.L(2) + sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:)));
2*sum((PHH.*bmat)*n01) == sum(PHH(:).*(2*repmat(-S.n(:,t),Sim.N,1).*bmat(:)))

umat        =(repmat(S.C(:,t)-P.beta*S.n(:,t),1,Sim.N)...
    -repmat((1-Sim.dt/P.tau_c)*S.C(:,t-1),1,Sim.N)').^2/Sim.dt;
M.u         = M.u + sum(PHH(:).*umat(:));

for m=1:Sim.M
    vmat    = (repmat(S.h(:,t,m)-S.n(:,t),1,Sim.N)...
        -repmat((1-Sim.dt/P.tau_h)*S.h(:,t-1,m),1,Sim.N)').^2/Sim.dt;
    M.v{m}  = M.v{m} + sum(PHH(:).*vmat(:));
end
[M.Q, M.Q, M.L, M.L], [M.v{1}; M.v{1}], %[M.u; M.v{1}]
