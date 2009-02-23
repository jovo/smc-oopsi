function M = smc_em_bern_moments_v1(Sim,S,R)

%%%%%%%%  compute moments
M.Ibar = sum(S.w_b.*S.I,1);
M.Ivar = sum((repmat(M.Ibar,Sim.N,1)-S.I).^2)/Sim.N;

M.Cbar = sum(S.w_b.*S.C,1);
M.Cvar = sum((repmat(M.Cbar,Sim.N,1)-S.C).^2)/Sim.N;

% M.pbar = sum(S.w_b.*S.p,1);
% M.pvar = sum((repmat(M.pbar,Sim.N,1)-S.p).^2)/Sim.N;
% 
% for l=1:Sim.M
%     M.hbar(l,:) = sum(S.w_b.*S.h(:,:,l),1);
%     M.hvar(l,:) = sum((repmat(M.hbar(l,:),Sim.N,1)-S.h(:,:,l)).^2)/Sim.N;
% end

%%%%%%%%  compute MSE's
M.mse.I=sum(S.w_b.*(S.I-repmat(R.I,Sim.N,1)).^2);
M.mse.C=sum(S.w_b.*(S.C-repmat(R.C,Sim.N,1)).^2);
end %function