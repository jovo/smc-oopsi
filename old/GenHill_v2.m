function F = GenHill_v2(P,C)
% generalized hill model
C(C<0)  = 0;
F       = P.alpha*(C-P.C_0).^P.n./((C-P.C_0).^P.n+P.k_d) + P.beta;