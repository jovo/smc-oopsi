function F = GenHill_v3(P,C)
% generalized hill model
C(C<0)  = 0;
F       = P.alpha*(C).^P.n./((C).^P.n+P.k_d) + P.beta;