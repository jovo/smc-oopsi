function F = GenHill_v1(P,C)
% generalized hill model

F = P.alpha*C.^P.n./(C.^P.n+P.k_d) + P.beta;