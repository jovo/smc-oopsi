function F = Hill_v1(P,C)
% generalized hill model
% code:
% C(C<0)  = 0;
% F       = C.^P.n./(C.^P.n+P.k_d);

C(C<0)  = 0;
F       = C.^P.n./(C.^P.n+P.k_d);