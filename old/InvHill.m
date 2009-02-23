function C = InvHill(F,n,k_d)
% function C = InvHill(F,n,k_d)
% F(F<0)=0;
% C = ((F*k_d)./(1-F)).^(1/n);

F(F<0)=0;
C = ((F*k_d)./(1-F)).^(1/n);