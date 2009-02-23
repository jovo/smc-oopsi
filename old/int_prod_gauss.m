% function [m v] = int_prod_gauss(mu1, sigma1, mu2, sigma2)
%
% v=sigma1.^2+sigma2.^2;
% m=mu2;

mu_ck   = [B.a*S.C(:,k-1)+B.beta B.a*S.C(:,k-1)];   %mean when spiked and not spiked

sig2_zk = repmat((1/B.sig2_c + 1./S.sig2_o(:,k)).^(-1),1,2);        %update sig2_zk
mu_zk   = sig2_zk.*(mu_ck./B.sig2_c + repmat(S.mu_o(:,k)./S.sig2_o(:,k),1,2));                     %update mu_zk

ln_Z    = 0.5*log(sig2_zk) - 0.5*(2*pi*repmat(S.sig2_o(:,k),1,2)*B.sig2_c) + 0.5*(mu_zk.^2./sig2_zk-repmat(S.mu_o(:,k).^2./S.sig2_o(:,k),1,2)-mu_ck.^2/B.sig2_c);  %update ln_z (ignore constant)

mu1=mu_ck;
sigma1_2=B.sig2_c;

mu2=repmat(S.mu_o(:,k),1,2);
sigma2_2=repmat(S.sig2_o(:,k),1,2);

%define mean and variance of product
v=(1./sigma1_2+1./sigma2_2).^(-1);
m=v.*(mu1./(sigma1_2)+mu2./(sigma2_2));

Z=sqrt(2*pi*v)./(2*pi*sqrt(sigma1_2.*sigma2_2)) .* exp(0.5*(m.^2./v -mu1.^2./sigma1_2-mu2.^2./sigma2_2));

v1=sigma1_2+sigma2_2;
m1=mu2;

norm(Z-1./sqrt(2*pi.*v1).*exp(-.5*(mu1-m1).^2./v1))

logZ=log(sqrt(2*pi*v)./(2*pi*sqrt(sigma1_2.*sigma2_2)) .* exp(0.5*(m.^2./v -mu1.^2./sigma1_2-mu2.^2./sigma2_2)));
logZ2=log(sqrt(2*pi*v)./(sqrt(2^2*pi^2*sigma1_2.*sigma2_2)) .* exp(0.5*(m.^2./v -mu1.^2./sigma1_2-mu2.^2./sigma2_2)));
logZ3=log(sqrt(2*pi*v./(2^2*pi^2*sigma1_2.*sigma2_2)) .* exp(0.5*(m.^2./v -mu1.^2./sigma1_2-mu2.^2./sigma2_2)));
logZ4=log(sqrt(v./(2*pi*sigma1_2.*sigma2_2)) .* exp(0.5*(m.^2./v -mu1.^2./sigma1_2-mu2.^2./sigma2_2)));
logZ5=0.5*log(v./(2*pi*sigma1_2.*sigma2_2)) + 0.5*(m.^2./v -mu1.^2./sigma1_2-mu2.^2./sigma2_2);

logZ6=log(1./sqrt(2*pi.*v1).*exp(-.5*(mu1-m1).^2./v1));
logZ7=-0.5*log(2*pi.*v1)-.5*(mu1-m1).^2./v1;
logZ8=-0.5*log(2*pi.*v1)-.5*(mu1-mu2).^2./v1;



%define normalization constant
ln_Z3 = 0.5*log(v) - 0.5*log(2*pi*sigma1.*sigma2) + 0.5*(m.^2./v -mu1.^2./sigma1.^2-mu2.^2./sigma2.^2);

ln_Z2=log(Z);
norm(ln_Z3-ln_Z2)




sig2_Ik  = (repmat(S.sig2_I(:,k),1,2)+B.sig2_c);
ln_G    = -0.5*log(2*pi*sig2_Ik) -0.5*(repmat(S.mu_o(:,k),1,2) - mu_ck).^2./sig2_Ik;




% %check if corrent
% prod1=1./sqrt(2*pi*sigma1.^2).*exp(-.5*(y-mu1).^2./(sigma1.^2))...
%     .*1./sqrt(2*pi*sigma2.^2).*exp(-.5*(y-mu2).^2./(sigma2.^2));
%
% prod2=Z./sqrt(2*pi.*v).*exp(-.5*(y-m).^2./v);
%
% norm(prod1 - prod2)
%
