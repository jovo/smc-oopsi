%% scalar gaussian
clear, clc

%define 2 gaussians
y=round(rand*10)+1;

mu1=round(rand*10)+1;
sigma1=round(rand*10)+1;

mu2=round(rand*10)+1;
sigma2=round(rand*10)+1;

%define mean and variance of product
v=(1/sigma1^2+1/sigma2^2)^(-1);
m=v*(mu1/(sigma1^2)+mu2/(sigma2^2));

%define normalization constant
Z=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(m^2/v -mu1^2/sigma1^2-mu2^2/sigma2^2));

%check if corrent
1/sqrt(2*pi*sigma1^2)*exp(-.5*(y-mu1)^2/sigma1^2)...
    *1/sqrt(2*pi*sigma2^2)*exp(-.5*(y-mu2)^2/sigma2^2)

Z/sqrt(2*pi*v)*exp(-.5*(y-m)^2/v)


%rewrite Z as a gaussian
v1=sigma1^2+sigma2^2;
m1=mu2;
k=1;

%check if correct
Z
k/sqrt(2*pi*v1)*exp(-.5*(mu1-m1)^2/v1)


%rewrite for affine functions of y
A=round(rand*10)+1;
b=round(rand*10)+1;

x=A*y+b;

1/sqrt(2*pi*v)*exp(-.5*(x-m)^2/v)

vy=v/A^2;
my=(m-b)/A;
k=1/A;

k/sqrt(2*pi*vy)*exp(-.5*(y-my)^2/vy)


%% matrices
%define 2 gaussians
y=round(rand(100,2)*10)+1;

mu1=round(rand(100,2)*10)+1;
sigma1=round(rand(100,2)*10)+1;

mu2=round(rand(100,2)*10)+1;
sigma2=round(rand(100,2)*10)+1;

%define mean and variance of product
v=(1./sigma1.^2+1./sigma2.^2).^(-1);
m=v.*(mu1./(sigma1.^2)+mu2./(sigma2.^2));

%define normalization constant
Z=sqrt(2*pi*v)./(2*pi*sigma1.*sigma2) .* exp(0.5*(m.^2./v -mu1.^2./sigma1.^2-mu2.^2./sigma2.^2));

%check if corrent
prod1=1./sqrt(2*pi*sigma1.^2).*exp(-.5*(y-mu1).^2./(sigma1.^2))...
    .*1./sqrt(2*pi*sigma2.^2).*exp(-.5*(y-mu2).^2./(sigma2.^2));

prod2=Z./sqrt(2*pi.*v).*exp(-.5*(y-m).^2./v);

norm(prod1 - prod2)

%rewrite Z as a gaussian
v1=sigma1.^2+sigma2.^2;
m1=mu2;
k=1;

%check if correct
norm(Z-k./sqrt(2*pi.*v1).*exp(-.5*(mu1-m1).^2./v1))


%rewrite for affine functions of y
A=round(rand(100,2)*10)+1;
b=round(rand(100,2)*10)+1;

x=A.*y+b;

g1=1./sqrt(2*pi.*v).*exp(-.5*(x-m).^2./v);

vy=v./A.^2;
my=(m-b)./A;
k=1./A;

g2=k./sqrt(2*pi*vy).*exp(-.5*(y-my).^2./vy)

norm(g1-g2)

%% other crap
% Z55=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(v^2*(mu1/sigma1^2+mu2/sigma2^2)^2/v -mu1^2/sigma1^2-mu2^2/sigma2^2))
% 
% Z56=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(v*(mu1/sigma1^2+mu2/sigma2^2)^2 -mu1^2/sigma1^2-mu2^2/sigma2^2))
% 
% Z57=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(v*((mu1/sigma1^2)^2+2*mu1*mu2/(sigma1^2*sigma2^2)+(mu2/sigma2^2)^2)...
% -mu1^2/sigma1^2-mu2^2/sigma2^2))
% 
% Z58=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(...
%     mu1^2*(v/sigma1^4-1/sigma1^2)+2*mu1*mu2*v/(sigma1^2*sigma2^2)+mu2^2*(v/sigma2^4-1/sigma2^2)))
% 
% Z59=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(...
%     mu1^2*(sigma1^2*sigma2^2/(sigma1^4*(sigma1^2+sigma2^2))-1/sigma1^2)...
%     +2*mu1*mu2*v/(sigma1^2*sigma2^2)...
%     +mu2^2*(sigma1^2*sigma2^2/(sigma2^4*(sigma1^2+sigma2^2))-1/sigma2^2)))
% 
% Z60=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(...
%     mu1^2*(sigma2^2/(sigma1^2*(sigma1^2+sigma2^2))-1/sigma1^2)...
%     +2*mu1*mu2*v/(sigma1^2*sigma2^2)...
%     +mu2^2*(sigma1^2/(sigma2^2*(sigma1^2+sigma2^2))-1/sigma2^2)))
% 
% Z61=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(...
%     mu1^2*((sigma2^2-sigma1^2-sigma2^2)/(sigma1^2*(sigma1^2+sigma2^2)))...
%     +2*mu1*mu2*v/(sigma1^2*sigma2^2)...
%     +mu2^2*((sigma1^2-sigma1^2-sigma2^2)/(sigma2^2*(sigma1^2+sigma2^2)))))
% 
% Z62=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(...
%     mu1^2*(-1/(sigma1^2+sigma2^2))...
%     +2*mu1*mu2*v/(sigma1^2*sigma2^2)...
%     +mu2^2*(-1/(sigma1^2+sigma2^2))))
% 
% Z64=sqrt(2*pi*v)/(2*pi*sigma1*sigma2) * exp(0.5*(...
%     mu1^2*(-1/(sigma1^2+sigma2^2))+2*mu1*mu2/(sigma1^2+sigma2^2)-mu2^2/(sigma1^2+sigma2^2)))
% 
% -2*log(Z64)
% logZ_2=mu1^2/(sigma1^2+sigma2^2)-2*mu1*mu2/(sigma1^2+sigma2^2)+mu2^2/(sigma1^2+sigma2^2)-2*log(sqrt(2*pi*v)/(2*pi*sigma1*sigma2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v1=solve('1/(sigma1^2+sigma2^2)=1/v1','v1')
% m1=solve('mu2/(sigma1^2+sigma2^2)=m1/v1','m1')
% k=solve('mu2^2/(sigma1^2+sigma2^2)-2*log(sqrt(2*pi*v)/(2*pi*sigma1*sigma1))=m1^2/v1-2*log(k/sqrt(2*pi*v1))','k')