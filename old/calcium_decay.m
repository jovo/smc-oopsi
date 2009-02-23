clear; clc;
dt=.001;
t=0:dt:10;
T=length(t);
tau_c = 1;
tau_h = 0.01;
C=zeros(1,T);
h=zeros(1,T);
n=zeros(1,T);
beta = 1;
gamma = 0;

for i=1:T-1
    if mod(i,50)==0 && i<T/8, n(i)=1;
    elseif mod(i,100)==0 && i<T/4, n(i)=1;
    elseif mod(i,500)==0 && i<T/2, n(i)=1;
    elseif mod(i,1000)==0, n(i)=1;
    end
    C(i+1) = (1-dt/tau_c)*C(i) + h(i);
%     h(i+1) = (1-dt/tau_h)*h(i) + beta*n(i);
    h(i+1) = (1-dt/tau_h+gamma*n(i))*h(i) + beta*n(i);
end

figure(1), clf, hold on
plot(t,C,'r')
plot(t,h,'b')
plot(t,n,'k')
legend('C','h','n')

% function Ca
% clear, clc
% dt=.01;
% t=0:dt:100;
% T=length(t);
% tau_c = 0.1;
% tau_h = 0.5;
% C=zeros(1,T);
% h=zeros(1,T);
% n=zeros(1,T);
% beta = 1;
% figure(1), clf, hold on
%
% for i=1:T-1
%     if mod(i,50)==0 && i<T/8, sub
%     elseif mod(i,100)==0 && i<T/4, sub
%     elseif mod(i,500)==0 && i<T/2, sub
%     elseif mod(i,1000)==0, sub,  end
%     C(i+1)= (1-dt/tau_c)*C(i) + beta*n(i);
%     %     C(i+1)= (1-dt/tau_c)*C(i) + h(i);
%     %     h(i+1)= (1-dt/tau_h)*h(i) + beta*n(i);
% end
%
%     function sub
%         n(i)=1;
%         plot(t(1:i),C(1:i),'r')
%         plot(t(1:i),h(1:i),'b')
%         stem(t(1:i),n(1:i),'k')
%         ind=find(n,2,'last');
%         h(ind(1)+1)-h(i)
%         keyboard
%     end
%
% end

% function Ca
% clear, clc
% dt=.01;
% t=0:dt:100;
% T=length(t);
% tau_c = 2;
% tau_h = 0.5;
% C=zeros(1,T);
% h=zeros(1,T);
% n=zeros(1,T);
% beta = 1;
% figure(1), clf, hold on
%
% for i=1:T-1
%     if mod(i,10)==0 && i<T/8, n(i)=1;
%     elseif mod(i,100)==0 && i<T/4, n(i)=1;
%     elseif mod(i,500)==0 && i<T/2, n(i)=1;
%     elseif mod(i,1000)==0, n(i)=1;  end
%     C(i+1)= (1-dt/tau_c)*C(i) + beta*n(i);
%     %     C(i+1)= (1-dt/tau_c)*C(i) + h(i);
%     %     h(i+1)= (1-dt/tau_h)*h(i) + beta*n(i);
% end
%
%     function sub
%         n(i)=1;
%         plot(t(1:i),C(1:i),'r')
%         plot(t(1:i),h(1:i),'b')
%         stem(t(1:i),n(1:i),'k')
%         ind=find(n,2,'last');
%         h(ind(1)+1)-h(i)
%         keyboard
%     end
%
% figure(1), clf, hold on
% plot(t,C,'r')
% plot(t,h,'b')
% plot(t,n,'k')
% legend('C','h','n')
%
% end