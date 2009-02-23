k=find(Sim.tvec==9.795)
figure, imagesc(S.w_f(:,k-9:k)), colorbar, title('w_f k-9:k')

a=find(S.C(:,k-1)>3.8)
figure
clf, hold on
Sim.tvec=0:Sim.dt:Sim.Nsec-Sim.dt; %time vector
plot(Sim.tvec(k-Sim.frac:k),S.C(:,k-Sim.frac:k)')
plot(Sim.tvec(k-Sim.frac:k),R.O(k-Sim.frac:k),'ok','LineWidth',2)
plot(Sim.tvec(k-Sim.frac:k),R.C(k-Sim.frac:k),'k','LineWidth',3)
plot(Sim.tvec(k-Sim.frac:k),sum(S.w_f(:,k-Sim.frac:k).*S.C(:,k-Sim.frac:k),1),'r','LineWidth',2)
plot(Sim.tvec(k-Sim.frac:k),S.C(a,k-Sim.frac:k)','b','LineWidth',2)

figure, imagesc(S.w_b(:,k-9:k)), colorbar, title('w_b k-9:k')

plot(Sim.tvec(k-Sim.frac:k),S.C(b,k-Sim.frac:k)','g','LineWidth',2)
plot(Sim.tvec(k-Sim.frac:k),S.C(c,k-Sim.frac:k)','r','LineWidth',2)

figure
clf, hold on
Sim.tvec=0:Sim.dt:Sim.Nsec-Sim.dt; %time vector
plot(Sim.tvec(k-10:k+10),S.C(:,k-10:k+10)')
plot(Sim.tvec(k-10:k+10),R.O(k-10:k+10),'ok','LineWidth',2)
plot(Sim.tvec(k-10:k+10),R.C(k-10:k+10),'k','LineWidth',3)
plot(Sim.tvec(k-10:k+10),sum(S.w_f(:,k-10:k+10).*S.C(:,k-10:k+10),1),'r','LineWidth',2)


figure, clf, hold on
plot(S.sig2_I(:,1:k)','r'), hold on, plot(S.sig2_I(:,1:k)','-.k')
plot(S.sig2_c(:,1:k)','r'), hold on, plot(S.sig2_c(:,1:k)','-.k')
plot(S.sig2_o(:,1:k)','r'), hold on, plot(S.sig2_o(:,1:k)','-.k')
plot(S.mu_o(:,1:k)','r'), hold on, plot(S.mu_o(:,1:k)','-.k')
plot(S.mu_I(:,1:k)','r'), hold on, plot(S.mu_I(:,1:k)','-.k')
hold off

figure(3), clf, hold on
plot(S.sig2_I(:,k+1:s)')
plot(S.mu_o(:,k+1:s)','r')
plot(S.sig2_c(:,k+1:s)')
plot(S.sig2_I(:,k+1:s)')
plot(S.mu_I(:,k+1:s)')

figure(3), clf, hold on
plot(S.sig2_I','r')
plot(S.sig2_c','k')
plot(S.sig2_I','b')
plot(S.mu_I','g')

figure, plot(Sim.tvec,S.mu_o','r')
figure, plot(Sim.tvec,S.sig2_o','b')


M.Ibar = sum(S.w_b.*S.I,1);

for k=1:Sim.K
    meanp(k)=S.w_b(S.I(:,k)==1,k)'*S.I(S.I(:,k)==1,k);
end

ln_I    = [log(S.p(:,k)) log(1-S.p(:,k))];          %compute [log(spike) log(no spike)]
mu_ck   = [B.a*S.C(:,k-1)+B.beta B.a*S.C(:,k-1)];   %mean when spiked and not spiked

sig2_zk = repmat((1/B.sig2_c + 1./S.sig2_o(:,k)).^(-1),1,2);        %update sig2_zk
mu_zk   = sig2_zk.*(mu_ck/B.sig2_c + repmat(S.mu_o(:,k)./S.sig2_o(:,k),1,2));                     %update mu_zk
ln_Z    = log(sig2_zk./sqrt(2*pi*repmat(S.sig2_o(:,k),1,2)*B.sig2_c)) + 0.5*(mu_zk.^2./sig2_zk-repmat(S.mu_o(:,k).^2./S.sig2_o(:,k),1,2)-mu_ck.^2/B.sig2_c);  %update ln_z (ignore constant)

sig2_Ik  = (repmat(S.sig2_I(:,k),1,2)+B.sig2_c);
ln_G    = -0.5*log(2*pi*sig2_Ik) -0.5*(repmat(S.mu_o(:,k),1,2) - mu_ck).^2./sig2_Ik;
