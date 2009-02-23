%% generate smoothing kernel and compute errors between true and estimated
% spike trains
gauss_smooth    = exp(-((linspace(-5,5,Sim.T)/.1).^2)); %smoothing gaussian kernel
smooth_real_n   = conv2(R.n,gauss_smooth,'same');       %tru spike train smoothed

% get spikes only using dF/F and threshold
dF              = diff(R.F);                            %get dF
[sdF ind]       = sort(dF);                             %sort dF
sp_ind          = ind(round(Sim.T*sp_thres):end);       %find top (1-sp_thres)% of spikes
thr_n           = zeros(1,Sim.T);                       %extize thrial guess of where spikes are
thr_n(sp_ind)   = 1;                                    %make vector of 1's where spikes are

% smooth thresholded spike train and compute error
smooth_thr_n    = conv2(thr_n,gauss_smooth,'same');     %smooth spike train with gaussian
thr_errs        = (smooth_real_n-smooth_thr_n).^2;      %compute error between real and thr spike train
thr_mean        = mean(thr_errs);                       %compute mse
thr_var         = var(thr_errs)/20;                     %compute var

% smooth ppr estimated spike train and compute error
smooth_ppr_n    = conv2(O.ppr_n,gauss_smooth,'same');
ppr_errs        = (smooth_real_n-smooth_ppr_n).^2;
ppr_mean        = mean(ppr_errs);                       %compute mse
ppr_var         = var(ppr_errs)/20;                     %compute var

%% plot actual and inferred spikes and smoothed trains
figure(1), clf, nrows1=6; gray=[0.75 0.75 0.75];

% calcium
subplot(nrows1,1,1), hold on,
plot(Sim.tvec, R.C,'color',gray)
ppr_ca=conv(O.A*exp(-Sim.tvec/O.tau),O.ppr_n);    %convolve exponential with spike times
plot(Sim.tvec, ppr_ca(1:Sim.T),'g')
ylabel('[Ca^2^+]_t'), axis('tight'), %legend('real','ppr','Location','EastOutside')

%fluorescence
subplot(nrows1,1,2), hold on,
plot(Sim.tvec(Sim.freq:Sim.freq:Sim.T),R.F(Sim.freq:Sim.freq:Sim.T),'.','color','k'), ylabel('F'), axis('tight')
plot(Sim.tvec,R.F-ppr_ca(1:Sim.T),'g'), %legend('F','resid','Location','EastOutside')

%dF and threshold
dF          = diff(R.F);    %get dF
[sdF ind]   = sort(dF);     %sort dF
subplot(nrows1,1,3), hold on, plot(Sim.tvec(Sim.freq:Sim.freq:Sim.T-1),dF(Sim.freq:Sim.freq:Sim.T-1),'.-','color','k'), ylabel('dF'), axis('tight')
subplot(nrows1,1,3), plot(Sim.tvec,sdF(round(Sim.T*sp_thres)),'k'), hold off, axis('tight')

%actual and smoothed spike train
subplot(nrows1,1,4), hold on, bar(Sim.tvec,R.n,'FaceColor',gray,'EdgeColor',gray), ylabel((sum(R.n))),
subplot(nrows1,1,4), plot(Sim.tvec, smooth_real_n,'color',gray), hold off, axis('tight')
title('tru spikes and smoothed')

%dF/F and threshold guess of spike train (and smoothed)
subplot(nrows1,1,5), hold on, bar(Sim.tvec,thr_n,'FaceColor','r','EdgeColor','r'), ylabel((sum(thr_n)))
subplot(nrows1,1,5), plot(Sim.tvec, smooth_thr_n,'r'), hold off, axis('tight')
title(['dF/F and threshold spikes and mse of smoothed = ', num2str(thr_mean), '+/-', num2str(thr_var)])

% %dF/F and threshold guess of spike train (and smoothed)
% subplot(nrows1,1,5), hold on, bar(Sim.tvec,O3.ppr_n,'FaceColor','r','EdgeColor','r'), ylabel((sum(O3.ppr_n)))
% subplot(nrows1,1,5), plot(Sim.tvec, smooth_ppr_n3,'r'), hold off, axis('tight')
% title(['dF/F and threshold spikes and mse of smoothed = ', num2str(ppr_mean3), '+/-', num2str(ppr_var3)])

%ppr guess of spike train (and smoothed)
subplot(nrows1,1,6), hold on, bar(Sim.tvec,O.ppr_n,'FaceColor','g','EdgeColor','g'), ylabel((sum(O.ppr_n)))
subplot(nrows1,1,6), plot(Sim.tvec, smooth_ppr_n,'g'), hold off, axis('tight')
title(['projection pursuit spikes and mse of smoothed = ', num2str(ppr_mean), '+/-', num2str(ppr_var)])

% %residuals
% subplot(nrows1,1,6), plot(Sim.tvec,resid(1:Sim.T)), title('residuals')
% 
%% plot real and estimated kernel
figure(2), clf, hold on,
plot(Sim.tvec, O.A*exp(-Sim.tvec/O.tau),'g')
plot(Sim.tvec,P.A*exp(-Sim.tvec./P.tau_c),'k')
hold off, legend('ppr','real')
