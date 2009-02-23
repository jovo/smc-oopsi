function MCvar = GetMCvar(T,N,X,w)

k=10;
n=100;
U_resamp  = rand(T,n,N); %resampling matrix
MCmean = zeros(n,T);
MCvar  = zeros(1,T);
i=1;
for t=1:T
    for i=1:n
        [foo,ind]   = histc(U_resamp(t,:,i),[0  cumsum(w(:,t))']);
        MCmean(i,t) = mean(X(ind+(t-1)*N));
    end
    MCvar(t)=var(MCmean(:,t));
end