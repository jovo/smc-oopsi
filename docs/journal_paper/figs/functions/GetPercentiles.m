function ptiles = GetPercentiles(p,w,x)

sz      = size(x);
T       = sz(2);
n       = length(p);
ptiles  = zeros(length(p),T);

for t=1:T;
    [s ind]    = sort(x(:,t));
    cumw        = cumsum(w(ind,t));
    for i=1:n
        ptiles(i,t) = s(find(cumw>p(i),1));
    end
end