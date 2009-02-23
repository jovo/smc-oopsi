function [no y] = binocheck(x,n,p)
        
   k = find(x >= 0  &  x == round(x)  &  x <= n);  
   nk = gammaln(n(k) + 1) - gammaln(x(k) + 1) - gammaln(n(k) - x(k) + 1);
   y = nk + x(k).*log( p(k)) + (n(k) - x(k)).*log(1 - p(k));
        
   nk = gammaln(n + 1) - gammaln(x + 1) - gammaln(n - x + 1);
   no = nk + x.*log( p) + (n - x).*log(1 - p);
