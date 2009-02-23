function finv = InvGenHill(B,F)

finv    = ((B.k_d*(F-B.beta))./(B.alpha-F+B.beta)).^(1/B.n); %initialize search with f^{-1}(o)
