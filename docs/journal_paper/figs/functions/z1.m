function x = z1(y)
% linearly normalize to between 0 and 1
x = (y-min(y))/(max(y)-min(y))+eps;