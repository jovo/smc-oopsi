clc
fprintf('\nbackward step..........')
for t=12000:-1:1
    if mod(t,100)==0 && t>=9900
        fprintf('\b\b\b\b\b%d',t)
    elseif mod(t,100)==0 && t>=900
        fprintf('\b\b\b\b%d',t)
    elseif mod(t,100)==0 
        fprintf('\b\b\b%d',t)
    end
end