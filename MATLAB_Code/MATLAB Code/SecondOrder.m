function [zeta, wn] = SecondOrder(ts,Mp)
A = (log(Mp/100))^2;
zeta = sqrt(A/(A+pi^2)); 
wn = 4/(zeta*ts); 
end

