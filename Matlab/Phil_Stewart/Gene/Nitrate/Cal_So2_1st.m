% function to calculate the oxyten concentration with 1st order kinetics
% with analytical formula for the denitrification study
% So2 = 1 (scaled bulk concentration )at x = x1, and d So2/dx = 0 at x = x2
% here 0 <= x <= x1, and x2 < x1, 
% D is the dimensionless diffusion coefficient
% dx is the spatial grid size
% So2 is the dimensionless oxygen concentration

function So2 = Cal_So2_1st(D, x1, x2,  dx)

x = 0:dx:x1;

a11 = exp(x1/sqrt(D));
a12 = exp(-x1/sqrt(D));
a21 = exp(x2/sqrt(D))/sqrt(D);
a22 = -exp(-x2/sqrt(D))/sqrt(D);

detA = a11*a22 - a21*a12;

c1 = a22/detA;

c2 = -a21/detA;

So2 = c1*exp(x/sqrt(D)) + c2*exp(-x/sqrt(D));

% pos = (x < x2);
% pos2 = (x == x2);
% 
% So2(pos) = exp(x(pos) - x2)*So2(pos2);

end



