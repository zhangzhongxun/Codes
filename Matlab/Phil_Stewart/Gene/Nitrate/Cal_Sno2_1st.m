% function to calculate the oxyten concentration with 1st order kinetics
% with analytical formula for the denitrification study
% Sno2 = SB (scaled bulk concentration )at x = x1, and So2 = SIn at x = x2
% here 0 <= x <= x1, and x3 < x2 < x1, 
% at x = x2, there is a flux balance (SB - SIn)*D/(x1 - x2) = D (d Sno2/dx)
% at x = x3, d Sno2/dx = 0 (no diffusion flux)
% D is the dimensionless diffusion coefficient
% dx is the spatial grid size
% Sno2 is the dimensionless oxygen concentration

function Sno2 = Cal_Sno2_1st(D, x1, x2, x3, SB, dx)

x = 0:dx:x1;

% N = length(x);

a11 = exp(x3/sqrt(D))/sqrt(D);
a12 = -exp(-x3/sqrt(D))/sqrt(D);
a21 = exp(x2/sqrt(D)) + (x1-x2)*exp(x2/sqrt(D))/sqrt(D);
a22 = exp(-x2/sqrt(D)) - (x1 - x2)*exp(-x2/sqrt(D))/sqrt(D);

detA = a11*a22 - a21*a12;

c1 = (-SB*a12)/detA;

c2 = (a11*SB)/detA;

Sno2 = c1*exp(x/sqrt(D)) + c2*exp(-x/sqrt(D));

% pos = (x < x2);
% 
% So2(pos) = .1*alpha;

pos2 = (x == x2);

k = (SB - Sno2(pos2))/(x1 - x2);

pos = (x > x2);

Sno2(pos) = SB + k*(x(pos) - x1);

% pos3 = (x == x3);
% pos = (x < x3);
% 
% Sno2(pos) = exp(x(pos) - x3)*Sno2(pos3);

end



