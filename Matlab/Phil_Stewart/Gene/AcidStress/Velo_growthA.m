% Function to calculate the rate of change of the biofilm thickness (L_dot)
% and the velocity due to growth
% Input:  L -> current thickness of the biofilm 
%         mu -> growth rate of biofilm
%         
% Output: L_dot -> the rate of change of L, a scalar
%         u  -> velocity, vector of length Nx+1

function [L_dot, u] = Velo_growthA(L, mu)

global Nx  dx  Sigma t0

u = zeros(1,Nx+1);

x = 0:dx:1;

% Integral of the growth rate over spatial domain which gives the
% velocity

for k = 2:Nx+1
    
    u(k) = (x(k-1)^2*u(k-1) + .5*dx*t0*(mu(k-1)*x(k-1)^2 + mu(k)*x(k)^2))/(x(k)^2);
    
end


% The rate of change of the biofilm thickness, combining detachment and
% growth
L_dot =  -Sigma*L*L + L*u(Nx+1);
 

end