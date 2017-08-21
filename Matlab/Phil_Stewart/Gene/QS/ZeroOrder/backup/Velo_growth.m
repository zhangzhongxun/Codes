% Function to calculate the rate of change of the biofilm thickness (L_dot)
% and the velocity due to growth
% Input:  L -> current thickness of the biofilm 
%         S -> substrate concentration at current time
%         
% Output: L_dot -> the rate of change of L, a scalar
%         u  -> velocity, vector of length Nx+1

function [L_dot, u] = Velo_growth(L, S, x0)

global Nx  dx  Mu_max Sc Sigma

x = 0:dx:1;

pos = find(x > x0);

u = zeros(Nx+1,1);

s_u = (S > Sc);

% Integral of the growth rate over spatial domain which gives the
% velocity
% 
% for k = 2:Nx+1
%     
%     u(k) = u(k-1) + Mu_max*s_u(k-1)*dx;
%     
% end

u(pos) = Mu_max*(x(pos) - max(0,x0));

% The rate of change of the biofilm thickness, combining detachment and
% growth
L_dot =  -Sigma*L*L + L*u(Nx+1);
 

end