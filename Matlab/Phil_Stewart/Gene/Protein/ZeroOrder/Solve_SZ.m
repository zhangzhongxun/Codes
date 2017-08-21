% Function to calculate the solution of the differential
% equation for s (nutrient) with zero order kinetics using analytical
% solution (which is a quadratic function of x). 
% Dirichlet BC at x = 1, no flux boundary condition at
% x0 where s(x0) = 0. Here x0 is solved also.
%        
% Input: L -> current biofilm thickness (scalar)
%        s0 -> Dirichlet BC at x = 1

% Output: s -> solution


function s = Solve_SZ(L, s0)

global  dx  Nx  Ds % dx = 1/Nx

x = 0:dx:1;  % the coordinates of spatial grid points in the scaled domain

x_off = x(1:Nx) + .5*dx; % offset grid points where s is defined

D = Ds/(L*L);

b = (-1 + sqrt(2*D))/D;  % coefficient of the linear term of the analytical solution

c = 1 - 1/(2*D) - b;   % coefficient of the constant term of the analytical solution

x0 = 1 - sqrt(2*D); % position where s becomes zero, must have 0 <= x0 <= 1

s = zeros(size(x_off)); % s = 0 for x < x0

pos = find(x_off >= x0);

s(pos) = (1/(2*D))*x_off(pos).^2 + b*x_off(pos) + c;

end