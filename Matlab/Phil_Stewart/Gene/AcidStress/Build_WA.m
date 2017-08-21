% Function to build the system of linear equations from discretizing the
% equation for W (gene), the matrix is tridiagonal with size Nx
%  no-flux boundary condition at x = 0 and Dirichlet BC at x = 1
%       
% Input: r -> reaction term, vector of length Nx, 
%        mu -> biofilm growth rate at current time 
%        wn -> gene concentration at current time step n
%        L_dot -> rate of change of biofilm thickness
%        L   -> biofilm thickness
%        u   -> velocity at current time step


% Output: d -> the main diagonal of the matrix (vector of length Nx)
%         e -> the subdiagonal of the matrix (vector of length Nx - 1)
%         f -> the superdiagonal of the matrix (vector of length Nx - 1)
%         rhs -> the right hand side of the linear system (vector of length
%         Nx)


function [d,e,f,rhs] = Build_WA(r, mu, wn, L_dot, L, u)

global  dx  dt Nx t0  Sigma Kg flag_kg% dx = 1/Nx

x = 0:dx:1;  % the coordinates of spatial grid points in the scaled domain

%x_off = x(1:Nx) + .5*dx; % offset grid points where s is defined

n = Nx+1;  % size of the tridiagonal linear system

d = zeros(n, 1);
e1 = zeros(n, 1);
f1 = zeros(n, 1);
rhs = r;

sig = sign(-L_dot);  % used for upwind scheme

d(2:n-1) = 1 + (dt/dx)*sig*(-x(2:n-1)*L_dot/L) + (dt/dx)*u(2:n-1) + flag_kg*Kg*dt;
e1(2:n-1) = (dt/dx)*((1 + sig)/2)*(L_dot/L)*x(2:n-1) - (dt/dx)*(x(1:n-2).^2./x(2:n-1).^2).*u(1:n-2);
f1(2:n-1) = (dt/dx)*((1 - sig)/2)*(-L_dot/L)*x(2:n-1);
rhs(2:n-1) = wn(2:n-1) + dt*r(2:n-1);

%  impose the BC at x = 0 (j = 0)
d(1) = 1 + dt*t0*mu(1) + flag_kg*Kg*dt;
rhs(1) = wn(1) + dt*r(1);

%  impose the  BC at x = 1 (j = n)
d(n) = 1 + dt*t0*mu(n) + (dt/dx)*Sigma*L + flag_kg*Kg*dt;
e1(n) = -(dt/dx)*Sigma*L;
rhs(n) = wn(n) + dt*r(n);


e = e1(2:n);

f = f1(1:n-1);

end