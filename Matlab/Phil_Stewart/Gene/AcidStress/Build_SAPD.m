% Function to build the system of linear equations from discretizing the
% equation for s (nutrient), the matrix is tridiagonal with size Nx
%  no-flux boundary condition at x = 0 and Dirichlet BC at x = 1
%       
% Input: r -> reaction term, vector of length Nx, 



% Output: d -> the main diagonal of the matrix (vector of length Nx)
%         e -> the subdiagonal of the matrix (vector of length Nx - 1)
%         f -> the superdiagonal of the matrix (vector of length Nx - 1)
%         rhs -> the right hand side of the linear system (vector of length
%         Nx)


function [d,e,f,rhs] = Build_SAPD(r)

global  dx  Nx % dx = 1/Nx

x = 0:dx:1;  % the coordinates of spatial grid points in the scaled domain

x_off = x(1:Nx) + .5*dx; % offset grid points where s is defined

n = Nx;  % size of the tridiagonal linear system

d = zeros(n, 1);
e1 = zeros(n, 1);
f1 = zeros(n, 1);
rhs = r;

d(2:n) = -((x_off(1:n-1).^2 + x_off(2:n).^2)./x(2:n).^2)/(dx*dx);
e1(2:n) = (x_off(1:n-1).^2)./(x(2:n).^2)/(dx*dx);
f1(2:n) = (x_off(2:n).^2)./(x(2:n).^2)/(dx*dx);

%  impose the natural BC at x = 0 (j = 0)
d(1) = -6/(dx*dx);
f1(1) = 6/(dx*dx);

%  impose the Dirichlet BC at x = 1 (j = n)
rhs(n) = rhs(n) - f1(n)*1;

e = e1(2:n);

f = f1(1:n-1);

end