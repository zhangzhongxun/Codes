% Function to build the system of linear equations from discretizing the
% equation for s (nutrient), the matrix is tridiagonal with size Nx
%  no-flux boundary condition at x = 0 and x = 1
%        Ds -> diffusion coefficient of nutrient
% Input: L -> current biofilm thickness (scalar)
%        s0 -> Dirichlet BC at x = 1



% Output: d -> the main diagonal of the matrix (vector of length Nx)
%         e -> the subdiagonal of the matrix (vector of length Nx - 1)
%         f -> the superdiagonal of the matrix (vector of length Nx - 1)
%         rhs -> the right hand side of the linear system (vector of length
%         Nx)


function [d,e,f,rhs] = Build_S(L, s0)

global  dx  Nx  Ds % dx = 1/Nx

x = 0:dx:1;  % the coordinates of spatial grid points in the scaled domain

n = Nx;  % size of the tridiagonal linear system

%sv = sign(L_dot); % Sign of the convection term, used for upwind scheme

eta = Ds/(dx*dx*L*L);

d = (1+2*eta)*ones(n, 1);
e1 = -eta*ones(n, 1);
f1 = -eta*ones(n, 1);
rhs = zeros(n,1);



f1(1) = -2*eta;  %  impose the natural BC at x = 0 (j = 1)

d(n) = 1 + 4*eta; 

e1(n) = -4*eta/3;  %  impose the Dirichlet BC at x = 1 (j = n)

rhs(n) = (8/3)*eta*s0;

e = e1(2:n);

f = f1(1:n-1);


end