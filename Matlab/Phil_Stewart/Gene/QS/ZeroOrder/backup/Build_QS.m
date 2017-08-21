% Function to build the system of linear equations from discretizing the
% equation for QS (Quorum sensing molecue concentration), 
% the matrix is tridiagonal with size Nx
% no-flux boundary condition at x = 0 and zero-Dirichlet at x = 1

% Input: L -> current biofilm thickness (scalar)
%        L_dot -> the time derivative of current biofilm thickness (scalar)
%        S -> substrate concentration (vector of length Nx) at current time
%             step
%        Dq -> diffusion coefficient of QS 
%        Qn -> the solution (chemical species concentration) at current
%              time step n (vector of length Nx)
%        W_L -> concentration of lasI
%        u  -> velocity due to growth, vector of length Nx+1 (onsite grid)
%        Q_BC: the value of Dirichlet BC at x = 1, default is 0


% Output: d -> the main diagonal of the matrix (vector of length Nx)
%         e -> the subdiagonal of the matrix (vector of length Nx - 1)
%         f -> the superdiagonal of the matrix (vector of length Nx - 1)
%         rhs -> the right hand side of the linear system (vector of length
%         Nx)


function [d,e,f,rhs] = Build_QS(L, L_dot, S, W_L, u, Q_BC)

global Alpha_1 Mu_max dx  Nx Dq Sc Kg_Q flag_kg_Q % dx = 1/Nx

x = 0:dx:1;  % the coordinates of spatial grid points in the scaled domain

n = Nx;  % size of the tridiagonal linear system

d = zeros(n, 1);
e1 = zeros(n, 1);
f1 = zeros(n, 1);
rhs = zeros(n,1);


s_u = (S > Sc);  % only has protein production when S > S_c, zero order kinetics


   % eta = .5*Dq*dt/(dx*dx*L*L);
    
    eta = Dq/(dx*dx*L*L);
    
    for j = 1:n

        d(j) = flag_kg_Q*Kg_Q + 2*eta + (.5/dx)*(u(j+1) - u(j));  % the main diagonal

        e1(j) = -eta - .5*u(j)/dx + (0.5/dx)*(x(j) + .5*dx)*(L_dot/L); % the subdiagonal

        f1(j) = -eta + .5*u(j+1)/dx - (0.5/dx)*(x(j) + .5*dx)*(L_dot/L); % the superdiagonal

        rhs(j) = Alpha_1*Mu_max*W_L(j)*s_u(j);

    end
    
    f1(1) = f1(1) + e1(1);  %  impose the natural BC at x = 0 (j = 1)
    
    e1(n) = e1(n) + f1(n)/3;  %  impose the Dirichlet BC at x = 1 (j = n)
    d(n) = d(n) - 2*f1(n);
    rhs(n) = rhs(n) - 8*f1(n)*Q_BC/3;
       

e = e1(2:n);

f = f1(1:n-1);

end



