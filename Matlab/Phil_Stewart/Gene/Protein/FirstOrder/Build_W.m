% Function to build the system of linear equations from discretizing the
% equation for w (GFP), the matrix is tridiagonal with size Nx
%  no-flux boundary condition at x = 0 and x = 1

% Input: L -> current biofilm thickness (scalar)
%        L_dot -> the time derivative of current biofilm thickness (scalar)
%        S -> substrate concentration (vector of length Nx) at current time
%             step
%        Dw -> diffusion coefficient of w 
%        Wn -> the solution (chemical species concentration) at current
%              time step n (vector of length Nx)
%        fi -> inducer concentration
%        u  -> velocity due to growth, vector of length Nx+1 (onsite grid)
%        flag -> 1: with diffusion, 2: without diffusion


% Output: d -> the main diagonal of the matrix (vector of length Nx)
%         e -> the subdiagonal of the matrix (vector of length Nx - 1)
%         f -> the superdiagonal of the matrix (vector of length Nx - 1)
%         rhs -> the right hand side of the linear system (vector of length
%         Nx)


function [d,e,f,rhs] = Build_W(L, L_dot, S, Wn, fi, u, flag)

global Alpha Mu_max dx dt Nx Dw % dx = 1/Nx

x = 0:dx:1;  % the coordinates of spatial grid points in the scaled domain

n = Nx;  % size of the tridiagonal linear system

d = zeros(n, 1);
e1 = zeros(n, 1);
f1 = zeros(n, 1);
rhs = zeros(n,1);

%sv = sign(L_dot); % Sign of the convection term, used for upwind scheme

if flag == 1 % with diffusion
    
    eta = .5*Dw*dt/(dx*dx*L*L);
    
    for j = 1:n
        
        d(j) = 1 + 2*eta + (.5*dt/dx)*(u(j+1) - u(j));  % the main diagonal
        
        e1(j) = -eta - .5*dt*u(j)/dx + 0.5*(dt/dx)*(x(j) + .5*dx)*(L_dot/L); % the subdiagonal
        
        f1(j) = -eta + .5*dt*u(j+1)/dx - 0.5*(dt/dx)*(x(j) + .5*dx)*(L_dot/L); % the superdiagonal
        
        if j == 1  %  no-flux BC at x = 0 (j = 1)
            
            rhs(j) = eta*Wn(j+1) + (1 - 2*eta)*Wn(j) + eta*Wn(j+1) + dt*Alpha*Mu_max*fi*S(j);
            
        elseif j == n    % no-flux BC at x = 1 (j = n)
            
            rhs(j) = eta*Wn(j-1) + (1 - 2*eta)*Wn(j) + eta*Wn(j-1) + dt*Alpha*Mu_max*fi*S(j);
            
        else
            
            rhs(j) = eta*Wn(j-1) + (1 - 2*eta)*Wn(j) + eta*Wn(j+1) + dt*Alpha*Mu_max*fi*S(j);
            
        end
        
    end
    
    f1(1) = f1(1) + e1(1);  %  impose the natural BC at x = 0 (j = 1)
    
    e1(n) = e1(n) + f1(n);  %  impose the natural BC at x = 1 (j = n)
    
elseif flag == 2 % without diffusion
    
    for j = 2:n-1
        
        
        d(j) = 1 + (.5*dt/dx)*(u(j+1) - u(j));  % the main diagonal
        
        e1(j) =  - .5*dt*u(j)/dx + 0.5*(dt/dx)*(x(j) + .5*dx)*(L_dot/L); % the subdiagonal
        
        f1(j) =  .5*dt*u(j+1)/dx - 0.5*(dt/dx)*(x(j) + .5*dx)*(L_dot/L); % the superdiagonal
        
        rhs(j) = Wn(j) + dt*Alpha*Mu_max*fi*S(j);
        
    end
    
    d(1) = 1 + dt*Mu_max*S(1);
    rhs(1) = Wn(1) + dt*Alpha*Mu_max*fi*S(1);
    
    d(n) = 1 + dt*Mu_max*S(n);
    rhs(n) = Wn(n) + dt*Alpha*Mu_max*fi*S(n);
    
else
    
    disp('flag must be 1 or 2 in Build_W.m');
    
end
    

e = e1(2:n);

f = f1(1:n-1);

end



