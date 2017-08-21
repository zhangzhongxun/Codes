% Function to set the initial condition for GFP concentration and the biofilm thickness
% 
% Output: L -> initial thickness of the biofilm, a dimensionless scalar
%         WL -> initial lsaI concentration
%         WR -> initial rsaL concentration
%         QS -> Quorum sensing molecule concentration
%         u -> velocity
%         length Nx, all are dimensionaless
function [L, WL, WR, QS, S u] = Initial_QS()

  global dx Nx Q_BC
  
  x = .5*dx:dx:1-.5*dx;
  
  L = 0.1; % 0.35; %1.39;
  
  WL = zeros(Nx,1);
  
  WR = zeros(Nx,1);
  
  QS =  zeros(Nx,1); %Q_BC*cos(pi*x); 
  
  S = zeros(Nx,1);
  
  u = zeros(Nx+1,1);
  
end
  
  