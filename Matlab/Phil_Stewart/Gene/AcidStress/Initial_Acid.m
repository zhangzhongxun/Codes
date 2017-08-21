% Function to set the initial condition for GFP concentration and the biofilm thickness
% 
% Output: L -> initial thickness of the biofilm, a dimensionless scalar
%         W -> initial gene concentration
%         P -> initial lactate concentration
%         S -> initial glucose
%         u -> velocity
%         length Nx, all are dimensionaless
function [L, W, P, S, u] = Initial_Acid()

  global Nx 
  
  L = 0.4; % 0.35; %1.39;
  
  W = zeros(Nx+1,1);
  
  P =  zeros(Nx+1,1); 
  
  S = zeros(Nx+1,1);
  
  u = zeros(Nx+1,1);
  
end
  
  