% Function to set the initial condition for GFP concentration and the biofilm thickness
% 
% Output: L -> initial thickness of the biofilm, a dimensionless scalar
%         w -> initial GFP concentration, vector of
%         s -> nutrient concentration
%         u -> velocity
%         length Nx, all are dimensionaless
function [L, w, s, u] = Initial_GFP()

  global Nx;
  
  L = 1.39;
  
  w = zeros(Nx,1);
  
  s = zeros(Nx,1);
  
  u = zeros(Nx+1,1);
  
end
  
  