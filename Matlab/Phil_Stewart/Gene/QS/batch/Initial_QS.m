% Function to set the initial condition for GFP concentration and the biofilm thickness
% 
% Output: y(1) -> Initial oxyten
%         y(2) -> Initial Bacteria cell
%         y(3) -> Quorum sensing molecule concentration
%         y(4) -> initial lsaI concentration
%         y(5) -> initial rsaL concentration


function y = Initial_QS()

y(1) = 1; % Oxygen

y(2) = .5e-6; % Bacteria cell concentration

y(3:5) =  0;

y = y(:);
  
end
  
  