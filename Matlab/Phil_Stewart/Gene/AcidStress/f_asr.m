%function to calculate the growth factor for gene
% Input: p -> lactate concentration, dimensionless (mg/L)
%        
% Output: f -> growth factor of gene (unitless)
function f = f_asr(p) 

  global p0   % p0 -> reference lactate concentration (mg/L)

  f = ones(size(p));
  
  pos = find(p < 500/p0);
  
  f(pos) = (p(pos)/(500/p0)).^3;
  
end