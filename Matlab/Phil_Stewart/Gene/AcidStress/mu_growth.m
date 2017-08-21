%function to calculate the dimnesional growth rate mu
% Input: p -> lactate concentration, dimensionless (mg/L)
%        s  -> glucose concentration
%       
% Output: g -> growth rate of biofilm, unit = 1/s

function g = mu_growth(p, s)

 global mu_max p0  % p0 -> reference lactate concentration (mg/L)
 
 g = zeros(size(p));
 
 pos = (p < 3000/p0); % positive growth only when glucose is positive and lactate is less than 3000 mg/l
 
 g(pos) = mu_max*s(pos).*(1 - p(pos)/(3000/p0));
 
end

