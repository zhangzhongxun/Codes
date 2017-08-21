%function to interpolate a vector by linear interpolation
%Input : u -> original vector, n -> number of interpolation points between
%two data points
%Output: uI -> interpolated vector
function uI = L_InterP(u, n)

   N = length(u);
   NI = N + (N-1)*n;
   
   I = 1:N;
   
   IS = 1: n+1 : NI; % index corresponding to original value
   
   IL = IS(1:N-1);
   
   IR = IS(2:N); 
   
   uI = zeros(1,NI);
   
   uI(IS) = u;
   
   for j = 1:n
       
       uI(IL + j) = u(1:N-1) + (u(2:N) - u(1:N-1))*(j/(n+1));
       
   end
   
end