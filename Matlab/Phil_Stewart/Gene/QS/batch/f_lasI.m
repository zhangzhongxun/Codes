% function for the growth of lsaI due to the presence of QS moledule

function y = f_lasI(Q, max_rate, Q_b)

   k = (0.1 - 0.001)/Q_b;
   
   pos1 = (Q <= Q_b);
   pos2 = (Q > Q_b);
   
   y = zeros(size(Q));
   y(pos1) = k*Q(pos1) + 0.001;
   y(pos2) = 0.1;
   
   y = max_rate*y;
   
end
   
   