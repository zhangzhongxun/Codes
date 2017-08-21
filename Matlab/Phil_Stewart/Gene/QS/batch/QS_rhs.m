% function to build the RHS of the ode system
% y(1) = S -> oxygen
% y(2) = X  -> cell density
% y(3) = Q -> QS molecule
% y(4) = WL -> lasI
% y(5) = WR -> rsaL


function r = QS_rhs(t,y)

global Mu_max Alpha_1 Alpha_2 Alpha_4 max_rate_WL Q_bar g_x g_l g_r Kg_Q flag_kg_Q K_ex

r(1) = -Mu_max*g_x*y(2) + K_ex*(1 - y(1));

r(2) = Mu_max*y(2); %multiply by y(1), oxygen?

r(3) = Alpha_1*Mu_max*y(4)*y(2) - flag_kg_Q*Kg_Q;

r(4) = exp(-Alpha_2*y(5))*f_lasI(y(3), max_rate_WL, Q_bar)*g_l*Mu_max;

r(5) = exp(-Alpha_4*y(1))*f_lasI(y(3), max_rate_WL, Q_bar)*g_r*Mu_max;

r = r(:);

end
