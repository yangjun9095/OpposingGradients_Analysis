% P_bound for 1A1R (r=0) for fitting
function P_bound=P_bound_1A1R_NoRepressor(x0)

%Starting conditions
p = x0(1);
w_ap = x0(2);
Ka= x0(3);  % scaling factor for Activator

w_rp = 0.01; % Useless in this case, since r=0 already.
% Default
a = linspace(10^(-2),10^2, 10^4);
r=0;

a = a/Ka;

Denominator= p + a.*p.*w_ap + r.*p.*w_rp +...
                r.*a.*p.*w_rp.*w_ap;
Numerator = 1 + p + a + r + a.*p.*w_ap + p.*r.*w_rp + ...
                a.*r + r.*a.*p.*w_rp.*w_ap;

P_bound = Denominator./Numerator;

end