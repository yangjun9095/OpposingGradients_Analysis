% P_bound for 2A (two activator binding sites)
function P_bound=P_bound_2A(x0)

%Starting conditions
p = x0(1);
w_ap = x0(2);
Ka= x0(3);  % scaling factor for Activator
w_aa=x0(4);

% Default
a = linspace(10^(-2),10^2, 10^4);
                    
a = a/Ka;

Numerator= p + 2*a.*p.*w_ap + ...
                a.^2.*p.*w_ap^2.*w_aa;
Denominator = 1 + p + 2*a + a.^2.*w_aa + 2*a.*p.*w_ap + ...
                a.^2.*p.*w_ap.^2.*w_aa;

P_bound = Numerator./Denominator;

end