function Prediction=r3_Hill_initial_rate_fit(x0, x, BcdData, RuntData)

%Gives the Prediction in order to do a fit with
%lsqnonlin

%Starting conditions
r_basal = x(1);
r = x(2);
K_a = x(3);
N = x(4);

K_r = x0(1);
r_R1 = 0; %x0(2);
r_R2 = 0; %x0(3);
w_R = x0(4);
r_R3 = 0; %x0(5);

Numerator = (r_basal*ones(size(BcdData)) + r*(BcdData./K_a).^N + ...
                3*r_R1*(BcdData./K_a).^N.*RuntData./K_r +...
                3*r_R2*(BcdData./K_a).^N.*(RuntData./K_r).^2*w_R + ...
                r_R3*(BcdData./K_a).^N.*(RuntData./K_r).^3*w_R^2);
            
Denominator = (1 + (BcdData./K_a).^N + 3*RuntData./K_r + 3*(RuntData./K_r).^2*w_R + (RuntData./K_r).^3*w_R^2 +...
                3*(BcdData./K_a).^N.*RuntData./K_r + 3*(BcdData./K_a).^N.*(RuntData./K_r).^2*w_R + ...
                (BcdData./K_a).^N.*(RuntData./K_r).^3*w_R^2);

Prediction = Numerator ./ Denominator;

end