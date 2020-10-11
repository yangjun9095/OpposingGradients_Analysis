function Prediction=r2_Hill_Thermo_initial_rate_fit(x0, x, BcdData, RuntData)

%Gives the Prediction in order to do a fit with
%lsqnonlin

%Starting conditions
r_basal = x(1);
r = x(2);
K_a = x(3);
N = x(4);

K_r = x0(1);
r_R1 = 0; %x0(2);
r_R2 = 0; %x0(2);
w_R = x0(4);

Numerator = (r_basal*ones(size(BcdData)) + r*(BcdData./K_a).^N + ...
                2*r_R1*(BcdData./K_a).^N.*RuntData./K_r + r_R2*(BcdData./K_a).^N.*(RuntData./K_r).^2*w_R);
            
Denominator = (1 + (BcdData./K_a).^N + 2*RuntData./K_r + (RuntData./K_r).^2*w_R +...
                2*(BcdData./K_a).^N.*RuntData./K_r + (BcdData./K_a).^N.*(RuntData./K_r).^2*w_R);

Prediction = Numerator ./ Denominator;
%Prediction = Rate_r0 ./ (1 + RuntData./K_r)


end