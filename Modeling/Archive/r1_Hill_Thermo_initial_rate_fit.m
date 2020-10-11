function Prediction=r1_Hill_Thermo_initial_rate_fit(x0, x, BcdData, RuntData)

%Gives the Prediction in order to do a fit with
%lsqnonlin

%Starting conditions
r_basal = x(1);
r = x(2);
K_a = x(3);
N =x(4);

K_r = x0(1);
% Assumption
x0(2) = 0;
r_R1 = x0(2);

Numerator = (r_basal*ones(size(BcdData)) + r*(BcdData./K_a).^N +r_R1*(BcdData./K_a).^N.*RuntData./K_r);
Denominator = (1 + (BcdData./K_a).^N + RuntData./K_r + (BcdData./K_a).^N.*RuntData./K_r);

Prediction = Numerator ./ Denominator;
%Prediction = Rate_r0 ./ (1 + RuntData./K_r)

end