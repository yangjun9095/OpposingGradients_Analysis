function Prediction=fractionon1_Hill_fit(x0, x, BcdData, RuntData)

%Gives the Prediction in order to do a fit with
%lsqnonlin

%Starting conditions
f_basal = x(1);
f = x(2);
K_a = x(3);
N=6;

K_r = x0(1);
%r_R = x0(2);

fraction_r0 = f_basal*ones(size(BcdData)) +  (BcdData./K_a).^N  * f ./ (1 + (BcdData./K_a).^N);
Prediction = f_basal*ones(size(BcdData)) + fraction_r0 ./ (1 + RuntData./K_r)

% chi2 = (Prediction - InitialRateData).^2;

end