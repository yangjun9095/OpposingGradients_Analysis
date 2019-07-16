function Prediction=r2_Hill_initial_rate_fit(x0, x, BcdData, RuntData)

%Gives the Prediction in order to do a fit with
%lsqnonlin

%Starting conditions
r_basal = x(1);
r = x(2);
K_a = x(3);
N=6;

K_r = x0(1);
%r_R = x0(2);

Rate_r0 = (r_basal*ones(size(BcdData)) + (BcdData./K_a).^N  * r) ./ (1 + (BcdData./K_a).^N);
Prediction = Rate_r0 ./ (1 + (RuntData./K_r).^2);

% chi2 = (Prediction - InitialRateData).^2;

end