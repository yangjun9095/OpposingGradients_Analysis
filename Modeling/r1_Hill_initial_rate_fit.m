function chi2=r1_Hill_initial_rate_fit(x0, x, BcdData, RuntData, InitialRateData)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
r_basal = x(1);
r = x(2);
K_a = x(3);

K_r = x0(1);
r_R = x0(2);

numerator = (r_basal*ones(size(BcdData)) + (BcdData./K_a).^6  * r + (RuntData./K_r).*(BcdData./K_a).^6*r_R);
denominator =  (1 + (BcdData./K_a).^6 + (RuntData./K_r).*(BcdData./K_a).^6 + (RuntData./K_r) );

Prediction = numerator ./ denominator;

chi2 = (Prediction - InitialRateData).^2;

end