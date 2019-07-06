function chi2=r0_Hill_initial_rate_fit(x0, BcdData,InitialRateData)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
r_basal = x0(1);
r = x0(2);
K_a = x0(3);

Prediction = (r_basal*ones(size(BcdData)) + (BcdData./K_a).^6  * r) ./ (1 + (BcdData./K_a).^6);

chi2 = (Prediction - InitialRateData).^2;

end