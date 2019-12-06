function Prediction=r0_Hill_initial_rate_fit(x0, BcdData)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
f_basal = x0(1);
f = x0(2);
K_a = x0(3);
%N = x0(4);
N=6;

Prediction = f_basal*ones(size(BcdData)) + (BcdData./K_a).^N  * f ./ (1 + (BcdData./K_a).^N);

end