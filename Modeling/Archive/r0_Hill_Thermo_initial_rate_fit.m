function Prediction=r0_Hill_Thermo_initial_rate_fit(x0, BcdData)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
r_basal = x0(1);
r = x0(2);
K_a = x0(3);
%N = x0(4);
N=6; % Let's have N (gMM constant) as a free parameter.

Prediction = (r_basal*ones(size(BcdData)) + r.*(BcdData./K_a).^N )./ (1 + (BcdData./K_a).^N);

end