function Prediction=rate_r0_Hill_V3(x0, Act)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

% Parameters
K_a = x0(1);
r = x0(2);
%p = x0(2);
%wAP = x0(3);
%N = x0(4);
N=x0(3);

Prediction = r *(Act./K_a).^N  ./...
                (1 + (Act./K_a).^N);

end