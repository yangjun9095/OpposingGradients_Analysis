function Prediction=rate_r0_Hill(x0, Act)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions c bv
K_a = x0(1);
p = x0(2);
wAP = x0(3);
r = x0(4);
%N = x0(4);
N=6;

Prediction = r * (p*ones(size(Act)) + p*(Act./K_a).^N  *wAP) ./...
                (1 + (Act./K_a).^N + p + p*(Act./K_a).^N*wAP);

end