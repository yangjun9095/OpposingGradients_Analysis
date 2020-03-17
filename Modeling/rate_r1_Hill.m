function Prediction=rate_r1_Hill(x0, y0, Act, Rep)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions c bv
K_a = x0(1);
p = x0(2);
wAP = x0(3);
r = x0(4);
%N = x0(4);
N=6;
K_r = y0(1);
wRP = y0(2);

Pbound =  (p*ones(size(Act)) + p*(Act./K_a).^N  *wAP + ...
                p*(Rep./K_r)*wRP + p*(Act./K_a).^N .*(Rep./K_r)*wRP*wAP) ./...
          (1 + (Act./K_a).^N + p + (Rep./K_r)+ (Act./K_a).^N .*(Rep./K_r)+...
                p*(Act./K_a).^N*wAP + p*(Rep./K_r)*wRP + p*(Act./K_a).^N .*(Rep./K_r)*wRP*wAP );
Prediction = r*Pbound

end