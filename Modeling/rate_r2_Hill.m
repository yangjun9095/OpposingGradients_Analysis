function Prediction=rate_r2_prediction(x0, y1, y2, Act, Rep)

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

% New parameters (repressor cooperativity?)

Pbound =  (p + p*(Act./K_a).^N  *wAP + ...
                2*p*(Rep./K_r)*wRP + p*(Rep./K_r).^2*wRP.^2 +...
                2*p*(Act./K_a).^N .*(Rep./K_r)*wRP*wAP) +...
                p*(Act./K_a).^N .*(Rep./K_r).^2*wRP^2*wAP)./...
          (1 + (Act./K_a).^N + p + 2*(Rep./K_r)+ (Rep./K_r).^2 +...
          2*(Act./K_a).^N .*(Rep./K_r)+ (Act./K_a).^N .*(Rep./K_r).^2+...
           p*(Act./K_a).^N*wAP + 2*p*(Rep./K_r)*wRP + p*(Rep./K_r).^2*wRP^2+...
           2*p*(Act./K_a).^N .*(Rep./K_r)*wRP*wAP + ...
           p*(Act./K_a).^N .*(Rep./K_r).^2*wRP^2*wAP);
Prediction = r*Pbound;

end