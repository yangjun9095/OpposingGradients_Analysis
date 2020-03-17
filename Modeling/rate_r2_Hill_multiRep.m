function Prediction=rate_r2_Hill_multiRep(x0, y0, Act, Rep)
%% A funciton to predict the Txn output (initial rate)
% with previously estimated parameters : 
% Also, with an assumption that 

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

%Starting conditions
K_a = x0(1);
p = x0(2);
wAP = x0(3);
r = x0(4);
%N = x0(4);
N=6;
K_r = y0(1);
wRP1 = y0(2);
wRP2 = y0(3);

% New parameters (repressor cooperativity?)

Pbound =  (p + p*(Act./K_a).^N  *wAP + ...
                p*(Rep./K_r)*(wRP1+wRP2) + p*(Rep./K_r).^2*wRP1*wRP2 +...
                p*(Act./K_a).^N .*(Rep./K_r)*(wRP1+wRP2)*wAP +...
                p*(Act./K_a).^N .*(Rep./K_r).^2*wRP1*wRP2*wAP)./...
          (1 + (Act./K_a).^N + p + 2*(Rep./K_r)+ (Rep./K_r).^2 +...
          2*(Act./K_a).^N .*(Rep./K_r)+ (Act./K_a).^N .*(Rep./K_r).^2+...
           p*(Act./K_a).^N*wAP + p*(Rep./K_r)*(wRP1+wRP2) + p*(Rep./K_r).^2*wRP1*wRP2+...
           p*(Act./K_a).^N .*(Rep./K_r)*(wRP1+wRP2)*wAP + ...
           p*(Act./K_a).^N .*(Rep./K_r).^2*wRP1*wRP2*wAP);
Prediction = r*Pbound;

end