function Prediction=rate_r1_Hill_V3(x0, y0, Act, Rep)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

% Parameters
K_a = x0(1);
r = x0(2);
N=x0(3);

K_r = y0(1);
w_AR = y0(2);

Numerator =  (Act/K_a).^N + (Act/K_a).^N .*(Rep/K_r)*w_AR ;
Denominator = 1 + (Act/K_a).^N + (Rep/K_r) +(Act/K_a).^N .*(Rep/K_r)*w_AR;
Pbound = Numerator ./ Denominator;
Prediction = r*Pbound;

end