function Prediction=rate_r2_prediction(x0, y1, y2, Act, Rep)

%Gives the chi square of the data to the fit form in order to do a fit with
%lsqnonlin

% Parameters
K_a = x0(1);
r = x0(2);
N=x0(3);

% repression
K_r1 = y1(1);
w_AR1 = y1(2);

K_r2 = y2(1);
w_AR2 = y2(2);

Numerator =  (Act/K_a).^N + (Act/K_a).^N .*(Rep/K_r1)*w_AR1 + ...
             (Act/K_a).^N .*(Rep/K_r2)*w_AR2 + (Act/K_a).^N .*(Rep/K_r1)*w_AR1.*(Rep/K_r2)*w_AR2;
Denominator = 1 + (Act/K_a).^N + (Rep/K_r1) + (Rep/K_r2) + (Rep/K_r1).*(Rep/K_r2) +...
                (Act/K_a).^N .*(Rep/K_r1)*w_AR1 + (Act/K_a).^N .*(Rep/K_r2)*w_AR2 +...
                (Act/K_a).^N .*(Rep/K_r1)*w_AR1.*(Rep/K_r2)*w_AR2;
Pbound = Numerator ./ Denominator;
Prediction = r*Pbound;


end