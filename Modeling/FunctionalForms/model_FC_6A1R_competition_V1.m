function [FC, Z_b, Z_bp, Z_br, Z_brp] =...
                    model_6A1R_competition_V1(Bcd, Runt,params) 

%% Definition of parameters
% [Kb, Kr, w_a, w_ap, w_ar, p] = params;
Kb = params(1);
Kr = params(2);
w_b = params(3);
w_bp = params(4);
w_br = params(5);
% w_rp = 1;
p = params(6);
% R_max = params(8);

b = Bcd./Kb;
r = Runt./Kr;

% Calculate the partition function
Z_b = (1-1/w_b) + 1/w_b * (1+w_b*b).^6;
Z_bp = p*(1-1/w_b) + p/w_b*(1+w_b*b*w_bp).^6;
Z_br = r*(1-1/w_b) + r./w_b .*(1+w_br*w_b*b).^6;
Z_brp = r*p*(1-1/w_b) + r*p/w_b.*(1+w_br*w_b*b*w_bp).^6;

% Calculate the P_bound
P_bound = (Z_bp + Z_brp)./ (Z_b + Z_bp + Z_br + Z_brp);

P_bound_null = (Z_bp)./(Z_b + Z_bp);

FC = P_bound./P_bound_null;
end