function [output, Z_b, Z_bp, Z_br, Z_brp] = ...
                model_FC_6A1R_combination_all_V1(Bcd, Runt,params) 

%% Definition of parameters
% [Kb, Kr, w_b, w_bp, w_br, w_rp, w_brp, p] = params;
Kb = params(1);
Kr = params(2);
w_b = params(3);
w_bp = params(4);

% repression
w_br = params(5);
w_rp = params(6);
w_brp = params(7);

% basal parameters
p = params(8);
% R_max = params(9);

b = Bcd./Kb;
r = Runt./Kr;

% Calculate the partition function
Z_b = (1-1/w_b) + 1/w_b * (1+w_b*b).^6;
Z_bp = p*(1-1/w_b) + p/w_b*(1+w_b*b*w_bp).^6;
Z_br = r*(1-1/w_b) + r/w_b .*(1+w_b*b*w_br).^6;
Z_brp = r*p*w_rp*(1-1/w_b) + r*p*w_rp/(w_b).*(1+w_b*b*w_bp*w_br*w_brp).^6;

% Calculate the P_bound
P_bound = (Z_bp + Z_brp)./ (Z_b + Z_bp + Z_br + Z_brp);

P_bound_null = (Z_bp)./(Z_bp + Z_b);

output = P_bound./P_bound_null;

%FC = P_bound./P_bound_null;
end