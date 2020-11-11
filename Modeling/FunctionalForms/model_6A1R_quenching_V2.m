function [output] = model_6A1R_quenching_V2(params, TF) 

%% Definition of parameters
% [Kb, Kr, w_b, w_bp, w_brp, p, R_max] = params;
Kb = params(1);
Kr = params(2);
w_b = params(3);
w_bp = params(4);
% w_ar = 1;
% w_rp = 1;
w_brp = params(5);
p = params(6);
R_max = params(7);

% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r = Runt./Kr;

% Calculate the partition function
Z_b = (1-1/w_b) + 1/w_b * (1+w_b*b).^6;
Z_bp = p*(1-1/w_b) + p/w_b*(1+w_b*b*w_bp).^6;
Z_br = r*(1-1/w_b) + r/w_b .*(1+w_b*b).^6;
Z_brp = r*p*(1-1/w_b) + r*p/(w_b).*(1+w_b*b*w_bp*w_brp).^6;

% Calculate the P_bound
P_bound = (Z_bp + Z_brp)./ (Z_b + Z_bp + Z_br + Z_brp);

%P_bound_null = (Z_ap)./(Z_a + Z_ap);
Rate = R_max*P_bound;

output = Rate;

%FC = P_bound./P_bound_null;
end