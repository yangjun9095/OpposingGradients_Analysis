function [Output] = model_6A1R_direct_repression_V3_R_p_constrained(params,TF) 

%% Definition of parameters
% [Kb, Kr, w_a, w_ap, w_rp, R_max, R_min] = params;
% Note that the p and R_max are now constrained with other parameters and
% data, especially with R_max and R_min
Kb = params(1);
Kr = params(2);
w_b = params(3);
w_bp = params(4);
w_rp = params(5);

% dependent parameters
R_max = params(6);
R_min = params(7);


% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r = Runt./Kr;

R = max(r); % repressor value at the b->0 limit.

% calculate the p based on other parameters (to constrain it)
p = (1+R)/(1+R*w_rp).*R_min./(R_max - R_min);

% Calculate the partition function
Z_b = (1-1/w_b) + 1/w_b * (1+w_b*b).^6;
Z_bp = p.*(1-1/w_b) + p./w_b.*(1+w_b*b*w_bp).^6;
Z_br = r*(1-1/w_b) + r./w_b .*(1+w_b*b).^6;
Z_brp = r.*p.*(1-w_rp/w_b) + w_rp*r.*p./w_b.*(1+w_b*b*w_bp).^6;

% Calculate the P_bound
P_bound = (Z_bp + Z_brp)./ (Z_b + Z_bp + Z_br + Z_brp);

Output = P_bound*R_max;
end