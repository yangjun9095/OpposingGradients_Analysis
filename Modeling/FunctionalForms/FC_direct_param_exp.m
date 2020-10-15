function FC = FC_direct_param_exp(params)

b = params(1);
r = params(2);
p = params(3);
w_b = params(4);
w_bp = params(5);
w_rp = params(6);

% mini-Partition functions
Z_b = 1 + 1/w_b*((1+b*w_b)^6-1);
Z_br = r+r/w_b*((1+b*w_b)^6-1);
Z_bp = p + p/w_b*((1+b*w_b*w_bp)^6 -1);
Z_brp = r*p*w_rp + r*p*w_rp/w_b*((1+b*w_b*w_bp)^6 - 1);

P_bound = (Z_bp + Z_brp) / (Z_b + Z_br + Z_bp + Z_brp);

% P_bound where there's no Runt protein
P_bound_null = Z_bp / (Z_b + Z_bp); 

FC = P_bound/P_bound_null;
end