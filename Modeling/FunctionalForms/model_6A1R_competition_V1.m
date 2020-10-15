function [output] = model_6A1R_direct_repression_V1(Bcd, Runt,...
                        params) 

%% Definition of parameters
% [Kb, Kr, w_a, w_ap, w_ar, p] = params;
Kb = params(1);
Kr = params(2);
w_a = params(3);
w_ap = params(4);
w_ar = params(5);
% w_rp = 1;
p = params(6);
R_max = params(7);

a = Bcd./Kb;
r = Runt./Kr;

% Calculate the partition function
Z_a = (1-1/w_a) + 1/w_a * (1+w_a*a).^6;
Z_ap = p*(1-1/w_a) + p/w_a*(1+w_a*a*w_ap).^6;
Z_ar = r*(1-1/w_a) + r./w_a .*(1+w_ar*w_a*a).^6;
Z_arp = r*p*(1-1/w_a) + r*p/w_a.*(1+w_ar*w_a*a*w_ap).^6;

% Calculate the P_bound
P_bound = (Z_ap + Z_arp)./ (Z_a + Z_ap + Z_ar + Z_arp);

% P_bound_null = (Z_ap)./(Z_a + Z_ap);
% 
% FC = P_bound./P_bound_null;
Rate = R_max * P_bound;
output = Rate;
end