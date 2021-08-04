function [Rate] = model_6A2R_HillModel_V3_compete_fixed_Kr_Bcd_RNAP_w_rp_params(params,TF, params_fixed) 

%% model function for 6A2R_HillV3_direct with fixed Bcd/RNAP/Kr,w_rp params
% (cooperativities considered, for Runt-Runt cooperativity, w_rr, and/or
% higher-order cooperativity, w_ho)

%% Definition of parameters

% params_fixed = [Kb, w_bp, p, R_max, K_r, omega_rp1, omega_rp2]
Kb = params_fixed(1);
w_bp = params_fixed(2);
p = params_fixed(3);
R_max = params_fixed(4);

% repressor parameters
K_r = params_fixed(5);
omega_br1 = params_fixed(6);
omega_br2 = params_fixed(7);

% params
omega_rr = params(1);
omega_ho = params(2);

% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r1 = Runt./K_r;
r2 = r1;

% Calculate the P_bound, and Rate
Z = 1 + p + b.^6 + r1 + r2 +...
    b.^6.*(r1*omega_br1 + r2*omega_br2 + p*w_bp) + r1*p + r2*p + r1.*r2.*omega_rr +...
    r1.*r2*p.*omega_rr + ...
    b.^6*p*w_bp.*(r1*omega_br1 + r2*omega_br2) + b.^6.*r1.*r2.*omega_rr.*omega_br1.*omega_br2.*omega_ho +...
    b.^6.*r1.*r2*p*w_bp*omega_br1*omega_br2.*omega_rr.*omega_ho;

Z_p = p +  b.^6*p*w_bp + r1*p + r2*p +...
        r1.*r2*p.*omega_rr +...
        b.^6*p*w_bp.*(r1*omega_br1 + r2*omega_br2)+...
         b.^6.*r1.*r2*p*w_bp*omega_br1*omega_br2.*omega_rr.*omega_ho;

P_bound = Z_p./Z;
        
Rate = R_max* P_bound;

%% For a dissection of each term's effect
end