function [Rate] = model_6A2R_HillModel_V3_quenching_fixed_Kr_Bcd_RNAPparams(params,TF, params_fixed) 

%% Definition of parameters
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params;
Kb = params_fixed(1);
w_bp = params_fixed(2);
p = params_fixed(3);
R_max = params_fixed(4);

% repressor parameters
K_r = params(1);
omega_rp1 = params(2);
omega_rp2 = params(3);


% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r1 = Runt./K_r;
r2 = r1;

% Calculate the P_bound, and Rate
Z = 1 + p + b.^6 + r1 + r2 +...
    b.^6.*(r1 + r2 + p*w_bp) + r1*p*omega_rp1 + r2*p*omega_rp2 + r1.*r2+...
    r1.*r2*p*omega_rp1*omega_rp2 + b.^6*p*w_bp.*(r1*omega_rp1 + r2*omega_rp2) + b.^6.*r1.*r2+...
    b.^6.*r1.*r2*p*w_bp*omega_rp1*omega_rp2;

Z_p = p +  b.^6*p*w_bp + r1*p*omega_rp1 + r2*p*omega_rp2 +...
        r1.*r2*p*omega_rp1*omega_rp2 + b.^6*p*w_bp.*(r1*omega_rp1 + r2*omega_rp2)+...
         b.^6.*r1.*r2*p*w_bp*omega_rp1*omega_rp2;

P_bound = Z_p./Z;
        
Rate = R_max* P_bound;

%% For a dissection of each term's effect
end