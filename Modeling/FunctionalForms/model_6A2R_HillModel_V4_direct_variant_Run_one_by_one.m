function [Rate] = model_6A2R_HillModel_V4_direct(params,TF) 

%% Hill.V4 assumes the saturation of Run sites, thus no need for K_r parameters.

%% Definition of parameters
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params;
Kb = params(1);
w_bp = params(2);
p = params(3);
R_max = params(4);

% repressor parameters
% Kr1 = params(5);
% Kr2 = params(6);
w_rp1 = params(5);
w_rp2 = params(6);


% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
% r1 = Runt./Kr1;
% r2 = Runt./Kr2;

% Calculate the P_bound, and Rate
Z = 3 + b.^6.*3 +...
    p*w_rp1 + p*w_rp2 +...
    p*w_rp1*w_rp2 + b.^6*p*w_bp.*(w_rp1 + w_rp2) +...
    b.^6.*p*w_bp*w_rp1*w_rp2;

Z_p = p*w_rp1 + p*w_rp2 +...
         p*w_rp1*w_rp2 + b.^6*p*w_bp.*(w_rp1 + w_rp2)+...
         b.^6.*p*w_bp*w_rp1*w_rp2;

P_bound = Z_p./Z;
        
Rate = R_max* P_bound;

%% For a dissection of each term's effect
end