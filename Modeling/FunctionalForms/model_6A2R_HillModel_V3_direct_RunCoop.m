function [Rate] = model_6A2R_HillModel_V3_direct_RunCoop(params,TF) 

%% Definition of parameters
% [Kb, w_bp, p, R_max, Kr1, Kr2, w_rp1, w_rp2] = params;
Kb = params(1);
w_bp = params(2);
p = params(3);
R_max = params(4);

% repressor parameters
Kr1 = params(5);
Kr2 = params(6);
w_rp1 = params(7);
w_rp2 = params(8);
% Run-Run cooperativity
w_rr = params(9);


% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r1 = Runt./Kr1;
r2 = Runt./Kr2;

% Calculate the P_bound, and Rate
Z = 1 + p + b.^6 + r1 + r2 +...
    b.^6.*(r1 + r2 + p*w_bp) + r1*p*w_rp1 + r2*p*w_rp2 + r1.*r2*w_rr+...
    r1.*r2*p*w_rp1*w_rp2*w_rr + b.^6*p*w_bp.*(r1*w_rp1 + r2*w_rp2) + b.^6.*r1.*r2*w_rr+...
    b.^6.*r1.*r2*p*w_bp*w_rp1*w_rp2*w_rr;

Z_p = p +  b.^6*p*w_bp + r1*p*w_rp1 + r2*p*w_rp2 +...
        r1.*r2*p*w_rp1*w_rp2*w_rr + b.^6*p*w_bp.*(r1*w_rp1 + r2*w_rp2)+...
         b.^6.*r1.*r2*p*w_bp*w_rp1*w_rp2*w_rr;

P_bound = Z_p./Z;
        
Rate = R_max* P_bound;

%% For a dissection of each term's effect
end