function [Rate] = model_6A1R_HillModel_V3_competition(params,TF) 

%% Definition of parameters
% [Kb, Kr, w_bp, w_rp, p, R_max] = params;
Kb = params(1);
Kr = params(2);
% w_b = params(3);
w_bp = params(3);
w_br = params(4);
p = params(5);
R_max = params(6);
% R_min = params(5);
% w_rp = 1;

% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r = Runt./Kr;

% Calculate the P_bound, and Rate
P_bound = p*(1+b.^6*w_bp + r + b.^6.*r*w_bp*w_br)./...
            (1 + p + b.^6 + r + b.^6.*r*w_br + b.^6*p*w_bp + r*p + b.^6.*r*p*w_bp*w_br);
        
Rate = R_max* P_bound;

%% For a dissection of each term's effect
end