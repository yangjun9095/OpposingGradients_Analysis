function [Rate] = model_6A0R_HillModel_V2_p_constrained(params,TF) 

%% Definition of parameters
% [Kb, w_bp, p, R_max] = params;
Kb = params(1);
w_bp = params(2);
p = params(3);
R_max = params(4);
R_min = params(5);

% repressor parameters
% Kr1 = params(5);
% Kr2 = params(6);
% w_rp1 = params(7);
% w_rp2 = params(8);


% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
% Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
% r1 = Runt./Kr1;
% r2 = Runt./Kr2;

% constraint on the "p" using the Bcd->0 limit, as in Liz&Jonathan's paper.
p = R_min/(R_max - R_min); 

% Calculate the P_bound, and Rate
Z = 1 + p + b.^6 + b.^6.*p*w_bp;

Z_p = p + b.^6*p*w_bp;

P_bound = Z_p./Z;
        
Rate = R_max* P_bound;

%% For a dissection of each term's effect
end