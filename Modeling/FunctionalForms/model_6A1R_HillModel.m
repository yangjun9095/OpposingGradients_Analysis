function [Rate] = model_6A1R_HillModel(params,TF) 

%% Definition of parameters
% [Kb, Kr, w, R_max, R_min] = params;
Kb = params(1);
Kr = params(2);
% w_b = params(3);
% w_bp = params(4);
% w_rp = params(5);
% p = params(6);
w = params(3);
R_max = params(4);
R_min = params(5);

% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r = Runt./Kr;

% Calculate the Rate
Rate = (R_min + R_max*b.^6.*(1 + r.*w))./(1 + b.^6+ r + b.^6.*r.*w );

end