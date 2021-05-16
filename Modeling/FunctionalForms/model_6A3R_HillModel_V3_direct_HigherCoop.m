function [Rate] = model_6A3R_HillModel_V3_direct_HigherCoop(params,TF) 

%% Definition of parameters
% assumptions : 
% 1) strong Bicoid-Bicoid cooperativity (w_b >>1), thus all-or-none.
% 2) K_{r} for the Runt binding sites are all equivalent.

% [Kb, w_bp, p, R_max,...
%   Kr, w_rp1, w_rp2, w_rp3,...
%   w_rr1, w_rr2, w_rr3,...
%   w_ho1, w_ho2, w_ho3, w_hoho] = params;
Kb = params(1);
w_bp = params(2);
p = params(3);
R_max = params(4);

% repressor parameters
Kr = params(5);
w_rp1 = params(6);
w_rp2 = params(7);
w_rp3 = params(8);
% Run-Run cooperativity
w_rr1 = params(9);
w_rr2 = params(10);
w_rr3 = params(11);

% higher-order cooperativity
w_ho1 = params(12);
w_ho2 = params(13);
w_ho3 = params(14);

% another higher-order term
w_hoho = params(15);

% TF inputs : Read a matrix of TF, each column represent different TFs
Bcd = TF(:,1);
Runt = TF(:,2);

% scale with the dissociation constant
b = Bcd./Kb;
r = Runt./Kr;

% Calculate the P_bound, and Rate
Z = 1 + p + b.^6 + 3*r + ...
    b.^6.*p*w_bp + 3*r.*b.^6 + r*p*(w_rp1 + w_rp2 +w_rp3) + r.^2*(w_rr1+w_rr2+w_rr3) + ... % 5C2
    b.^6.*r*p*w_bp*(w_rp1 + w_rp2 +w_rp3) + b.^6.*r.^2*(w_rr1+w_rr2+w_rr3) +...
    r.^2*p*(w_rp1*w_rp2*w_rr3 + w_rp2*w_rp3*w_rr1 + w_rp3*w_rp1*w_rr2) + r.^3*w_rr1*w_rr2*w_rr3 + ...  % 5C3
    r.^3*p*(w_rp1*w_rp2*w_rp3*w_rr1*w_rr2*w_rr3*w_ho1*w_ho2*w_ho3*w_hoho) +...
    b.^6.*r.^3*w_rr1*w_rr2*w_rr3 + ...
    r.^2.*b.^6*p*w_bp*(w_rp1*w_rp2*w_rr3 + w_rp2*w_rp3*w_rr1 + w_rp3*w_rp1*w_rr2) +...
    b.^6.*r.^3.*p*(w_rp1*w_rp2*w_rp3*w_rr1*w_rr2*w_rr3*w_ho1*w_ho2*w_ho3*w_hoho);

Z_p = p +  b.^6*p*w_bp + r*p*(w_rp1 + w_rp2 +w_rp3) + b.^6.*r*p*w_bp*(w_rp1 + w_rp2 +w_rp3) +...
        r.^2*p*(w_rp1*w_rp2*w_rr3 + w_rp2*w_rp3*w_rr1 + w_rp3*w_rp1*w_rr2) + ...
        r.^3*p*(w_rp1*w_rp2*w_rp3*w_rr1*w_rr2*w_rr3*w_ho1*w_ho2*w_ho3*w_hoho) + ...
        r.^2.*b.^6*p*w_bp*(w_rp1*w_rp2*w_rr3 + w_rp2*w_rp3*w_rr1 + w_rp3*w_rp1*w_rr2) +...
        b.^6.*r.^3.*p*(w_rp1*w_rp2*w_rp3*w_rr1*w_rr2*w_rr3*w_ho1*w_ho2*w_ho3*w_hoho);
        
        

P_bound = Z_p./Z;
        
Rate = R_max* P_bound;

%% For a dissection of each term's effect
end