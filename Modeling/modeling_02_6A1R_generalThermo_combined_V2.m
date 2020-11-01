%% Thermodynamic modeling script V2.
% Difference from V1 models: 
% 1) generalized thermodynamic model, thinking of each Bcd binging case
% independently.
% 2) Use lsqcurvefit instead of lsqnonlin to use nlpredci (calculating the
% CI of the output).
% 3) For 3 datasets ([001],[010],[100]), for different types of models, we
% perform non-linear regression fitting (chi-squared), then extract fitted
% parameters with their CI, as well as the CI of the output.

% 4) input TF data should be read into a matrix with different columns for
% different TFs (Bcd, Runt, and Runt nulls)





%% 
options = optimoptions('lsqnonlin','Display','iter', 'Algorithm', 'trust-region-reflective');

[x,resnorm,residual,exitflag,output,lambda,jacobian] =...
                            lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options) 
                        
                        
                        
                        
                        
                        
function yfit = subfun_6A1R_combination(params,X)

Bcd = X(:,1); % Bcd
Run = X(:,2); % Runt
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
omegaDP = params(nSets + 3);
offset = params(nSets + 4);

yfit = amplitude.*((n+(x./KD(dsid)).*n.*omegaDP)./(1+x./KD(dsid)+n+(x./KD(dsid)).*n.*omegaDP))+offset;

end


function yfit = subfun_hill(params,X)

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
offset = params(nSets + 3);

yfit = amplitude.*(((x./KD(dsid)).^n)./(1+((x)./KD(dsid)).^n))+offset;

end


function yfit = subfun_simplebinding_weak(params,X)

%simplebinding in the weak promoter limit.

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
omegaDP = params(nSets + 2);
offset = params(nSets + 3);

yfit = amplitude.*(((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP))+offset;

end                        