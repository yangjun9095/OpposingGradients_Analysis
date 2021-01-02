function [params_fit, Res, Jacobian, CI, STD, Ypred, delta] = ...
            global_fit_construct_generalThermoV2(inputTF, data, mdl0, mdl, params0, lb, ub)
%% model is given as input arguments. (both mdl0 and mdl)

% Set the parameter bounds and initial value for the query

% options = optimoptions('FunctionTolerance', 10^(-4),'Display','iter', 'Algorithm', 'trust-region-reflective');

optimoptions = optimset('TolFun',1E-8, 'MaxIter', 1E8, 'MaxFunEvals', 1E5);
                   
% fit using the lsqcurvefit
[params_fit,~,Res,~,~,~,Jacobian] =...
                    lsqcurvefit(mdl0, params0, inputTF, data, lb, ub, optimoptions);

%% Calculate the CI of the model fit and parameters
% First, get the CI for the fitted parameters
CI = nlparci(params_fit, Res, 'jacobian', Jacobian);

% Calculate the confidence interval of each parameter
STD = (CI(:,2) - CI(:,1)) /2;

% Second, calculate the CI for the predicted fit using nlpredci
[Ypred,delta] = nlpredci(mdl,inputTF, params_fit, Res, 'Jacobian', full(Jacobian));

end