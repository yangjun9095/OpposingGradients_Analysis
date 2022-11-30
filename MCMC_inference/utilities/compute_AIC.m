function [AIC] = compute_AIC(model, data, data_sigma,num_params)
%% Description
% AIC: This function computes the Akaike Information Criterion from the result
% of the Markov Chain Monte Carlo parameter estimation.
% AIC = -2*log(Likelihood) + 2*(number of parameters)

% We computed the log-likelihood by using sum of squares function, thus we
% will just modify the result of that function slightly to compute the AIC.

%% input
% model: expected values from the model with a set of parameters (from MCMC
% parameter estimation).
% data: experimental data (observed values). Note that the model and data
% should have the same dimension, 1xN vector.
% data_sigma: standard deviation from each data point.
% num_params = number of parameters

%% compute the log-likelihood
Log_likelihood = -1/2*(nansum((model-data).^2./(data_sigma.^2)) + nansum(log(2*pi*data_sigma.^2)))

%% compute the AIC
AIC = -2*Log_likelihood + 2*num_params;
end