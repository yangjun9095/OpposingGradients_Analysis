function [HPD] = HighestPostDensity(trace, mass_frac)
%% Description
% Bayesian analysis - finding the Highest Posterior Density (HPD) from a
% given chain(trace) within mass_frac range (the default is 0.95, which is
% 95% HPD).

% Parameters
% trace : MCMC chain for a single variable.
% mass_frac : float with 0<mass_frac<=1
% The fraction of the prob. to be included in the HPD.
% For example, mass_frac = 0.95 gives a 95% HPD.

% Returns
% HPD

% Sort the trace (in ascending order)
trace_sorted = sort(trace);

% length of the trace (chain)
n_trace = length(trace_sorted);

% get number of samples that should be included in HPD
n_samples = floor(mass_frac*n_trace);

% get width of 
end