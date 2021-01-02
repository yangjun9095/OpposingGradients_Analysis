function SS = SumofSquaresFunction_DirectRepression_V1(data, params)
% Log likelihood probability function for MCMC fitting.
% Currently the model infers values for . 
% This is for use in the MCMCstat package and returns the sum-of-squares function. 

% Parameters
% ----------
% construct: string.
% String specifying the construct used.
% data: structure
% data (MCMCdata) has fields like below. Note that they are already
% pre-processed as the AP range given in the main script.
% xdata : [APbins APbins]
% ydata : [Rate_null Rate_WT]
% Bcd : [Bcd Bcd]
% Run : [RunNull Run]

% params: array
% Vector of fit parameters (Kb, Kr, w_b, w_bp, p, R_max, w_br, w_rp, w_brp)

% Return
% ------
% SS: sum-of-squares residual

% Note that this assumes our sampling assumes Gaussian error with fixed
% variance/standard deviation for each measurement.

%% Extract fit parameters and data

% Create the interpolated AP bins for fitting
APbins = data.xdata;
% APbins = data.APbins;
% dx = median(diff(APbins))/10; % make 1/10 AP binning, thus 0.25% per bin
% APbins_interp = APbins(1):dx:APbins(end);

%Extract the Rate values
% We will use the pre-defined range of xdata and ydata for the fitting.
Bcd = data.Bcd;
Run = data.Run;

Rate = data.ydata;

Kb = params(1); % dissociation constant of Bcd
Kr = params(2); % dissociation constant of Run
w_b = params(3); % interaction between Bcd molecules
w_bp = params(4); % interaction between Bicoid and RNAP
w_rp = params(5); %
p = params(6); % 
R_max = params(7); %


% repression terms
% w_br = 1; %
% w_brp = 1; %

% Construct a parameter_set for the model input : this should be
% streamlined further to keep consistency in parameter order later.
% [Kb, Kr, w_a, w_ap, w_rp, p, R_max] = params;
params_set = [Kb Kr w_b w_bp w_rp p R_max];
%% Set up posterior function for MCMC
%Calculate the simulated Rate (initial slope)
Rate_sim = model_6A1R_direct_repression_V1(Bcd, Run, params_set);

%% Compute sum-of-squares
% Calculate the residuals of experimental and theoretical predictions
residuals = Rate- Rate_sim;

% Compute the sum-of-squares function
SS = nansum(residuals.^2);
end