function SS = SumofSquaresFunction_InitialSlopeMCMC(data,params)
% Log likelihood probability function for MCMC fitting.
% Currently the model infers values for . 
% This is for use in the MCMCstat package and returns the sum-of-squares function. 

% Parameters
% ----------
% construct: string.
%   String specifying the construct used.
% data: structure
%   Structure containing data 
% (the first column is AP position, the second column is averaged initial rate over embryos.
% params: array
%   Vector of fit parameters (Kb, Kr, w_b, w_bp, p, R_max, w_br, w_rp, w_brp)

% Return
% ------
% SS: sum-of-squares residual

% Note that this assumes our sampling assumes Gaussian error with fixed
% variance/standard deviation for each measurement.

%% Extract fit parameters and data

% %Create interpolated time for model fitting
% t = data.xdata; %Time vector
% dt = mean(t(2:end)-t(1:(end-1))); %Average time resolution of the dataset for model usage
% t_interp = t(1):dt:t(end); %Interpolated time vector with even time resolution

% Create the interpolated AP bins for fitting
APbins = data.APbins;
dx = median(diff(APbins))/10; % make 1/10 AP binning, thus 0.25% per bin
APbins_interp = APbins(1):dx:APbins(end);

%Extract the Rate values
Rate_null = data.Rate_null;
Rate_WT = data.Rate_WT;

Kb = params(1); % dissociation constant of Bcd
Kr = params(2); % dissociation constant of Run
w_b = params(3); % interaction between Bcd molecules
w_bp = x(4); % interaction between Bicoid and RNAP
p = x(5); % 
R_max = x(6); %
% repression terms
w_br = x(7); %
w_rp = x(8); %
w_brp = x(9); %

%% Set up posterior function for MCMC
%Calculate the simulated Rate (initial slope)
PolPos = ConstantElongationSim(v,ton,R_full,t_interp); %Simulated polymerase positions
Rate = model

%Interpolate the simulated fluorescence signals back to experimental time
%resolution
MS2 = interp1(t_interp,MS2,t);
PP7 = interp1(t_interp,PP7,t);
fluorSim = [MS2,PP7]; %Simulated fluorescences

%% Compute sum-of-squares
% Calculate the residuals of experimental and theoretical predictions
residuals = fluorExp - fluorSim;

% Compute the sum-of-squares function
SS = nansum(residuals.^2);
end