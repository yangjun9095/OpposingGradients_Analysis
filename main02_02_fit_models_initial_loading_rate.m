function main02_02_fit_models_initial_loading_rate
%% Part2. Fitting the initial loading rates (mean) with theoretical models
% Here, I will use Mean12, Mean13, and Mean14 for fitting with some
% expression for different constructs with different numbers of binding
% sites. Note that the Mean12,13, and 14 means that initial rate of RNAP
% loading fit from individual embryos, then averaged.

% Load the Bcd, Runt data
Bicoid = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed/Bcd-Averaged');
Runt = load('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\2018-05-24-Runt-JB3-MCP-mCherry-vasa-eGFP1\CompiledNuclei');

% First, we should build a prediction matrix using handful of parameters.
% Use Rate_r0.m for optimizing the parameters, r_max, r_basal, and K_d
% using lsqnonlin.

% Second, we need to decide which nc we're going to model.
% One good starting place would be nc13, although we still don't have a
% good predictive power on the whole timecourse. Or, the beginning of nc14
% might be good as well.

% Now, I will start from 
BcdFluo = Bicoid.MeanVectorAP(Bicoid.nc14+20,:);
RuntFluo = Runt.MeanVectorAP(Runt.nc14+20,:);

% plot to check
hold on
plot(BcdFluo)
plot(RuntFluo)
title('Bcd and Runt profile over AP')
xlabel('AP')
ylabel('Protein concentration (AU)')
legend('Bcd','Runt')

%% limit the AP bins for fitting
APbinstart = 11;
APbinend = 26;

% Initial condition
R_max=500; % maximum
R_bas=100; % basal, probably due to ectopic expression...
Kd = 0.1;
% fitting for r0 first,
p0=[0.1,R_max,R_bas];
%lb=[0,R_max,R_bas];
lb = [0, 0, 0];
ub=[1,R_max,R_bas];
options = optimset('Display','iter');

f= @(p)Rate_r0(p,BcdFluo(APbinstart:APbinend))-Mean14(1,APbinstart:APbinend);

P = lsqnonlin(f, p0, lb, ub, options);



end