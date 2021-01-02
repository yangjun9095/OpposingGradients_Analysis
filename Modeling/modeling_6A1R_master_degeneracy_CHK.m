function modeling_6A1R_master_function(varargin)
% Description : 
% 
% Last updated : 11/25/2020, by YJK

% Assumptions : 

%% Step0. Set up the directories
% file to save
SavePath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo';
%% Step1. Get the actual input TF and output data
% Here, I'll use the time-averaged Bcd and Runt profiles processed by
% main01_04/05 scripts.
% Actually, the time-averaging "time window" doesn't matter that much,
% since the gradient scales nicely over time. But, we use the 0-10 min into
% nc14.

FilePath = 'S:\YangJoon\Dropbox\OpposingGradient\OpposingGradients_ProcessedData\AveragedDatasets_Feb2020';

Data = load([FilePath, filesep, 'compiledData.mat']);

compiledData = Data.compiledData;

DataTypes = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};

%% If the data pre-processing steps are done, then just import the output
TFData = load([FilePath, filesep, 'TFinput.mat']);
TFdata = TFData.TFinput;

Bicoid = TFdata(:,1);
Runt = TFdata(:,2);
RuntNull = TFdata(:,3);

% Make a matrix whose each column is each TF
TF(:,1) = Bicoid;
TF(:,2) = Runt;
TF(:,3) = RuntNull;

%% Model-type1 : direct repression

%% Pre-process the data (both input TF and output rate(slope)) such that we can fit to a given model.

% Import the real data from the compiledData
% Then, pick one dataset as our starting point.(let's start with one Run site, for example, [001])                
% Define the AP range
APaxis = 0:0.025:1;

% Pick the construct for fitting
% N-th element in the DataTypes nomenclature structure.
data_model_input = struct;
k=1; %counter
for construct = [2,5,6]
    % initialize the variables
    data = [];
    input_combined =[];
    data_interp = [];
    input_combined_interp = [];

    Rate = compiledData{construct+1,9};
    Rate_SEM = compiledData{construct+1,10};

    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};

    % Define the range of fitting
    APaxis = 0:0.025:1;
    
    if construct == 2
        APpos1 = 20;% [% of embryo length]
        APpos2 = 42.5;% [% of embryo length]
    else
        APpos1 = 20;% [% of embryo length]
        APpos2 = 45;% [% of embryo length]
    end

    APbin1 = APpos1/2.5 + 1;
    APbin2 = APpos2/2.5 + 1;

    fitRange = APbin1:APbin2;

    % Trim the data of WT and Null into one vector
    data = [Rate_null(fitRange); Rate(fitRange)];

    input_combined(:,1) = [Bicoid(fitRange); Bicoid(fitRange)];
    input_combined(:,2) = [RuntNull(fitRange); Runt(fitRange)];
    

    % save into a structure
    data_model_input(k).construtName = DataTypes(construct);
    data_model_input(k).APbins = APaxis(fitRange);
    data_model_input(k).input = input_combined;
    data_model_input(k).output = data;
    k=k+1;

end   

    
%% Fitting process using the global_fit_construct_generalThermoV2.m
% Model : Direct repression
% Inputs : xdata, ydata, and model for the fitting

% Pick a model from different modes
% Define anonymous function from the given model input. (either direct,
% competition, quenching, etc.)
mdl0 = @model_6A1R_direct_repression_V2;
mdl = @(params, TF) model_6A1R_direct_repression_V2(params, TF);

% Define the parameter bounds for the fitting
% params = [Kb, Kr, w_a, w_ap, w_rp, p, R_max];
lb = [10^(-1) 10^(-1) 1 1 10^(-3) 0 50];
ub = [10^4 10^3 10^3 10^3 1 1 1000];

% Initialize the FitResult structure
FitTemp = struct;

construct = 5; % [001] construct (r1-close)
k=2; % indexing from the data_model_input structure.
% define the inputs for the fitting
inputTF = data_model_input(k).input;
data = data_model_input(k).output;

n_sim = 10;

for i = 1:n_sim
    
    % define the initial parameter guess : using rand
    % as it seems that the lsq fitting is dependent on the initial value
    % for some reason...
%     params0 = [100*rand, 100*rand, 10^2*rand, 10^2*rand, rand, rand, 350];
    if i==1
        params0 = [10, 5, 2, 2, 0.2, 0.001, 300]; % example
    else
        params0 = [10*rand, 5*rand, 100*rand, 1000*rand, 0.2, 0.001, 300];
    end
    
    % initialize the variables
    params_fit = [];
    Res = [];
    Jacobian = [];
    CI = [];
    STD = [];
    Ypred = [];
    delta = [];
    

    
    [params_fit, Res, Jacobian, CI, STD, Ypred, delta] = ...
                global_fit_construct_generalThermoV2(inputTF, data, mdl0, mdl, params0, lb, ub);
            
    % save the fit results (parameters and Jacobian, etc.) into the FitResult structure.        
    FitTemp(i).constructName = DataTypes(construct);
    FitTemp(i).params0 = params0;
    FitTemp(i).params_fit = params_fit;
    FitTemp(i).Res = Res;
    FitTemp(i).Jacobian = Jacobian;
    FitTemp(i).CI = CI;
    FitTemp(i).STD = STD;
    FitTemp(i).Ypred = Ypred;
    FitTemp(i).delta = delta;
   
end

%% pick a couple of fits and params for demonstration of degeneracy
% index = [1, 4, 5, 12, 13, 14, 15, 21, 51];
k=1; % counter
params_deg = [];
fit = [];
fit_SE = [];

for i = 1:n_sim
    params_deg(i,:) = FitTemp(i).params_fit;
    fit(i,:) = FitTemp(i).Ypred;
    fit_SE(i,:) = FitTemp(i).delta;
%     k=k+1;
end

%% generate plots of fits and params
% first, fits

construct = 5; %[001]
% initialize the variables
Rate = [];
Rate_SEM = [];
Rate_null = [];
Rate_null_SEM = [];

% Pull the relevant data from the compiledData
Rate = compiledData{construct+1,9};
Rate_SEM = compiledData{construct+1,10};

Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

% pull the fitted result and error bar (Ypred and delta) from the
% lsqcurvefit fitting.
APrange = 9:19;

hold on
% Runt nulls
errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
% Runt WT
errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))

xlim([0.2 0.6])
xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 400])
yticks([0 100 200 300 400])

% legend([h(1) h(2)],'Runt null','Runt WT')
xlabel('embryo length')
ylabel({'initial RNAP', 'loading rate (AU/min)'})
box on
StandardFigure(gcf,gca)

% plot different fits from different sets of parameters
colorindex = [2, 3, 5, 8];
k=1;  % counter

% Pick the fits that looked good. Other ones seem to be stuck in the local
% minima.
for i = [1,6,7,8]%1:n_sim%length(index)
    plot(APaxis(APrange), fit(i,1:length(APrange)),'Color',ColorChoice(colorindex(k),:) ,'LineWidth',2)
    plot(APaxis(APrange), fit(i,1+length(APrange):end), 'Color',ColorChoice(colorindex(k),:) ,'LineWidth',2)
    k=k+1;
%     pause
end
% shadedErrorBar(APrange, Ypred(1:length(APrange)), delta(1:length(APrange)),'lineProps',{'markerfacecolor',ColorChoice(4,:)})
% shadedErrorBar(APrange, Ypred(length(APrange)+1:end), delta(length(APrange)+1:end),'lineProps',{'markerfacecolor',ColorChoice(1,:)})


%   %Save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo\Direct_Degeneracy_test';
saveas(gcf,[FigPath,filesep,'raw_fits_degenerate_',constructNames{construct},'.tif']); 
saveas(gcf,[FigPath,filesep,'raw_fits_degenerate_',constructNames{construct},'.pdf']); 



%% Plot the degenerate parameters
% names = {'Kb','Kr','w_b','w_{bp}','w_{rp}','p','R_{max}'};

% plot different sets of parameters (degenerate ones)
colorindex = [2, 3, 5, 8];
k=1;  % counter

figure
hold on
for i=[1,6,7,8]
    plot(1:7, params_deg(i,:),'o',...
        'Color', ColorChoice(colorindex(k),:),...
        'LineWidth', 2)
    % 'MarkerFaceColor',ColorChoice(colorindex(k),:),...
    k=k+1;
end

legend('set1','set2','set3','set4','Location','SouthWest')

set(gca, 'YScale', 'log')

xlim([0 8])
xticklabels({'','K_{b}','K_{r}','\omega_{b}','\omega_{bp}','\omega_{rp}','p','R_{max}',''})
xlabel('parameters')
ylabel('fitted values')

StandardFigure(gcf,gca)

% Save the plot
FigPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput\ModelingV2_generalizedThermo\Direct_Degeneracy_test';
saveas(gcf,[FigPath,filesep,'fittedParams_degenerate_',constructNames{construct},'.tif']); 
saveas(gcf,[FigPath,filesep,'fittedParams_degenerate_',constructNames{construct},'.pdf']); 











%% Extract the fitted parameters (this is, in some sense, similar to the Markov chain).
% we will plot in a way that we make cornerplots in MCMC.
% Also, this will show us how degenerate the model is.

for i=1:n_sim
    params_fitted(i,:) = FitTemp(i).params_fit;
end

















%% generate a corner plot
m = [params_fitted(:,1), params_fitted(:,2), params_fitted(:,3),...
        params_fitted(:,4), params_fitted(:,5), params_fitted(:,6), params_fitted(:,7)];
corner = figure;
% names = {'Kb','Kr','w_b','w_{bp}','w_{rp}','p','R_{max}'};
ecornerplot(m,'names',names);


%% histogram

for i=1:length(params_fit)
    clf
    histogram(params_fitted(:,i))
    xlabel(names{i})
    pause
end
%% check the fitted parameters by pluggin back into the model,
% to see whether the fit matches well with the data

%% plot the fitting results
APaxis = 0:0.025:1;

% counter
k=1;
for construct = [2,5,6]
    clf
    
    % initialize the variables
    Rate = [];
    Rate_SEM = [];
    Rate_null = [];
    Rate_null_SEM = [];
    
    % Pull the relevant data from the compiledData
    Rate = compiledData{construct+1,9};
    Rate_SEM = compiledData{construct+1,10};

    Rate_null = compiledData{construct+1+8,9};
    Rate_null_SEM = compiledData{construct+1+8,10};
    
    % pull the fitted result and error bar (Ypred and delta) from the
    % lsqcurvefit fitting.
    APrange = data_model_input(k).APbins;
    Ypred = FitTemp(k).Ypred;
    delta = FitTemp(k).delta;
    k=k+1;
    
    hold on
    % Runt nulls
    errorbar(APaxis, Rate_null, Rate_null_SEM,'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))
    shadedErrorBar(APrange, Ypred(1:length(APrange)), delta(1:length(APrange)),'lineProps',{'markerfacecolor',ColorChoice(4,:)})
    % Runt WT
    errorbar(APaxis, Rate, Rate_SEM,'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
    shadedErrorBar(APrange, Ypred(length(APrange)+1:end), delta(length(APrange)+1:end),'lineProps',{'markerfacecolor',ColorChoice(1,:)})


    xlim([0.2 0.6])
    xticks([0.2 0.3 0.4 0.5 0.6])
    ylim([0 400])
    yticks([0 100 200 300 400])

    xlabel('embryo length')
    ylabel({'initial RNAP', 'loading rate (AU/min)'})

    StandardFigure(gcf,gca)
    pause
%   %Save the plot
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,'fit_direct_V3_',constructNames{construct},'.pdf']); 
end

%% Save the fitting results into structures
save(['S:\YangJoon\Dropbox\OpposingGradientsFigures\PipelineOutput', filesep,...
    'ModelingV2_generalizedThermo\DirectRepression\001_degeneracy_Check.mat'],...
        'data_model_input', 'FitTemp')

end