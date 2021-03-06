% generate plots of initial slopes, fold-change, etc.
% Name : make_InitialSlope_fig

clear 
close all
addpath('../utilities')
% set ID variables
DropboxFolder = 'S:/YangJoon/Dropbox/OpposingGradient';
FigureRoot = 'S:/YangJoon/Dropbox/OpposingGradientsFigures/PipelineOutput';

% load data structure
load([DropboxFolder,filesep,'OpposingGradients_ProcessedData',filesep,'AveragedData.mat'])

FigPath = [FigureRoot, filesep, 'MS2Traces'];
mkdir(FigPath)

%% Data info
DataTypesForFit = {'r0-new','r1-new','r2-new','r3-new','r1-close','r1-mid','r2_1+2','r2_1+3',...
                    'r0_RuntNull','r1_RuntNull','r2_RuntNull','r3_RuntNull',...
                    'r1_close_RuntNull','r1_mid_RuntNull','r2_close_RuntNull','r2_far_RuntNull'};

% Name of constructs so that we can label plots and files
constructNames = {'000','100','011','111','001','010','110','101',...
                    '000, null','100, null','011, null','111, null','001, null','010, null','110, null','101, null'};
%% Color module
% This is defining the line color
% We have 8 distinct datasets, with or without Runt protein.
% I think selecting 8 distinguishable color sets, then changing the
% brightness by either adding/subtracting white would be a better idea than
% selecting 16 different color sets.

colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.purple = [171,133,172]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

%colorDict.magenta = [208,109,171]/255;
%colorDict.lightBlue = [115,142,193]/255;
colorDict.lightgreen = [205,214,209]/255;
colorDict.pink = [232,177,157]/255;
colorDict.thickpink = [132,27,69]/255;

% Define a color matrix, 8 colors right now.
ColorChoice = [colorDict.blue; colorDict.green;...
                colorDict.yellow; colorDict.red; colorDict.brown;...
                colorDict.purple; colorDict.darkgreen; colorDict.thickpink]; 

% For now, I'll add white (color+[1 1 1])/2 to make thinner color (for the
% Runt nulls)
%% Plot module
% define the figure handle, fig_name
fig_name = figure; 
hold on
errorbar(X, Y, Z, 'LineWidth',2,'Color',ColorChoice(1,:))

% xTicks, yTicks
xticks([0.2 0.3 0.4 0.5 0.6 0.7 0.8])

set(gca,'yticklabel',[])

% no title, no-caps on the axis labels
xlabel('')
ylabel('')

legend('','','Location','SouthWest')

box on

StandardFigure(fig_name, fig_name.CurrentAxes)

% Save the plot
figPath = '';
saveas(gcf,[figPath,filesep,'name','.tif']); 
saveas(gcf,[figPath,filesep,'name','.pdf']); 

%% MS2 time traces (averaged over multiple embryos) for each construct with/without Runt protein
APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_MS2 = figure;

% Choose the APbin that we want to look at
APbin = 13; % 30% of the embryo length
for construct=1:8 % total number of constructs (enhancers)
    clf
    vars = {'Data','Data_null','Time','Time_null','Fluo','Fluo_null',...
                'Fluo_SEM','Fluo_SEM_null','NC14','NC14_null'};
    clear (vars{:})
    
    % Extract the NC14 traces out of the master structure
    Data = AveragedData{construct+1,2};
    Data_null = AveragedData{construct+1+8,2};
    
    Time = Data.ElapsedTime;
    Fluo = Data.MeanVectorAP;
    Fluo_SEM = Data.SEVectorAP;
    NC14 = Data.nc14;
    
    Time_null = Data_null.ElapsedTime;
    Fluo_null = Data_null.MeanVectorAP;
    Fluo_SEM_null = Data_null.SEVectorAP;
    NC14_null = Data_null.nc14;
    

    hold on
    errorbar(Time(NC14:end) - Time(NC14), Fluo(NC14:end,APbin) ,...
                Fluo_SEM(NC14:end,APbin),'LineWidth',2,'Color',ColorChoice(construct,:))
    errorbar(Time_null(NC14_null:end) - Time_null(NC14_null), Fluo_null(NC14_null:end,APbin) ,...
                Fluo_SEM_null(NC14_null:end,APbin),'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2 )

    % xTicks, yTicks
    %xlim([0 40])
    ylim([0 1200])
    xticks([0 10 20 30 40])
    %yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('time into nc14 (min)')
    ylabel('mean fluorescence (AU)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_MS2, fig_MS2.CurrentAxes)
%    pause(1)

    % Save the plot
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
%     saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% Check the synchronization for individual DataTypes
% % [101]
% Data = AveragedData{17,2};
% 
% Time = Data.ElapsedTime;
% Fluo = Data.MeanVectorAP;
% Fluo_SEM = Data.SEVectorAP;
% Fluo_ind = Data.MeanVectorAP_individual;
% nc14 = Data.nc14;
% 
% hold on
% errorbar(Time(nc14:end) - Time(nc14), Fluo(nc14:end, APbin), Fluo_SEM(nc14:end, APbin))
% 
% for i=1:length(Fluo_ind(1,1,:))
%     plot(Time(nc14:end), Fluo_ind(nc14:end, APbin, i))
% end

%% averaged MS2 traces over ALL nuclei(averaged over multiple embryos) for each construct with/without Runt protein

% Make a new directory
FigPath = [FigureRoot, filesep, 'MS2Traces_ALLNuclei'];
mkdir(FigPath)

APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_MS2 = figure;

% Choose the APbin that we want to look at
APbin = 13; % 30% of the embryo length
for construct=1:8 % total number of constructs (enhancers)
    clf
    vars = {'Data','Data_null','Time','Time_null','Fluo','Fluo_null',...
                'Fluo_SEM','Fluo_SEM_null','NC14','NC14_null'};
    clear (vars{:})
    
    % Extract the NC14 traces out of the master structure
    Data = AveragedData{construct+1,2};
    Data_null = AveragedData{construct+1+8,2};
    
    Time = Data.ElapsedTime;
    Fluo = Data.MeanVectorAP;
    Fluo_SEM = Data.SEVectorAP;
    NC14 = Data.nc14;
    Fraction_inst = Data.FractionON;

    Time_null = Data_null.ElapsedTime;
    Fluo_null = Data_null.MeanVectorAP;
    Fluo_SEM_null = Data_null.SEVectorAP;
    NC14_null = Data_null.nc14;
    Fraction_inst_null = Data_null.FractionON;
    
    % for [000], [111] Runt nulls, we don't have Histone marker, thus
    % having the FractionON info as # of spots/APbin area. We need to
    % scale this to [0-1] scale as other datasets.
    if construct==1 || construct==4 % [000]/[111], Runt nulls 
        [row,col] = size(Fraction_inst_null);
        scale = prctile(reshape(Fraction_inst_null,[1,row*col]),99);
        Fraction_inst_null_rescaled = Fraction_inst_null./scale;
        Fraction_inst_null = Fraction_inst_null_rescaled;
%         Fraction_inst_SEM_null_rescaled = Fraction_inst_SEM_null./scale;
    else
    end

    hold on
    % Runt WT
    errorbar(Time(NC14:end) - Time(NC14),...
                Fluo(NC14:end,APbin).*Fraction_inst(NC14:end,APbin) ,...
                Fluo_SEM(NC14:end,APbin).*Fraction_inst(NC14:end,APbin),...
                'LineWidth',2,'Color',ColorChoice(construct,:))
    % Runt null        
    errorbar(Time_null(NC14_null:end) - Time_null(NC14_null),...
                Fluo_null(NC14_null:end,APbin).*Fraction_inst_null(NC14_null:end,APbin) ,...
                Fluo_SEM_null(NC14_null:end,APbin).*Fraction_inst_null(NC14_null:end,APbin),...
                'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2 )

    % xTicks, yTicks
    %xlim([0 40])
    ylim([0 1200])
    xticks([0 10 20 30 40])
    %yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('time into nc14 (min)')
    ylabel('mean fluorescence (AU)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_MS2, fig_MS2.CurrentAxes)
%    pause(1)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% Make fraction of active nuclei (instantaneous) plot
% Make a new directory
FigPath = [FigureRoot, filesep, 'InstantFractionON'];
mkdir(FigPath)

APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_MS2 = figure;

% Choose the APbin that we want to look at
APbin = 13; % 30% of the embryo length
for construct=1:8 % total number of constructs (enhancers)
    clf
    vars = {'Data','Data_null','Time','Time_null','Fluo','Fluo_null',...
                'Fluo_SEM','Fluo_SEM_null','NC14','NC14_null',...
                'frac_individual_rescaled','Fraction_inst','Fraction_inst_null'};
    clear (vars{:})
    
    % Extract the NC14 traces out of the master structure
    Data = AveragedData{construct+1,2};
    Data_null = AveragedData{construct+1+8,2};
    
    % WT
    [~,~,nEmbryos] = size(Data.MeanVectorAP_individual);
    NC14 = Data.nc14;
    Time = Data.ElapsedTime(NC14:end) - Data.ElapsedTime(NC14);
%     Fluo = Data.MeanVectorAP(NC14:end,:);
%     Fluo_SEM = Data.SEVectorAP(NC14:end,:);
    
    Fluo_ind = Data.MeanVectorAP_individual(NC14:end,:,:);
    Fluo = nanmean(Fluo_ind, 3);
    Fluo_SEM = nanstd(Fluo_ind, 0,3)./sqrt(nEmbryos);
    
    Fraction_inst = Data.FractionON(NC14:end,:);
    Fraction_inst_SEM = nanstd(Data.FractionON_individual,0,3)./sqrt(length(Data.FractionON_individual(1,1,:)));
    Fraction_inst_SEM = Fraction_inst_SEM(NC14:end,:);

    % Runt null
    [~,~,nEmbryos_null] = size(Data_null.MeanVectorAP_individual);
    NC14_null = Data_null.nc14;
    Time_null = Data_null.ElapsedTime(NC14_null:end) - Data_null.ElapsedTime(NC14_null);
%     Fluo_null = Data_null.MeanVectorAP;
%     Fluo_SEM_null = Data_null.SEVectorAP;

    Fluo_ind_null = Data_null.MeanVectorAP_individual(NC14:end,:,:);
    Fluo_null = nanmean(Fluo_ind_null, 3);
    Fluo_SEM_null = nanstd(Fluo_ind_null, 0,3)./sqrt(nEmbryos_null);
    
    Fraction_inst_null = Data_null.FractionON;
    Fraction_inst_SEM_null = nanstd(Data_null.FractionON_individual,0,3)./sqrt(length(Data_null.FractionON_individual(1,1,:)));
    
    % for [000], [111] Runt nulls, we don't have Histone marker, thus
    % having the FractionON info as # of spots/APbin area. We need to
    % scale this to [0-1] scale as other datasets.
    if construct==1 || construct==4 % [000]/[111], Runt nulls 
        [row,col] = size(Fraction_inst_null);
        scale = prctile(reshape(Fraction_inst_null,[1,row*col]),99);
        Fraction_inst_null_rescaled = Fraction_inst_null./scale;
        Fraction_inst_null = Fraction_inst_null_rescaled;
        
        % individual fraction on scaling independently
        frac_individual = Data_null.FractionON_individual;
        numEmbryos = length(frac_individual(1,1,:));
        for i=1:numEmbryos
            frac_individual_rescaled(:,:,i) = frac_individual(:,:,i)./max(max(frac_individual(:,:,i)));
        end
        Fraction_inst_SEM_null = nanstd(frac_individual_rescaled,0,3)./sqrt(numEmbryos);
    else
    end

    hold on
    % Runt WT
    errorbar(Time(NC14:end) - Time(NC14),...
                Fraction_inst(NC14:end,APbin) ,...
                Fraction_inst_SEM(NC14:end,APbin),...
                'LineWidth',2,'Color',ColorChoice(construct,:))
    % Runt null        
    errorbar(Time_null(NC14_null:end) - Time_null(NC14_null),...
                Fraction_inst_null(NC14_null:end,APbin) ,...
                Fraction_inst_SEM_null(NC14_null:end,APbin),...
                'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2 )

    % xTicks, yTicks
    xlim([0 40])
    ylim([0 1.4])
    xticks([0 10 20 30 40])
    %yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('time into nc14 (min)')
    ylabel({'fraction of'; 'active nuclei'})

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_MS2, fig_MS2.CurrentAxes)
%    pause(1)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

%% averaged MS2 traces over "competent" nuclei(averaged over multiple embryos) for each construct with/without Runt protein

% Make a new directory
FigPath = [FigureRoot, filesep, 'MS2Traces_competentNuclei'];
mkdir(FigPath)

APaxis = 0:0.025:1;
% 2, 10th rows are [000], WT and Runt null
fig_MS2 = figure;

% Choose the APbin that we want to look at
APbin = 13; % 30% of the embryo length


for construct=1:8 % total number of constructs (enhancers)
    clf
    
    % Pull the data either for Runt WT or Runt null
    hold on
    for i=1:2
        
        % initialize the varialbes
        vars = {'Data','Time',...
            'Fluo_competent','Fluo_competent_SEM',...
            'NC14',...
            'Fluo_competent_ind', 'Fluo_ind','NParticles_ind',...
            'Ncompetent_ind','Data'};
        clear (vars{:})
        
        % Extract the NC14 traces out of the master structure
        if i==1
           Data = AveragedData{construct+1,2}; 
        else
            Data = AveragedData{construct+1+8,2};
        end
    
    
        % extract useful fields for plotting / process a bit
        [~,~,nEmbryos] = size(Data.MeanVectorAP_individual);
        NC14 = Data.nc14;
        Time = Data.ElapsedTime(NC14:end) - Data.ElapsedTime(NC14);

        % to calculate the averaged fluo over "competent" nuclei
        Fluo_ind = Data.MeanVectorAP_individual(NC14:end,:,:);
        NParticles_ind = Data.NParticlesAP_individual(NC14:end,:,:);
        % estimate the number of competent nuclei as 
        Ncompetent_ind = squeeze(max(NParticles_ind,[],1)); % APbins x embryos

        for embryo = 1:nEmbryos
            Fluo_competent_ind(:,:,embryo) = Fluo_ind(:,:,embryo).*...
                                            NParticles_ind(:,:,embryo)./...
                                            Ncompetent_ind(:,embryo)';
        end

        Fluo_competent = nanmean(Fluo_competent_ind,3);
        Fluo_competent_SEM = nanstd(Fluo_competent_ind,0,3)./sqrt(nEmbryos);
    
        if i==1
            % Runt WT
            errorbar(Time, Fluo_competent(:,APbin), Fluo_competent_SEM(:,APbin),...
                        'LineWidth',2,'Color',ColorChoice(construct,:))
        else
            % Runt null
            errorbar(Time, Fluo_competent(:,APbin), Fluo_competent_SEM(:,APbin),...
                        'LineWidth',2,'Color',(ColorChoice(construct,:)+[1 1 1])/2)
        end
    end


    % xTicks, yTicks
    xlim([0 40])
    ylim([0 1200])
    xticks([0 10 20 30 40])
    %yticks([0 100 200 300 400])

    %set(gca,'yticklabel',[])

    % no title, no-caps on the axis labels
    xlabel('time into nc14 (min)')
    ylabel('mean fluorescence (AU)')

    legend(constructNames{construct},constructNames{construct+8},'Location','NorthEast')

    box on

    StandardFigure(fig_MS2, fig_MS2.CurrentAxes)
%    pause(1)

    % Save the plot
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.tif']); 
    saveas(gcf,[FigPath,filesep,constructNames{construct},'.pdf']); 
end

% Note that the [000], Runt null should be pushed 2 minutes
%% plot for check (individual mean over competent nuclei vs averaged over ON nuclei)
% APbin = 13; % 30%
% hold on
% % for embryo = 1:nEmbryos
% %     plot(Time, Fluo_competent(:,APbin,embryo))
% % end
% errorbar(Time, Fluo(:,APbin), Fluo_SEM(:,APbin))
% errorbar(Time, Fluo_competent(:,APbin), Fluo_competent_SEM(:,APbin))
% 
% 
% xlabel('time into nc14 (min)')
% ylabel('fluorescence (AU)')
% 
% legend('instantaneous ON','Competent')
% 
% StandardFigure(gcf,gca)
%% averaged MS2 time traces from a single embryo (averaged over all nuclei in that APbin)
% Note. This is just to get a sense of how trend looks like, without
% worrying about the synchronization issue that could arise from a poor
% mitosis registration at the "CheckDivisionTimes.m"
