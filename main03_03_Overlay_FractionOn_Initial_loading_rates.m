function main02_06_Overlay_FractionOn_Initial_loading_rates
% DECRIPTION 
% This script is for generating an overlaid plot of
% Fraction ON and loading rates profile across AP. 

% Caveats
% 1) loading rates could be (1) initial loading rate (either asymmetric or
% the linear fit), (2) mean loading rate over one nuclear cycle
% 2) This is loading pre-saved data from other scripts so that we can
% change the data as we update things.

%% Load the datasets
% Load the Fraction ON
DataPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\OpposingGradients_ProcessedData';

FractionON = load([DataPath,filesep,'FractionON_r0123.mat']);

FractionON_Average = FractionON.FractionON_Average;
FractionON_SEM = FractionON.FractionON_SEM;
FractionON_individual = FractionON.FractionON_individual;


%% Load the initial loading rate (Asymmetric fit over ON nuclei)
InitialRate_Asymmetric = load([DataPath,filesep,'AveragedInitialRate_AsymmetricFit_FixedFittingScript.mat']);

% Assign the fields from the loaded (processed) data
InitialRate_r0 = InitialRate_Asymmetric.average_fittedRate_r0;
InitialRate_SEM_r0 = InitialRate_Asymmetric.SEM_fittedRate_r0;

InitialRate_r1 = InitialRate_Asymmetric.average_fittedRate_r1;
InitialRate_SEM_r1 = InitialRate_Asymmetric.SEM_fittedRate_r1;

InitialRate_r2 = InitialRate_Asymmetric.average_fittedRate_r2;
InitialRate_SEM_r2 = InitialRate_Asymmetric.SEM_fittedRate_r2;

InitialRate_r3 = InitialRate_Asymmetric.average_fittedRate_r3;
InitialRate_SEM_r3 = InitialRate_Asymmetric.SEM_fittedRate_r3;

%% Optional : Normalization of the Initial rate
% What's the maximum value to divide with? I'll pick up the maximum of the
% initial rate among NC13 and NC14 of all constructs.

%% Color definition
% This is defining the line color
colorDict = struct();
colorDict.blue = [115,143,193]/255; %[115,143,170]/255;
colorDict.red =  [213,108,85]/255; %[200,108,85]/255;
colorDict.yellow = [234,194,100]/255;
colorDict.cyan = [108,188,233]/255;
colorDict.magenta = [208,109,171]/255;
colorDict.lightBlue = [115,142,193]/255;
colorDict.purple = [171,133,172]/255;
colorDict.green =  [122,169,116]/255; %[122,150,116]/255;
colorDict.brown = [179,155,142]/255;
colorDict.darkgreen = [126,157,144]/255;

ColorChoice = [colorDict.magenta; colorDict.lightBlue; colorDict.yellow; colorDict.red; colorDict.brown]; % 4 embryos max. it could be extended easily
%% Plot (Overlay)
% Plot NC13 and NC14 separately.
APaxis = 0:0.025:1;

for NC=13:14
    FractionON_figure(NC-12) = figure(NC-12);
    hold on
    % Loop over different constructs (r0,1,2,3,3')
    yyaxis left
    for i=1:length(FractionON_individual) 
        % Extract the Fraction ON from the structure
        FractionON_individual_temp = FractionON_individual{i};
        FractionON_Average_temp = FractionON_Average{i};
        FractionON_SEM_temp = FractionON_SEM{i};
        
        [~,~,numEmbryos] = size(FractionON_individual_temp);

        errorbar(0:0.025:1,FractionON_Average_temp(:,NC-12),FractionON_SEM_temp(:,NC-12),'Color',ColorChoice(i,:))
    end
    title(['Overlaid ','Fraction ON and Initial loading rate',' @ NC ',num2str(NC)])
    xlabel('AP (EL)')
    ylabel('Fraction ON')
    xlim([0.2 0.8])
    ylim([0 1.2])
    
    pause
    yyaxis right
    nc = NC-12;
    % r0,1,2,3 altogether
    hold on
    errorbar(0:0.025:1,InitialRate_r0(:,nc),InitialRate_SEM_r0(:,nc),'--','Color',ColorChoice(1,:))
    errorbar(0:0.025:1,InitialRate_r1(:,nc),InitialRate_SEM_r1(:,nc),'--','Color',ColorChoice(2,:))
    errorbar(0:0.025:1,InitialRate_r2(:,nc),InitialRate_SEM_r2(:,nc),'--','Color',ColorChoice(3,:))
    errorbar(0:0.025:1,InitialRate_r3(:,nc),InitialRate_SEM_r3(:,nc),'--','Color',ColorChoice(4,:))

    xlim([0.15 0.8])
    ylim([0 400])

    legend('r0','r1','r2','r3')
    xlabel('AP Position')
    ylabel('Initial rate (AU/min)')
    
    hold off
    StandardFigure(FractionON_figure(NC-12),FractionON_figure(NC-12).CurrentAxes)

end
end