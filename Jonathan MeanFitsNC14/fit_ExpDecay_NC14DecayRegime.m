%% Script to fit an exponential fit to the accumulated mRNA in the decay regime of NC14
function fit_ExpDecay_NC14DecayRegime(varargin)
% This script is basically to regenerate the FigS4 in Garcia, 2013 paper,
% where they attempted to extract the turn-off time by fitting an
% exponential decay to the accumulated mRNA from the Time_peak.
% Last updated : 9/1/2020, Yang Joon Kim

% INPUT : Prefix
% OUTPUT : mRNA_max, mRNA_max_SD, Tau, Tau_SD, T_peak saved in
% MeanFitsAsymmetric.mat

% OPTIONS
% 'TimeWindow',[Time1 Time2]

% Workflow
% 0) Load the CompiledParticles (MeanVectorAP, etc.) using the Prefix
% 1) Find the timepoint when the fluorescence peaks.
% 2) Calculate the accmulated mRNA from that time point onward.
% 3) Fit that curve with the equation, mRNA(t) = mRNA_max * (1- exp(-(t-t_peak)/Tau))
%   where the parameters are mRNA_max (maximum accumulated mRNA) and the decay constant, Tau.

% Let's use the same language as in FitMeanAPAsymmetric.m such that we can
% estimate the standard deviation from this fitting. (lsqnonlin)
% Also, save the output (inferred parameters) : mRNA_max, mRNA_max_SD, Tau,
% Tau_SD, in MeanFitsAsymmetric.mat file in the Result folder as well, so
% that we can easily load initial slope, and Tau (as well as the mRNA_max) 

close all


%% Get the default folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders;

% Options
checkFits = false;

for i=1:length(varargin)
    Prefix=varargin{1};
    if strcmpi(varargin{i},'TimeWindow')
        TimeWindow = varargin{i+1}; %[t1 t2] in nc14
    elseif strcmpi(varargin{i},'CheckFitting')
        checkFits = true;
    end
end

%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);

        


%% Load the complied particles and the division information                                    
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'])

if exist([DropboxFolder,filesep,Prefix,'\APDivision.mat'], 'file')
    load([DropboxFolder,filesep,Prefix,'\APDivision.mat'], 'APDivision')
else
    error('Could not load APDivision.mat. Make sure to have done the manual check of division.')
end


%%  Extract the fields from the cell structure (This is for fields like MeanVectorAP
% that are saved inside {}.
channel = 1;

if iscell(MeanVectorAP)
    MeanVectorAP = MeanVectorAP{channel};
    SDVectorAP = SDVectorAP{channel};
    NParticlesAP = NParticlesAP{channel};
end

%% Load the MeanFitsAsymmetric.mat file from FitMeanAPAsymmetric.m
if exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat']);
    if isempty(FitResults)
        warning('MeanFitsAsymmetric not found. Have you done the fitting?')
    end
elseif exist([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV2.mat'])
    load([DropboxFolder,filesep,Prefix,filesep,'MeanFitsV2.mat']);
    if isempty(FitResults)
        warning('MeanFitsAsymmetric not found. Have you done the fitting?')
    end
else
    warning('MeanFitsAsymmetric not found. Have you done the fitting?')
end


%% Extract info about whether each AP bins is approved or not.

Approved = zeros(length(APbinID),1);

for i=1:length(APbinID)
    Approved(i,1) = FitResults(i,3).Approved;
end

ApprovedAPbins = find(Approved==1);
%% Go through all (approved) AP bins and follow the workflow described above.

% First, find the first, and the last non-zero APbins
% APbinStart=min(find(sum(NParticlesAP)));
% APbinEnd=max(find(sum(NParticlesAP)));

numAPbins = length(APbinID);

% Initialize the parameters as NaNs 
% maximum value of integrated fluo
IntFluo_max = nan(numAPbins,1);
IntFluo_max_SD = nan(numAPbins,1);
% decay constant (tau)
Tau = nan(numAPbins,1);
Tau_SD = nan(numAPbins,1);

% time window and integrated fluo
%tWindow = nc14:length(ElapsedTime);
%int_fluo = zeros(numAPbins, length(tWindow));

% Time point at which the fluorescence peaks
Time_peak = nan(numAPbins,1);
Time_peak_indices = nan(numAPbins,1);

for AP = ApprovedAPbins(1):ApprovedAPbins(end) %APbinStart:APbinEnd
    
    % nc14
    NC14 = APDivision(14,AP);
    
    % define the Fluo and Time
    Time = ElapsedTime(NC14:length(ElapsedTime)) - ElapsedTime(NC14); % initialize to 0 min
    Fluo = MeanVectorAP(NC14:length(ElapsedTime), AP)';
    %Fluo_SEM = SDVectorAP(nc14:length(ElapsedTime), AP)./NParticlesAP(nc14:length(ElapsedTime), AP);
    
    % filter out NaNs from the Fluo
    Fluo(isnan(Fluo)) = 0; % making as zero is fair as it doesn't contribute to the integrated amount.
    
    if sum(Fluo)~=0
    
        % find the peak of the Fluo (within 5min-15min into nc14, just to ignore some weird spikes before or after)
        index_peak = find(Fluo == max(Fluo(8:17))); %8, 25 are hard-coded to account for 40sec/frame. This shoudl be relaxed to be more flexible.
        Time_peak(AP,1) = Time(index_peak);
        Time_peak_indices(AP,1) = index_peak;

        % calculate the integrated fluo from this T_peak (decay regime)
        % first, initialize the inf_fluo as the time length can be different
        % for APbins (due to the fact that we're using the APDivision.mat)
        int_fluo = zeros(1,length(Time));
        for t = index_peak+1:length(Time)
            int_fluo(t) = trapz(Time(index_peak:t), Fluo(index_peak:t));
        end

        % fit with an exponential curve to extract 1)maximum value, 2)decay
        % constant
        % x = [IntFluo_max Tau]

        % Note that I'm fitting only from the T_peak either 20 min of NC14 (or
        % end of our measurement)
        time_interval = median(diff(ElapsedTime));
        tEnd = TimeWindow(2);%20; % min

        tEnd_index = min(floor(tEnd/time_interval), length(Time)); % counted from the beginning of nc14
        tWindow_fit_index = [index_peak, tEnd_index];

        % interpolation for a smoother fitting
        t_fold = 10;
        Time_fit = Time(index_peak):time_interval/t_fold:Time(tEnd_index);
        int_fluo_fit = interp1(Time, int_fluo, Time_fit);

        % fit
        fun = @(x) x(1)*(1 - exp(-(Time_fit-Time(index_peak))/x(2))) - int_fluo_fit;
        %fun = @(x) x(1)*(1 - exp(-(Time(tWindow_fit_index)-Time(index_peak))/x(2))) - int_fluo(tWindow_fit_index,AP)';
        x0 = [10000, 5];
    %     lb = [max(int_fluo(:,AP)), 0];
    %     ub = [3*max(int_fluo(:,AP)),30];
        [xFit,resnorm,residual,exitflag,output,lambda,jacobian] =  lsqnonlin(fun,x0); %, lb, ub);

        IntFluo_max(AP,1) = xFit(1);
        Tau(AP,1) = xFit(2);

        % Calculate the Standard deviation using the confidence interval
        try
            CI=nlparci(xFit,residual,'jacobian',jacobian,'alpha',0.68);

            IntFluo_max_SD(AP,1) = (CI(1,2)-CI(1,1))/2;
            Tau_SD(AP,1) = (CI(2,2)-CI(2,1))/2;
        catch
            warning([num2str((AP-1)*2.5),'% bin is not returning the CI. Check the plot'])
        end
    else
    end
end

%% Save the fitted result into a structure, then save into the MeanAPAsymmetric.mat
% structure : NC14DecayRegimeFitResults
% fields : IntFluo_max, IntFluo_max_SD, Tau, Tau_SD
NC14DecayRegimeFitResults.IntFluo_max = IntFluo_max;
NC14DecayRegimeFitResults.IntFluo_max_SD = IntFluo_max_SD;

NC14DecayRegimeFitResults.Tau = Tau;
NC14DecayRegimeFitResults.Tau_SD = Tau_SD;

NC14DecayRegimeFitResults.Time_peak = Time_peak; % peak time

% Approval/Disapproval
% initialize the approval/disapproval states as zeros (so that this NC14
% decay regime fit is independent of initial slope fitting. For duration
% calculation, we'd need to pick embryos, APbins where both are approved.)

% for manual checking, we need to initialize this
% Approved = zeros(numAPbins,1); 
NC14DecayRegimeFitResults.Approved = Approved; % default is from the FitMeanAPAsymmetric.m result

% %Save the information
% save([DropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat'],...
%     'FitResults','NC14DecayRegimeFitResults')
% display('MeanFitsAsymmetric.mat updated')   

%% Plot for checking if the fitting was reasonable
% generate plot
% for AP=ApprovedAPbins(1):ApprovedAPbins(end)
%     clf
%     
%     NC14 = APDivision(14,AP);
%     Time =  ElapsedTime(NC14:length(ElapsedTime)) - ElapsedTime(NC14);
%     Fluo = MeanVectorAP(NC14:length(ElapsedTime), AP)';
%     Fluo_SEM = SDVectorAP(NC14:length(ElapsedTime), AP)./NParticlesAP(NC14:length(ElapsedTime), AP);
%     
%     % Convert the NaNs in fluo to zeros
%     Fluo(isnan(Fluo)) = 0;
%     
%     if sum(Fluo)~=0 
%         % Calculate the integrated fluo again, as the length is different for
%         % different AP bins
%         index_peak = Time_peak_indices(AP,1);
%         int_fluo = zeros(1,length(Time));
%         for t = index_peak+1:length(Time)
%             int_fluo(t) = trapz(Time(index_peak:t), Fluo(index_peak:t));
%         end
% 
%         hold on
%         % MS2 spot fluo plot
%         yyaxis left
%         h(1) = errorbar(Time, Fluo, Fluo_SEM)
%         h(2) = xline(Time(Time_peak_indices(AP)),'--')
%         ylim([0 max(Fluo)*1.4])
%         ylabel('mean fluorescence (AU)')
% 
%         % Integrated fluo plot
%         yyaxis right
%         % integrate fluo(data)
%         h(3) = plot(Time, int_fluo)
%         % fit
%         h(4) = plot(Time(Time_peak_indices(AP):end), IntFluo_max(AP)*(1 - exp(-(Time(Time_peak_indices(AP):end)-Time(Time_peak_indices(AP)))/Tau(AP))))
%         % marking the Tau from T_peak
%         h(5) = xline(Time(Time_peak_indices(AP)) + Tau(AP),'--')
%         ylabel('integrated fluorescence (AU*min)')
%         % xTicks, yTicks
%         xlim([0 60])
%         %ylim([0 IntFluo_max(AP)*1.4])
%         try
%             ylim([0 max(int_fluo)*1.4])
%         catch
%             ylim([0 IntFluo_max(AP)*1.4])
%         end
%         %xticks([0 10 20 30 40 50])
% 
%         % set(gca,'yticklabel',[])
% 
%         % no title, no-caps on the axis labels
%         xlabel('time into nc14 (min)')
%         title(['AP=',num2str((AP-1)*2.5),'%'])
%         % xticks([Time(Time_peak_indices(AP)) Time(Time_peak_indices(AP)) + Tau(AP)])
%         % xticklabels({'T_{peak}','T_{off}'})
% 
%         legend([h(1) h(3) h(4)],'MS2','integrated','Fit','Location','NorthEast')
% 
%         % StandardFigurePBoC(fig_name, fig_name.CurrentAxes)
%         StandardFigurePBoC([h(1) h(3) h(4)],gca)
%         pause(0.5)
%     else
%     end
% end
%% Make a while loop to approve or disapprove the fitted result for each AP bin
% initialize the conditions
cc=1;
% Start with the first AP bin
AP=ApprovedAPbins(1);

% file path
figPath = [DropboxFolder,filesep,Prefix,filesep,'NC14DecayRegimeFits']
mkdir(figPath)

% if CheckFits is true, then do a while loop for Approval/Disapproval of
% fits
if checkFits
    while (cc~=13)
        %% figure
        clf

        % set tje background color for approval/disapproval
        if Approved(AP)==-1
            set(gcf,'Color','r')
        elseif Approved(AP)==1
            set(gcf,'Color','g')
        else
            set(gcf,'Color','default')
        end
        %% plot to check the integrated fluo and fit (Tau)

        %for AP=ApprovedAPbins(1):ApprovedAPbins(end)
        % NC14 = APDivision(14,AP);
        Time =  ElapsedTime(nc14:length(ElapsedTime)) - ElapsedTime(nc14);
        Fluo = MeanVectorAP(nc14:length(ElapsedTime), AP)';
        Fluo_SEM = SDVectorAP(nc14:length(ElapsedTime), AP)./NParticlesAP(nc14:length(ElapsedTime), AP);

        % generate plot
        hold on
        % MS2 spot fluo plot
        yyaxis left
        h(1) = errorbar(Time, Fluo, Fluo_SEM)
        h(2) = xline(Time(Time_peak_indices(AP)),'--')
        ylim([0 max(Fluo)*1.4])
        ylabel('mean fluorescence (AU)')
        % Integrated fluo plot
        yyaxis right
        h(3) = plot(Time, int_fluo(:,AP))
        h(4) = plot(Time(Time_peak_indices(AP):end), IntFluo_max(AP)*(1 - exp(-(Time(Time_peak_indices(AP):end)-Time(Time_peak_indices(AP)))/Tau(AP))))
        h(5) = xline(Time(Time_peak_indices(AP)) + Tau(AP),'--')
        ylabel('integrated fluorescence (AU*min)')
        % xTicks, yTicks
        xlim([0 60])
        ylim([0 IntFluo_max(AP)*1.4])
        %xticks([0 10 20 30 40 50])

        % set(gca,'yticklabel',[])

        % no title, no-caps on the axis labels
        xlabel('time into nc14 (min)')
        title(['AP=',num2str((AP-1)*2.5),'%'])
        % xticks([Time(Time_peak_indices(AP)) Time(Time_peak_indices(AP)) + Tau(AP)])
        % xticklabels({'T_{peak}','T_{off}'})

        legend([h(1) h(3) h(4)],'MS2','integrated','Fit','Location','NorthEast')

        % StandardFigurePBoC(fig_name, fig_name.CurrentAxes)
        StandardFigurePBoC([h(1) h(3) h(4)],gca)

        % save the plot after approval/disapproval
        %end

        %% Keypad for options (move between APbins, Approval/Disapproval of fits, etc.)
        % To do : Add functionality to change the fitting range, or parameters,
        % or constrain these further.

        ct=waitforbuttonpress;
        cc=get(gcf,'currentcharacter');
        cm=get(gca,'CurrentPoint');

        %Move between AP positions
        if (ct~=0)&(cc=='.')&(AP<ApprovedAPbins(end)) & ~isnan(Time_peak_indices(AP))
            AP = AP+1;
        elseif (ct~=0)&(cc==',')&(AP>ApprovedAPbins(1))
            AP=AP-1;

        %Approve, disapprove fit
        elseif (ct~=0)&(cc=='q')
            if NC14DecayRegimeFitResults.Approved(AP)==0
                NC14DecayRegimeFitResults.Approved(AP)=1;
                % Save the plot
                title(['Approved, ','AP=',num2str((AP-1)*2.5),'%'])
                saveas(gcf,[figPath,filesep,num2str((AP-1)*2.5),'%','.tif']); 
    %             saveas(gcf,[figPath,filesep,num2str((AP-1)*2.5),'%','.pdf']); 
    %             saveas(gcf,[figPath,filesep,num2str((AP-1)*2.5),'%','.png']); 
            elseif NC14DecayRegimeFitResults.Approved(AP)==1
                NC14DecayRegimeFitResults.Approved(AP)=1;
            end


        %Disapprove, disapprove fit
        elseif (ct~=0)&(cc=='w')
            if NC14DecayRegimeFitResults.Approved(AP)==0
                NC14DecayRegimeFitResults.Approved(AP)=-1;
                % Save the plot
                title(['Disapproved, ','AP=',num2str((AP-1)*2.5),'%'])
                saveas(gcf,[figPath,filesep,num2str((AP-1)*2.5),'%','.tif']); 
    %             saveas(gcf,[figPath,filesep,num2str((AP-1)*2.5),'%','.pdf']); 
    %             saveas(gcf,[figPath,filesep,num2str((AP-1)*2.5),'%','.png']); 
            elseif NC14DecayRegimeFitResults.Approved(AP)==-1
                NC14DecayRegimeFitResults.Approved(AP)=-1;
            end

        %Save
        elseif (ct~=0)&(cc=='e')
            save([DropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat'],...
                'FitResults','NC14DecayRegimeFitResults')
        display('MeanFitsAsymmetric.mat updated')
        end
    end
else
    % just show the plots?
    for AP=ApprovedAPbins(1):ApprovedAPbins(end)
        clf

        NC14 = APDivision(14,AP);
        Time =  ElapsedTime(NC14:length(ElapsedTime)) - ElapsedTime(NC14);
        Fluo = MeanVectorAP(NC14:length(ElapsedTime), AP)';
        Fluo_SEM = SDVectorAP(NC14:length(ElapsedTime), AP)./NParticlesAP(NC14:length(ElapsedTime), AP);

        % Convert the NaNs in fluo to zeros
        Fluo(isnan(Fluo)) = 0;

        if sum(Fluo)~=0 
            % Calculate the integrated fluo again, as the length is different for
            % different AP bins
            index_peak = Time_peak_indices(AP,1);
            int_fluo = zeros(1,length(Time));
            for t = index_peak+1:length(Time)
                int_fluo(t) = trapz(Time(index_peak:t), Fluo(index_peak:t));
            end

            hold on
            % MS2 spot fluo plot
            yyaxis left
            h(1) = errorbar(Time, Fluo, Fluo_SEM)
            h(2) = xline(Time(Time_peak_indices(AP)),'--')
            ylim([0 max(Fluo)*1.4])
            ylabel('mean fluorescence (AU)')

            % Integrated fluo plot
            yyaxis right
            % integrate fluo(data)
            h(3) = plot(Time, int_fluo)
            % fit
            h(4) = plot(Time(Time_peak_indices(AP):end), IntFluo_max(AP)*(1 - exp(-(Time(Time_peak_indices(AP):end)-Time(Time_peak_indices(AP)))/Tau(AP))))
            % marking the Tau from T_peak
            h(5) = xline(Time(Time_peak_indices(AP)) + Tau(AP),'--')
            ylabel('integrated fluorescence (AU*min)')
            % xTicks, yTicks
            xlim([0 60])
            %ylim([0 IntFluo_max(AP)*1.4])
            try
                ylim([0 max(int_fluo)*1.4])
            catch
                ylim([0 IntFluo_max(AP)*1.4])
            end
            %xticks([0 10 20 30 40 50])

            % set(gca,'yticklabel',[])

            % no title, no-caps on the axis labels
            xlabel('time into nc14 (min)')
            title(['AP=',num2str((AP-1)*2.5),'%'])
            % xticks([Time(Time_peak_indices(AP)) Time(Time_peak_indices(AP)) + Tau(AP)])
            % xticklabels({'T_{peak}','T_{off}'})

            legend([h(1) h(3) h(4)],'MS2','integrated','Fit','Location','NorthEast')

            % StandardFigurePBoC(fig_name, fig_name.CurrentAxes)
            StandardFigurePBoC([h(1) h(3) h(4)],gca)
            pause(0.5)
        else
        end
    end
end
%% Save the information
save([DropboxFolder,filesep,Prefix,filesep,'MeanFitsAsymmetric.mat'],...
    'FitResults','NC14DecayRegimeFitResults')
display('MeanFitsAsymmetric.mat updated')   

close all
%% generate individual plots
% % select an AP bin
% AP = 9;
% 
% % NC14 = APDivision(14,AP);
% Time =  ElapsedTime(nc14:length(ElapsedTime)) - ElapsedTime(nc14);
% Fluo = MeanVectorAP(nc14:length(ElapsedTime), AP)';
% Fluo_SEM = SDVectorAP(nc14:length(ElapsedTime), AP)./NParticlesAP(nc14:length(ElapsedTime), AP);
% 
% traceFig = figure;
% hold on
% % MS2 spot fluo plot
% h(1) = errorbar(Time, Fluo, Fluo_SEM)
% h(2) = xline(Time(Time_peak_indices(AP)),'--')
% xlim([0 40])
% ylim([0 max(Fluo)*1.4])
% 
% yticks([0 200 400 600 800 1000])
% 
% ylabel('mean fluorescence (AU)')
% xlabel('time into nc14 (min)')
% hold off
% 
% legend([ h(1)],'MS2','Location','NorthEast')
% % StandardFigurePBoC(fig_name, fig_name.CurrentAxes)
% StandardFigurePBoC([h(1)],gca)
% % Save the plot
% figPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\Data\NC14_DecayRegime_fitting';
% saveas(traceFig,[figPath,filesep,'example_r2far_set1_20%_MS2trace','.tif']); 
% saveas(traceFig,[figPath,filesep,'example_r2far_set1_20%_MS2trace','.pdf']); 
% 
% 
% % Integrated fluo plot
% integratedFig = figure;
% hold on
% h(2) = xline(Time(Time_peak_indices(AP)),'--')
% h(3) = plot(Time, int_fluo(:,AP))
% h(4) = plot(Time(Time_peak_indices(AP):end), IntFluo_max(AP)*(1 - exp(-(Time(Time_peak_indices(AP):end)-Time(Time_peak_indices(AP)))/Tau(AP))),'--')
% h(5) = xline(Time(Time_peak_indices(AP)) + Tau(AP),'--')
% h(6) = yline(IntFluo_max(AP),'--')
% ylabel({'mRNA produced during'; 'decay regime (AU)'})
% % xTicks, yTicks
% xlim([0 40])
% ylim([0 IntFluo_max(AP)*1.4])
% 
% xlabel('time into nc14 (min)')
% hold off
% 
% yticks([0 2000 4000 6000 8000 10000])
% 
% legend([ h(3) h(4)],'integrated','Fit','Location','NorthEast')
% %StandardFigure(gcf, gca)
% StandardFigurePBoC([h(3) h(4)],gca)
% % Save the plot
% figPath = 'S:\YangJoon\Dropbox\OpposingGradientsFigures\Data\NC14_DecayRegime_fitting';
% saveas(integratedFig,[figPath,filesep,'example_r2far_set1_20%_integratedmRNA_fit','.tif']); 
% saveas(integratedFig,[figPath,filesep,'example_r2far_set1_20%_integratedmRNA_fit','.pdf']); 

%% Additional quality control scripts
%% plot the turn-on time and turn-off time
% Note that the turn-on time should be extracted from the
% MeanFitsAsymmetric.mat file.
% APaxis = 0:0.025:1;
% 
% hold on
% errorbar(APaxis, Time_peak + Tau,Tau_SD,'o')

%% relative error in fitting
% plot(APaxis, Tau_SD./Tau,'o')
% ylim([0 0.1])
% xlabel('embryo length')
% ylabel('relative error')

%%
end