function FitMeanAPTiltedTrapezoid_NC13Only(Prefix)
%Yang Joon and Paul, 7/1/2018
%This code is modified from Fit_Initial_Rates to fit the initialRates of
%the averaged DataSets, and also the tilted trapezoid as three different
%rates.

% Description 
% You have to click the min point then the maximum for initial rate
% Same commands than FitMeanAPSymetric, but 'q' changes from Approved to
% Disapproved, 'f' is for doing the fit again.

% Load the input dataset (averaged over multiple embryos)
%Change the following line to change DataSet : 
%Construct = 'r3';
%File = strcat('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\',Construct,'_FromNC12.mat');

%load(File);
%% Build The FitResults Structure

%Parameters:
MinParticles=1;     %Minimum number of particles in an AP bin
MinTimePoints=5;    %Minimum number of time points where we'll have at least
                    %the minimum number of particles.
ElongationRate=1.54;    %In kb/minutes.
GeneLength=5.296;       %Distance from the first MS2 site to the end of the
                        %TUB3'UTR in kb.
Delay=GeneLength/ElongationRate;    %Minutes for PolII to fall off after reaching
                                    %the first MS2 site.

                                    
close all


%Get the default folders
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders;

% if ~isempty(varargin)
%     Prefix=varargin{1};
%                
% else
%     FolderTemp=uigetdir(DropboxFolder,'Choose folder with files to analyze');
%     Dashes=strfind(FolderTemp,'\');
%     Prefix=FolderTemp((Dashes(end)+1):end);
% end

%Get the relevant folders now:
[SourcePath,FISHPath,DropboxFolder,MS2CodePath]=...
    DetermineLocalFolders(Prefix);
        


%Load the complied particles and the division information                                    
load([DropboxFolder,filesep,Prefix,'\CompiledParticles.mat'])

if exist([DropboxFolder,filesep,Prefix,'\APDivision.mat'])
    load([DropboxFolder,filesep,Prefix,'\APDivision.mat'])
else
    error('Could not load APDivision.mat. Make sure to have done the manual check of division.')
end

%File2 = strcat('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\',Construct,'_InitialFit.mat');

try
    load([DropboxFolder,filesep,Prefix,'\MeanFitsV2.mat'])

catch
    warning('No prior Fit Results exists')
    if (exist('FitResults', 'var' )) == 0
        for i = 1:41
            for j=1:3
                % Think about three rates, initial rise, tilted decrease,
                % and stop (RNAP falling off)
                FitResults(i,j).RateFit1 = [];
                FitResults(i,j).SDRateFit1 = [];
                FitResults(i,j).RateFit2 = [];
                FitResults(i,j).SDRateFit2 = [];
                FitResults(i,j).RateFit3 = [];
                FitResults(i,j).SDRateFit3 = [];
                
                FitResults(i,j).TimeStart = [];
                FitResults(i,j).SDTimeStart = [];
                FitResults(i,j).TimePeak = [];
                FitResults(i,j).SDTimePeak = [];
                FitResults(i,j).Approved = -1;
                FitResults(i,j).X=[0,0];
                FitResults(i,j).Y=[0,0];
            end
        end
    end
end

% YJK, 2018-12-00. I noticed that some datasets have MeanVectorAP,
% SDVectorAP, and NParticlesAP in a cell (1x1). For those, we should grab
% them out.

if iscell(MeanVectorAP) || iscell(SDVectorAP) || iscell(NParticlesAP)
    MeanVectorAP = cell2mat(MeanVectorAP);
    SDVectorAP = cell2mat(SDVectorAP);
    %SDVectorAP = cell2mat(SDVectorAP);
    NParticlesAP = cell2mat(NParticlesAP);
end


%% Initiations of parameters

%Parameters:
MinParticles=1;     %Minimum number of particles in an AP bin
MinTimePoints=5;    %Minimum number of time points where we'll have at least
%the minimum number of particles.


close all

CurrentNC=13;
X=zeros(41,3,2);
Y=zeros(41,3,2);
nc(1)=nc12;
nc(2)=nc13;
nc(3)=nc14;

%Go through each AP bin
FitFigure=figure;
i=min(find(sum(NParticlesAP))); %index of the first AP bin that has a non-zero number of particles
cc=1;


%% While Loop

while (cc~='x')
    
    figure(FitFigure)
    clf
    % Display of figures
    
    if FitResults(i,CurrentNC-11).Approved==-1
        set(gcf,'Color','r')
    elseif FitResults(i,CurrentNC-11).Approved==1
        set(gcf,'Color','g')
    else
        set(gcf,'Color','default')
    end
    
    % Definition of TimeWindows for different NC
    if CurrentNC==12
        FrameWindow=nc12:nc13;
    elseif CurrentNC==13
        FrameWindow=nc13:nc14;
    else
        FrameWindow=nc14:length(ElapsedTime);
    end
    
    
    %% Plots of the Data 
    % YJK : This part needs to be changed in that SEM as SDVectorAP /
    % sqrt(Number of embryos), which is saved in the SDVectorAP field.
    PlotHandle=errorbar(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
        MeanVectorAP(FrameWindow,i),...
        SDVectorAP(FrameWindow,i),'.-k');
    hold on
    
    PlotHandle(end+1)=plot(ElapsedTime(FrameWindow)-ElapsedTime(FrameWindow(1)),...
        MeanVectorAP(FrameWindow,i),'or');
    
    title([num2str(2.5*i),'% AP',...
        ', nc',num2str(CurrentNC)])
    
    %Set the limits on the x-axis for NC 14
    if CurrentNC==14
        xlim([0,60])
    end
    
    ylabel('Mean fluorescence nucleus')
    xlabel('Time into nc (min)')
    
    %% Ask to click four points (for three phases) to set the initial slope, approve.
    
    if (FitResults(i,CurrentNC-11).X==[0,0] & (sum(NParticlesAP(FrameWindow,i))~=0))
        [FitResults(i,CurrentNC-11).X,FitResults(i,CurrentNC-11).Y]=ginput(4);
        FitResults(i,CurrentNC-11).Approved=1;
        set(gcf,'Color','g');

        
    end
    
    %Change if someone took the top for the bottom, etc
    if FitResults(i,CurrentNC-11).X(2)<FitResults(i,CurrentNC-11).X(1)
        X1=FitResults(i,CurrentNC-11).X(2);
        FitResults(i,CurrentNC-11).X(2) = FitResults(i,CurrentNC-11).X(1);
        FitResults(i,CurrentNC-11).X(1) = X1;
        Y1=FitResults(i,CurrentNC-11).Y(2);
        FitResults(i,CurrentNC-11).Y(2) = FitResults(i,CurrentNC-11).Y(1);
        FitResults(i,CurrentNC-11).Y(1) = Y1;
    end
    
        
    Xplot=squeeze(FitResults(i,CurrentNC-11).X);
    Yplot=squeeze(FitResults(i,CurrentNC-11).Y);
    PlotHandle(end+1)=plot(Xplot,Yplot,'r','LineWidth',1.5);
    
     
    
    a = (Yplot(2)-Yplot(1))./(Xplot(2)-Xplot(1));
    b = Yplot(1)-a*Xplot(1);
    T_on = -(b./a);
    
    c = (Yplot(3)-Yplot(2))./(Xplot(3)-Xplot(2)); % rate at the tilted phase
    T_peak = Xplot(2);
    d = (Yplot(4)-Yplot(3))./(Xplot(4)-Xplot(3));
    
    
    FitResults(i,CurrentNC-11).RateFit1 = a;
    FitResults(i,CurrentNC-11).TimeStart = T_on;
    FitResults(i,CurrentNC-11).RateFit2 = c;
    FitResults(i,CurrentNC-11).TimePeak = T_peak;
    FitResults(i,CurrentNC-11).RateFit3 = d;
    %FitResults(i,CurrentNC-11).SDTimeStart = ? ;
    %FitResults(i,CurrentNC-11).SDRateFit = ? ; 
    % The points selected appear in red
    
    Range=FrameWindow;
    Range=Range((ElapsedTime(Range)-ElapsedTime(nc(CurrentNC-11)))>Xplot(1)-0.1);
    Range=Range((ElapsedTime(Range)-ElapsedTime(nc(CurrentNC-11)))<Xplot(2)+0.1);
    
    Fluo=MeanVectorAP(FrameWindow,i);
    Fluo=Fluo(Range-nc(CurrentNC-11)+1);
    
    PlotHandle(end+1)=plot(ElapsedTime(Range)-ElapsedTime(FrameWindow(1)),...
        Fluo,'or','MarkerFaceColor','r');
    
    legend(['tON = ',num2str(FitResults(i,CurrentNC-11).TimeStart),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimeStart)],...
        ['t_{peak} = ',num2str(FitResults(i,CurrentNC-11).TimePeak),' \pm ',num2str(FitResults(i,CurrentNC-11).SDTimePeak)],...
        ['Rate1 = ',num2str(FitResults(i,CurrentNC-11).RateFit1),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateFit1)],...
        ['Rate2 = ',num2str(FitResults(i,CurrentNC-11).RateFit2),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateFit2)],...
        ['Rate3 = ',num2str(FitResults(i,CurrentNC-11).RateFit3),' \pm ',num2str(FitResults(i,CurrentNC-11).SDRateFit3)],...
        'Location','SouthOutside')
    
    
    hold off
    ylabel('Mean fluorescence nucleus')
    xlabel('Time into nc (min)')
    

    
    
    
    
    
    %% Operations to change the figures
    
    
    figure(FitFigure)
    try
        ct=waitforbuttonpress;
    catch
        error('Fits not saved.');
    end
    cc=get(FitFigure,'currentcharacter');
    cm=get(gca,'CurrentPoint');
    
    %Move between AP positions
    if (ct~=0)&(cc=='.')&(i<41)
        i=i+1;
    elseif (ct~=0)&(cc==',')&(i>1)
        i=i-1;
        
        %Approve, disapprove fit
    elseif (ct~=0)&(cc=='q')
        if FitResults(i,CurrentNC-11).Approved==-1
            FitResults(i,CurrentNC-11).Approved=1;
        elseif FitResults(i,CurrentNC-11).Approved==1
            FitResults(i,CurrentNC-11).Approved=-1;
        end
        
        %Re-make the fit
    elseif (ct~=0)&(cc=='f')
        FitResults(i,CurrentNC-11).X=[0,0];
        FitResults(i,CurrentNC-11).Y=[0,0];
        FitResults(i,CurrentNC-11).Approved=-1;
        %Switch NCs
    elseif (ct~=0)&(cc=='m')&CurrentNC<14
        CurrentNC=CurrentNC+1;
    elseif (ct~=0)&(cc=='n')&CurrentNC>12
        CurrentNC=CurrentNC-1;  
        
        %Save
    elseif (ct~=0)&(cc=='v')
        
    Name = [DropboxFolder,filesep,Prefix,'\MeanFitsV2.mat']
    save(Name,'FitResults')
    display('MeanFitsV2.mat saved')
        
        %Debug mode
    elseif (ct~=0)&(cc=='9')
        keyboard
        
    end
    
end




%Name = strcat('E:\YangJoon\LivemRNA\Data\Dropbox\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\',Prefix,'_InitialFit.mat');
Name = [DropboxFolder,filesep,Prefix,'\MeanFitsV2.mat']
save(Name,'FitResults')
display('MeanFitsV2.mat saved')


close(FitFigure)
end