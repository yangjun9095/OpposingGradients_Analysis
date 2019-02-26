function Plot3D_BcdRunt_Output
% Description 
% The goal is to make 3D plots for Bcd, Runt, and output (either MS2 traces
% or loading rate, etc.)
%% Load the datasets
% For the inputs, the data are background-substracted and start from NC12
% (ElapsedTime(nc12) = 0, nc12=1) 

% Inputs
Bicoid = load('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed/Bcd-Averaged');
Runt = load('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\2018-05-24-Runt-JB3-MCP-mCherry-vasa-eGFP1\CompiledNuclei');

% Outputs
Construct(1)=load('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed/r0_From_NC12');
Construct(2)=load('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed/r1_From_NC12');
Construct(3)=load('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed/r2_From_NC12');
Construct(4)=load('E:\Paul-J\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed/r3_From_NC12');

%% Prepare the Datasets (input and output) - Interpolation

% Building two strucures: Input(NC-11), with 3 fields (Time, Bcd, Run)
%                       : Output(NbConstruct+1,NC-11), with two fields
%                                                   (MeanVectorAP,Time)
%
%All the Data are interpolated (10s) and cut so that they are all the same
%size
%
% We also normalize the inputs (so that Kd is a ratio, between 0 and 1)

MaxBcd = max(max(Bicoid.MeanVectorAP));
MaxRun = max(max(Runt.MeanVectorAP));

% Correct the Runt data because the nc14 previously determined was too
% early (no Histone channel)
Runt.nc14=83;
clear c

for NC = 1:3
    
    
    BcdInter=Interpolation_10s(Bicoid,NC,First,Last);
    RunInter=Interpolation_10s(Runt,NC,First,Last);
    
    %take the minimum of the two length
    if BcdInter.L < RunInter.L
        L=BcdInter.L;
    else
        L=RunInter.L;
    end
    
    
    for R=1:4
        DataInter(R)=Interpolation_10s(Construct(R),NC,First,Last);
        
        %Determining the length of DataSets
        if DataInter(R).L < L
            L=DataInter(R).L;
        end
        
    end
    
    Input(NC).Time = BcdInter.Time(1:L);
    Input(NC).Bcd = BcdInter.Data(1:L,:)./MaxBcd;
    Input(NC).Run = RunInter.Data(1:L,:)./MaxRun;
    
    for R=1:4
        Output(R,NC).MeanVectorAP=DataInter(R).Data(1:L,:);
        Output(R,NC).Time=BcdInter.Time(1:L);
        
        
%         figure(R*100+NC)
%         plot(Output(R,NC).Time,Output(R,NC).MeanVectorAP,'.')
    end
    
    OutputInter = Output;
    for R=1:4
        for t=3:length(Input(NC).Time)-2
            for AP=(First:Last)
                Output(R,NC).MeanVectorAP(t,AP)=nansum(OutputInter(R,NC).MeanVectorAP(t-2:t+2,AP))./5;
                nansum(OutputInter(R,NC).MeanVectorAP(t-2:t+2,AP))
              
            end
        end
    end
    
%     figure(NC)
%     plot(Input(NC).Time,Input(NC).Bcd,'.')
%     title('Bcd')
%     
%     figure(10+NC)
%     plot(Input(NC).Time,Input(NC).Run,'.')
%     title('Runt')
    
end


%% Plot the correlation for Dm
%ap = 11;
% Color code for time points
iStart=33; %Start time
iEnd=93; %End time
colormap(jet(256));
cmap=colormap ;
Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

hold on
for i=1:length(Dm)
    for tpoint = 33:15:93
        clear InferredProtein
        clear NBProteinFluo
        clear Z

        InferredProtein = mRNADiffusion(i).InferredProtein;
        NBProteinFluo = mRNADiffusion(i).NBProteinFluo;
        Z = Dm(i)*ones(size(InferredProtein));
        Z = log10(Z);
        scatter3(InferredProtein(tpoint,:),Z(tpoint,:),NBProteinFluo(tpoint,:))
        
        title('Prediction vs Measured protein')
        xlabel('Predicted Protein (AU)')
        ylabel('log(Dm)')
        zlabel('Measured Protein (AU)')
        %zlim([])
        pause
    end
end
%% Correlation for different Dp (protein diffusion coefficient)
Dm = 0;
Dp = [0 0.1 1 10 100];
Tm = 60;
Tp = 50;
rp = 2;
ProteinDiffusion={};

for i=1:length(Dp)
    clear InferredProtein
    clear NBProteinFluo
    [AccumulatedmRNA,Protein,RNAPLoadingRate3,InferredProtein, NBProteinFluo] = ...
                CorrelatemRNAProtein_function(Prefix1,Prefix2,Dm,Dp(i),Tm,Tp,rp)
    ProteinDiffusion(i).InferredProtein = InferredProtein;
    ProteinDiffusion(i).NBProteinFluo = NBProteinFluo;
    ProteinDiffusion(i).Dp = Dp(i);
    
end 

%% Plot the correlation for Dp
%ap = 12;
hold on
for i=1:length(Dp)
    for tpoint = 33:15:93
        clear InferredProtein
        clear NBProteinFluo
        clear Z

        InferredProtein = ProteinDiffusion(i).InferredProtein;
        NBProteinFluo = ProteinDiffusion(i).NBProteinFluo;
        Z = Dp(i)*ones(size(InferredProtein));
        Z = log10(Z);
        scatter3(InferredProtein(tpoint,:),Z(tpoint,:),NBProteinFluo(tpoint,:))
        
        title('Prediction vs Measured protein')
        xlabel('Predicted Protein (AU)')
        ylabel('log(Dp)')
        zlabel('Measured Protein (AU)')
        %zlim([])
        pause
    end
end
%% Correlation for different Tp (protein half-life)
Dm = 0;
Dp = 5;
Tm = 60;
Tp = [1,5,10,50,100,1000];
rp = 2;
ProteinLifeTime={};

for i=1:length(Tp)
    clear InferredProtein
    clear NBProteinFluo
    [AccumulatedmRNA,Protein,RNAPLoadingRate3,InferredProtein, NBProteinFluo] = ...
                CorrelatemRNAProtein_function(Prefix1,Prefix2,Dm,Dp,Tm,Tp(i),rp)
    ProteinLifeTime(i).InferredProtein = InferredProtein;
    ProteinLifeTime(i).NBProteinFluo = NBProteinFluo;
    ProteinLifeTime(i).Tp = Tp(i);
    
end 

%% Plot the correlation for Tp
%ap = 12;
hold on
for i=1:length(Tp)
    for tpoint = 33:15:93
        clear InferredProtein
        clear NBProteinFluo
        clear Z

        InferredProtein = ProteinLifeTime(i).InferredProtein;
        NBProteinFluo = ProteinLifeTime(i).NBProteinFluo;
        Z = Tp(i)*ones(size(InferredProtein));
        Z = log10(Z);
        scatter3(InferredProtein(tpoint,:),Z(tpoint,:),NBProteinFluo(tpoint,:))
        
        title('Prediction vs Measured protein')
        xlabel('Predicted Protein (AU)')
        ylabel('log(Tp)')
        zlabel('Measured Protein (AU)')
        %zlim([])
        pause
    end
end

%% save the result variables
%save(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-03-02-Hb-P2P-MS2V5-MCP-GFP\ComparePredictiontoMeasurement.mat'],...
%    'mRNADiffusion','ProteinDiffusion','ProteinLifeTime')
end