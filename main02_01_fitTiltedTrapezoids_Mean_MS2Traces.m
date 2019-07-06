function main02_01_fitTiltedTrapezoids_Mean_MS2Traces
%% Part3. Rates of RNAP loading for trapezium
% added by Yang Joon Kim, 9/1/2018
% The goal is to plot two different RNAP loading rates for trapezium.
% The fitting is done by FitTiltedTrapezoid_Mean.m script

%% Load the datasets
Data(1,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r0_InitialFit.mat');
Data(2,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r1_InitialFit.mat');
Data(3,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r2_InitialFit.mat');
Data(4,1) = load('E:\YangJoon\LivemRNA\Data\Dropbox\OpposingGradient\Data_Processed\FittingMS2Traces\r3_InitialFit.mat');

%% Fill the variables to plot

% T0(Construct, DataSet, AP, NC), same for the others
T0=nan(4,1,41,3);
Rate1=nan(4,1,41,3);
Rate2=nan(4,1,41,3);
Rate3=nan(4,1,41,3);
%SD_rate=nan(4,4,41,3);


%Nb_Approved_Rate=zeros(4,41,3);


for i = 1:4 %Construct r0 to r3 loop
    for j = 1:length(Data(i,:)) % Different number of dataSets for each construct
        for NC= 12:14
            for AP = 1:41
                if Data(i).FitResults(AP,NC-11).Approved == 1
                    
                    T0(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).TimeStart;
                    Rate1(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).RateFit1;
                    Rate2(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).RateFit2;
                    Rate3(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).RateFit3;
                    %SD_rate(i,j,AP,NC-11)=Data(i,j).FitResults(AP,NC-11).SDRateFit;
                    
                    %Nb_Approved_Rate(i,AP,NC-11)=Nb_Approved_Rate(i,AP,NC-11)+1;
                end
            end
        end
    end
end



end