%% make plots on repression from different contributing factors

%% Load the datasets
FC_spatial_mean
FC_spatial_error

FC_spatial_mean_mRNA
FC_spatial_error_mRNA

FC_spatial_mean_fracComp
FC_spatial_error_fracComp


construct = [1,2,5,6,3,7,8,4];
order = 1:8;

FC = FC_spatial_mean.*FC_spatial_error_fracComp;

FC_error = sqrt((FC_spatial_error./FC_spatial_mean).^2 + (FC_spatial_error_fracComp./FC_spatial_mean_fracComp).^2).*FC;

hold on
% errorbar(order, FC_spatial_mean(construct), FC_spatial_error(construct),'o')
% errorbar(order, FC_spatial_mean_fracComp(construct), FC_spatial_error_fracComp(construct),'o')

errorbar(order, FC, FC_error,'o')
errorbar(order, FC_spatial_mean_mRNA(construct), FC_spatial_error_mRNA(construct),'o')
yline(1,'--')
xticklabels({'000','100','001','010','011','110','101','111'})

xlabel('construct')
ylabel('repression')

StandardFigure(gcf,gca)

% saveas(gcf,[FigPath,filesep,'repression_20%-30%_averaged''.tif']); 
% saveas(gcf,[FigPath,filesep,'repression_20%-30%_averaged''.pdf']); 

%% Using AccmulatedData, compiledData

%% Pick an example DataType, check the contribution of each factor to the fold-change

% Let's start with either 001, 011, 101, 110, 111
%% 001
construct = 5;

FC_rate = compiledData{construct+1,9}./compiledData{construct+1+8,9};
fracError1 = compiledData{construct+1,10}./compiledData{construct+1,9};
fracError2 = compiledData{construct+1+8,10}./compiledData{construct+1+8,9};
FC_rate_error = sqrt(fracError1.^2 + fracError2.^2).*FC_rate;

if construct==4 || construct==1
    FC_frac= compiledData{construct+1,16}(:,3)./compiledData{construct+1+8,16};
    fracError1 = compiledData{construct+1,17}(:,3)./compiledData{construct+1,16}(:,3);
    fracError2 = compiledData{construct+1+8,17}./compiledData{construct+1+8,16};
    FC_frac_error = sqrt(fracError1.^2 + fracError2.^2).*FC_frac;
else
    FC_frac= compiledData{construct+1,16}(:,3)./compiledData{construct+1+8,16}(:,3);
    fracError1 = compiledData{construct+1,17}(:,3)./compiledData{construct+1,16}(:,3);
    fracError2 = compiledData{construct+1+8,17}(:,3)./compiledData{construct+1+8,16}(:,3);
    FC_frac_error = sqrt(fracError1.^2 + fracError2.^2).*FC_frac;
end

FC_predict = FC_rate.*FC_frac;
FC_predict_SEM = sqrt((FC_rate_error./FC_rate).^2 + (FC_frac_error./FC_frac).^2).*FC_predict;

FC_mRNA = AccumulatedData{construct+1,8}./AccumulatedData{construct+1+8,8};
fracError1 = AccumulatedData{construct+1,9}./ AccumulatedData{construct+1,8};
fracError2 = AccumulatedData{construct+1+8,9}./ AccumulatedData{construct+1+8,8};
FC_mRNA_SEM = sqrt(fracError1.^2 + fracError2.^2).*FC_mRNA;

hold on
errorbar(APaxis, FC_predict, FC_predict_SEM,'LineWidth',2,'CapSize',0)%,'Color',ColorChoice(construct,:))
errorbar(APaxis, FC_mRNA, FC_mRNA_SEM,'LineWidth',2,'CapSize',0)

ylim([0 1.2])
