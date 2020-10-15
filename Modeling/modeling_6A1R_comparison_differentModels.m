%% generate plots to compare different models
% competition vs. quenching vs. direct repression

% 1) Import Bcd and Runt datasets (averaged)
% 2) Import Txn output data (compiledData)
% 3) Import parameters from modeling_6A_1R_XXX scripts

% 4) For each set of parameters from a model, plug into the model_6A1R_MODE script
% to generate the model prediction, then plot altogether with the real
% data.

%%
%% real data
construct = 5; %[001]
Rate = compiledData{construct+1,9};
Rate_SEM = compiledData{construct+1,10};

Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

% calculate the fold-change, FC
FC_001 = Rate./Rate_null;
% calculate the joint error first, using the fractional error
fracError1 = Rate_SEM./Rate;
fracError2 = Rate_null_SEM./Rate_null;
FC_SEM_001 = sqrt(fracError1.^2 + fracError2.^2).*FC_001;

construct = 6; %[001]
Rate = compiledData{construct+1,9};
Rate_SEM = compiledData{construct+1,10};

Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

% calculate the fold-change, FC
FC_010 = Rate./Rate_null;
% calculate the joint error first, using the fractional error
fracError1 = Rate_SEM./Rate;
fracError2 = Rate_null_SEM./Rate_null;
FC_SEM_010 = sqrt(fracError1.^2 + fracError2.^2).*FC_010;

%% Part. Bicoid saturating regime
clear FC_quenching
clear FC_direct
clear FC_competition

P = logspace(-3,-1,10); % p
W_bp = logspace(0.3,2,10); % w_bp
R = logspace(-2,0,10); % r
W_brp = logspace(-3,-1,10); % w_brp
W_rp = logspace(-3,-1,10); % w_rp
W_br = logspace(-3,-1,10); % w_br
W_b = logspace(0.3,2,10); % w_b
b = 1000; % [Runt]/K_R
% b=0.1; % [Bcd]/Kd

for i=1:length(P)
    for j=1:length(W_bp)
        for k=1:length(R)
            for l=1:length(W_brp)
                for m=1:length(W_b)
                    % pull out each value
                    p = P(i);
                    w_bp = W_bp(j);
                    r = R(k);
                    % repression term
                    w_brp = W_brp(l);
                    w_rp = W_rp(l);
                    w_br = W_br(l);

                    w_b = W_b(m);
                    
                    % assemble the parameters for inputs
                    params_quench = [b, r, p, w_b, w_bp, w_brp];
                    params_direct = [b, r, p, w_b, w_bp, w_rp];
                    params_compete = [b, r, p, w_b, w_bp, w_br];
                    % Calculate the FC for each scenario using
                    % FC_mechanism_param_exp.m function
                    FC_quenching(i,j,k,l,m) = FC_quenching_param_exp(params_quench);
                    FC_direct(i,j,k,l,m) = FC_direct_param_exp(params_quench);
                    FC_competition(i,j,k,l,m) = FC_competition_param_exp(params_quench);
                end
            end
        end
    end
end

FC_quenching = reshape(FC_quenching, [length(P)*length(W_bp)*length(B)*length(W_brp)*length(W_b) 1]);
FC_direct = reshape(FC_direct, [length(P)*length(W_bp)*length(B)*length(W_brp)*length(W_b) 1]);
FC_competition = reshape(FC_competition, [length(P)*length(W_bp)*length(B)*length(W_brp)*length(W_b) 1]);


%% plot
hold on
plot(FC_quenching, ones(size(FC_quenching)),'o','Color',ColorChoice(2,:))%,'LineWidth',1.5)
plot(FC_quenching, 2*ones(size(FC_quenching)),'o','Color',ColorChoice(1,:))%,'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'
plot(FC_direct, 3*ones(size(FC_direct)),'o','Color',ColorChoice(4,:))%,'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'

% plot the real data (FC) at 20% on top
APpos = 20; % [%]
APbin = APpos/2.5 + 1;

errorbar(FC_001(APbin),4, FC_SEM_001(APbin),'o','horizontal','Color',ColorChoice(5,:),'LineWidth',1.5)
errorbar(FC_010(APbin),5, FC_SEM_010(APbin),'o','horizontal','Color',ColorChoice(6,:),'LineWidth',1.5)


xlabel('fold-change')
xlim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1 1.2])
ylim([0 7])
yticks([])

legend('competition','quenching','direct','001','010','Location','NorthWest')
box on

StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,'Bcd_sat_regime_FC_different_models','.tif']); 
% saveas(gcf,[FigPath,filesep,'Bcd_sat_regime_FC_different_models','.pdf']); 


%% Part. Repressor saturating regime
clear FC_quenching_post
clear FC_direct_post
clear FC_competition_post
% From the matrices, pull one value for each variable.
P = logspace(-3,-1,10); % p
W_bp = logspace(0,1,10); % w_bp
B = logspace(-3,-2,10); % r
W_brp = logspace(-3,-1,10); % w_brp
W_rp = logspace(-3,-1,10); % w_rp
W_br = logspace(-3,-1,10); % w_br
W_b = logspace(0,1,10); % w_b
r = 100; % [Runt]/K_R
% b=0.1; % [Bcd]/Kd

for i=1:length(P)
    for j=1:length(W_bp)
        for k=1:length(B)
            for l=1:length(W_brp)
                for m=1:length(W_b)
                    % pull out each value
                    p = P(i);
                    w_bp = W_bp(j);
                    b = B(k);
                    % repression term
                    w_brp = W_brp(l);
                    w_rp = W_rp(l);
                    w_br = W_br(l);

                    w_b = W_b(m);
                    
                    % assemble the parameters for inputs
                    params_quench = [b, r, p, w_b, w_bp, w_brp];
                    params_direct = [b, r, p, w_b, w_bp, w_rp];
                    params_compete = [b, r, p, w_b, w_bp, w_br];
                    % Calculate the FC for each scenario using
                    % FC_mechanism_param_exp.m function
                    FC_quenching_post(i,j,k,l,m) = FC_quenching_param_exp(params_quench);
                    FC_direct_post(i,j,k,l,m) = FC_direct_param_exp(params_quench);
                    FC_competition_post(i,j,k,l,m) = FC_competition_param_exp(params_quench);
                end
            end
        end
    end
end

FC_quenching_post = reshape(FC_quenching_post, [length(P)*length(W_bp)*length(B)*length(W_brp)*length(W_b) 1]);
FC_direct_post = reshape(FC_direct_post, [length(P)*length(W_bp)*length(B)*length(W_brp)*length(W_b) 1]);
FC_competition_post = reshape(FC_competition_post, [length(P)*length(W_bp)*length(B)*length(W_brp)*length(W_b) 1]);

%% plot
hold on
plot(FC_competition_post, ones(size(FC_competition_post)),'o','Color',ColorChoice(2,:))%,'LineWidth',1.5)
plot(FC_quenching_post, 2*ones(size(FC_quenching_post)),'o','Color',ColorChoice(1,:))%,'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'
plot(FC_direct_post, 3*ones(size(FC_direct_post)),'o','Color',ColorChoice(4,:))%,'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'

% plot the real data (FC) at 20% on top
APpos = 40; % [%]
APbin = APpos/2.5 + 1;
errorbar(FC_001(APbin),4, FC_SEM_001(APbin),'o','horizontal','Color',ColorChoice(5,:),'LineWidth',1.5)
errorbar(FC_010(APbin),5, FC_SEM_010(APbin),'o','horizontal','Color',ColorChoice(6,:),'LineWidth',1.5)

xlabel('fold-change')
xlim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1 1.2])
ylim([0 7])
yticks([])

legend('competition','quenching','direct','001','010','Location','NorthWest')
box on

StandardFigure(gcf,gca)

% Save the plot
% saveas(gcf,[FigPath,filesep,'Runt_sat_regime_FC_different_models','.tif']); 
% saveas(gcf,[FigPath,filesep,'Runt_sat_regime_FC_different_models','.pdf']); 