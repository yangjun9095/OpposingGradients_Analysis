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

%[100]
construct = 2; 
Rate = compiledData{construct+1,9};
Rate_SEM = compiledData{construct+1,10};

Rate_null = compiledData{construct+1+8,9};
Rate_null_SEM = compiledData{construct+1+8,10};

% calculate the fold-change, FC
FC_100 = Rate./Rate_null;
% calculate the joint error first, using the fractional error
fracError1 = Rate_SEM./Rate;
fracError2 = Rate_null_SEM./Rate_null;
FC_SEM_100 = sqrt(fracError1.^2 + fracError2.^2).*FC_100;

%[001]
construct = 5; 
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

%[010]
construct = 6; 
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

P = logspace(-3,-2,10); % p
W_bp = logspace(0,1,10); % w_bp
R = logspace(-2,1,10); % r
W_brp = logspace(-3,-1,10); % w_brp
W_rp = logspace(-3,-1,10); % w_rp
W_br = logspace(-3,-1,10); % w_br
W_b = logspace(0,1,10); % w_b
b = 100; % [Runt]/K_R
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
                    % Limits taken analytically
%                     FC_quenching(i,j,k,l,m) = (1 + p*w_bp^6 + r*w_brp^6 + r*p*w_bp^6*w_brp^6)/(1 + p*w_bp^6 + r + r*p*w_bp^6*w_brp^6);
%                     FC_direct(i,j,k,l,m) = (1+r*w_rp+p*w_bp^6+r*p*w_rp*w_bp^6)/(1+r+p*w_bp^6+r*p*w_rp*w_bp^6);
%                     FC_competition(i,j,k,l,m) = 1;
                end
            end
        end
    end
end

FC_quenching = reshape(FC_quenching, [length(P)*length(W_bp)*length(R)*length(W_brp)*length(W_b) 1]);
FC_direct = reshape(FC_direct, [length(P)*length(W_bp)*length(R)*length(W_brp)*length(W_b) 1]);
FC_competition = reshape(FC_competition, [length(P)*length(W_bp)*length(R)*length(W_brp)*length(W_b) 1]);


%% plot the FC (parameter exploration for anterior)
fig_anterior = figure;
hold on
plot(FC_competition, ones(size(FC_competition)),'o','Color',ColorChoice(2,:),'LineWidth',1.5)
plot(FC_quenching, 2*ones(size(FC_quenching)),'Color',ColorChoice(1,:),'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'
plot(FC_direct, 3*ones(size(FC_direct)),'Color',ColorChoice(4,:),'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'

% plot the real data (FC) at 20% on top
APpos = 20; % [%]
APbin = APpos/2.5 + 1;

errorbar(FC_001(APbin),4, FC_SEM_001(APbin),'o','horizontal','Color',ColorChoice(5,:),'LineWidth',1.5)
errorbar(FC_010(APbin),5, FC_SEM_010(APbin),'o','horizontal','Color',ColorChoice(6,:),'LineWidth',1.5)
errorbar(FC_100(APbin),6, FC_SEM_100(APbin),'o','horizontal','Color',ColorChoice(7,:),'LineWidth',1.5)


xlabel('fold-change')
xlim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1 1.2])
ylim([0 7])
yticks([])

legend('competition','quenching','direct','001','010','100','Location','NorthWest')
box on

StandardFigure(fig_anterior,fig_anterior.CurrentAxes)

% Save the plot
fig_anterior.Renderer='Painters';
% saveas(fig_anterior,[FigPath,filesep,'Bcd_sat_regime_FC_different_models','.tif']); 
% saveas(fig_anterior,[FigPath,filesep,'Bcd_sat_regime_FC_different_models','.pdf']); 


%% Part. Repressor saturating regime
clear FC_quenching_post
clear FC_direct_post
clear FC_competition_post
% From the matrices, pull one value for each variable.
P = logspace(-3,-2,10); % p
W_bp = logspace(0.3,1,10); % w_bp
B = logspace(-3,-2,10); % r
W_brp = logspace(-3,-0.5,10); % w_brp
W_rp = logspace(-3,-0.5,10); % w_rp
W_br = logspace(-3,-0.5,10); % w_br
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
clf
window = 1:length(FC_competition_post);
hold on
plot(FC_competition_post(window), ones(size(window)),'Color',ColorChoice(2,:),'LineWidth',1.5)
plot(FC_quenching_post(window), 2*ones(size(window)),'Color',ColorChoice(1,:),'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'
plot(FC_direct_post(window), 3*ones(size(window)),'Color',ColorChoice(4,:),'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'

% plot the real data (FC) at 20% on top
APpos = 40; % [%]
APbin = APpos/2.5 + 1;
errorbar(FC_001(APbin),4, FC_SEM_001(APbin),'o','horizontal','Color',ColorChoice(5,:),'LineWidth',1.5)
errorbar(FC_010(APbin),5, FC_SEM_010(APbin),'o','horizontal','Color',ColorChoice(6,:),'LineWidth',1.5)
errorbar(FC_100(APbin-1),6, FC_SEM_100(APbin-1),'o','horizontal','Color',ColorChoice(7,:),'LineWidth',1.5)

xlabel('fold-change')
xlim([0 1.2])
xticks([0 0.2 0.4 0.6 0.8 1 1.2])
ylim([0 7])
yticks([])

legend('competition','quenching','direct','001','010','100','Location','NorthWest')
box on


StandardFigure(gcf,gca)

% Save the plot
saveas(gcf,[FigPath,filesep,'Runt_sat_regime_FC_different_models','.tif']); 
saveas(gcf,[FigPath,filesep,'Runt_sat_regime_FC_different_models','.pdf']); 

%% generate predictions of FC over AP axis for several sets of parameters
% Question : How do I do this systematically over a large parameter space?

% Quenching
% params = [Kb, Kr, w_b, w_bp, w_brp, p, R_max]; 
% params_quench = [99.9964    0.1000    1.0576    8.7865    0.8447    0.0003  185.5756];
% 
% % FC_fit
% FC_FittedParams = model_FC_6A1R_competition_V1(Bcd, Runt,...
%                         params_quench(1:end-1));    
                    
                    

% Define the range of parameters to pull out, then generate the model
% prediction (to assess the model behavior)
K_B = logspace(-1,2,10); % K_b
K_R = logspace(-1,2,10); % K_r
W_b = logspace(0,1,10); % w_b
W_bp = logspace(0,1,10); % w_bp
W_brp = logspace(-3,-0.5,10); % w_brp
W_rp = logspace(-3,-0.5,10); % w_rp
W_br = logspace(-3,-0.5,10); % w_b     
P = logspace(-3,-1,10); % p

FC_quenching = nan(length(K_B), length(K_R), length(W_b), length(W_bp), length(P), length(W_brp),41);
FC_direct = nan(length(K_B), length(K_R), length(W_b), length(W_bp), length(P), length(W_brp),41);
FC_competition = nan(length(K_B), length(K_R), length(W_b), length(W_bp), length(P), length(W_brp),41);

for i=1:length(K_B)
    for j=1:length(K_R)
        for k=1:length(W_b)
            for l=1:length(W_bp)
                for m=1:length(P)
                    for n=1:length(W_brp) % for the repression terms, w_rp, w_brp, w_br
                    
                        % pull out each value
                        K_b = K_B(i);
                        K_r = K_R(j);
                        w_b = W_b(k);
                        w_bp = W_bp(l);
                        p = P(m);

                        % repression term
                        w_brp = W_brp(n);
                        w_rp = W_rp(n);
                        w_br = W_br(n);

                        % assemble the parameters for inputs
                        params_quench =  [K_b, K_r, w_b, w_bp, w_brp, p];
                        params_direct =  [K_b, K_r, w_b, w_bp, w_rp, p];
                        params_compete = [K_b, K_r, w_b, w_bp, w_br, p];
                        % Calculate the FC for each scenario using
                        % FC_mechanism_param_exp.m function
                        FC_quenching(i,j,k,l,m,n,:) = model_FC_6A1R_quenching_V1(Bcd, Runt, params_quench);
                        FC_direct(i,j,k,l,m,n,:) = model_FC_6A1R_direct_repression_V1(Bcd, Runt, params_direct);
                        FC_competition(i,j,k,l,m,n,:) = model_FC_6A1R_competition_V1(Bcd, Runt, params_compete);

                    end
                end
            end
        end
    end
end

% FC_quenching = reshape(FC_quenching, [length(P)*length(W_bp)*length(R)*length(W_brp)*length(W_b) 1]);
% FC_direct = reshape(FC_direct, [length(P)*length(W_bp)*length(R)*length(W_brp)*length(W_b) 1]);
% FC_competition = reshape(FC_competition, [length(P)*length(W_bp)*length(R)*length(W_brp)*length(W_b) 1]);


%% plot the FC (parameter exploration for anterior)
fig = figure;
hold on
for i=1:length(K_B)
    for j=1:length(K_R)
        for k=1:length(W_b)
            for l=1:length(W_bp)
                for m=1:length(P)
                    for n=1:length(W_brp) % for the repression terms, w_rp, w_brp, w_br

                    end
                end
            end
        end
    end
end

% plot(FC_competition, ones(size(FC_competition)),'o','Color',ColorChoice(2,:),'LineWidth',1.5)
% plot(FC_quenching, 2*ones(size(FC_quenching)),'Color',ColorChoice(1,:),'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'
% plot(FC_direct, 3*ones(size(FC_direct)),'Color',ColorChoice(4,:),'LineWidth',1.5)%,'o','MarkerSize',0.01,'MarkerFaceColor'

% plot the real data (FC) at 20% on top
% APpos = 20; % [%]
% APbin = APpos/2.5 + 1;
% 
% errorbar(FC_001(APbin),4, FC_SEM_001(APbin),'o','horizontal','Color',ColorChoice(5,:),'LineWidth',1.5)
% errorbar(FC_010(APbin),5, FC_SEM_010(APbin),'o','horizontal','Color',ColorChoice(6,:),'LineWidth',1.5)
% errorbar(FC_100(APbin),6, FC_SEM_100(APbin),'o','horizontal','Color',ColorChoice(7,:),'LineWidth',1.5)
% 
% 
% xlabel('fold-change')
% xlim([0 1.2])
% xticks([0 0.2 0.4 0.6 0.8 1 1.2])
% ylim([0 7])
% yticks([])
% 
% legend('competition','quenching','direct','001','010','100','Location','NorthWest')
% box on
% 
% StandardFigure(fig_anterior,fig_anterior.CurrentAxes)
% 
% % Save the plot
% fig_anterior.Renderer='Painters';
% saveas(fig_anterior,[FigPath,filesep,'Bcd_sat_regime_FC_different_models','.tif']); 
% saveas(fig_anterior,[FigPath,filesep,'Bcd_sat_regime_FC_different_models','.pdf']); 