%% Compare the MCMC inference results from different bounds/conditions
% The idea here is to check how solid/robust our inference protocol is,
% by comparing the inference results (inferred parameters and the resulting
% fit) from repeated MCMC runs with different bounds/conditions.

%% Import the structure of MCMC results and MCMCdata
% MCMCdata (for input TF and output data)
% MCMC_trials = MCMCTemp;

for i=1:length(MCMC_trials)
    % extract the MCMC results
    chain = MCMC_trials(i).chain;
    results = MCMC_trials.results;
    n_steps = length(chain);
    n_burns = 0.5*n_steps;
    
    % extract the model function
    model_MCMC = results.modelfun;
    
    % define the parameters using their means
    mean_Kb = nanmean(chain(n_burns:end,1));
    mean_Kb = nanmean(chain(n_burns:end,2));
    mean_w_b = nanmean(chain(n_burns:end,3));
    mean_w_bp = nanmean(chain(n_burns:end,4));
    mean_w_rp = nanmean(chain(n_burns:end,5));
    mean_p = nanmean(chain(n_burns:end,6));
    mean_R_max = nanmean(chain(n_burns:end,7));
    % define the parameters
    params_inferred =...
            [mean_Kb, mean_Kr, mean_w_b, mean_w_bp, mean_w_rp, mean_p, mean_R_max];
        
    output = model_6A1R_direct_repression_V2(params_inferred, MCMCdata.xdata);
    
    clf
    hold on
    % Runt null
    errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10}, 'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))

    % Runt WT
    errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10}, 'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))
    % Runt Null
    plot(APaxis(APbinRange), output(1:length(APbinRange)),'Color',ColorChoice(4,:))
    % Runt WT
    plot(APaxis(APbinRange), output(1+length(APbinRange):end),'Color',ColorChoice(1,:))
    pause
end


%% Plot

figure
hold on
% Runt null
errorbar(APaxis, compiledData{construct+1+8,9}, compiledData{construct+1+8,10}, 'o','Color',ColorChoice(4,:),'CapSize',0,'MarkerFaceColor',ColorChoice(4,:))

% Runt WT
errorbar(APaxis, compiledData{construct+1,9}, compiledData{construct+1,10}, 'o','Color',ColorChoice(1,:),'CapSize',0,'MarkerFaceColor',ColorChoice(1,:))