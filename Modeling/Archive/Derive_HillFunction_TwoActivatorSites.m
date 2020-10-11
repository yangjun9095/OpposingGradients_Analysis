% Derive_Hill_function_TwoActivatorSites

% This function is a theoretical exploration for an enhancer with two
% activator binding sites.

% Use P_bound_2A function to get P_bound for a given combination of
% parameters.

% Let's see how the P_bound looks like as we tune the parameters.
% Parameters for P_bound_2A(x0)
% p = x0(1);
% w_ap = x0(2);
% Ka= x0(3);  % scaling factor for Activator
% w_aa=x0(4);
a = linspace(10^(-3),10^3, 10^4);

%P = linspace(10^(-2),1,10);
P = [10^(-2) 10^(-1) 1 10 100];
%w_ap = linspace(1, 10^2, 10);
w_ap = [1 10 100];
%Ka = linspace(10^(-1), 10, 10);
Ka = [0.1 1 10];
%w_aa = linspace(1,10^3,10);
w_aa = [1 10 100 1000];

%Pbound = nan(length(P), length(w_ap), length(Ka), length(w_aa),);

% Parameter exploration
for i=1:length(P)
    for j=1:length(w_ap)
        for k=1:length(Ka)
            for l=1:length(w_aa)
                Pbound{i,j,k,l} = P_bound_2A([P(i),w_ap(j),Ka(k),w_aa(l)]);
            end
        end
    end
end

%% Plot to check
% 1) For a limit where a^2*w_aa > a
% w_aa >> 1, weak promoter, so p < 1
hold on
for l=1:length(w_aa)
    plot(a, (Pbound{1, 3, 2, l}))
    set(gca,'XScale','log')
    pause
end

title('P_{bound} for activator cooperativity')
xlabel('Activator (AU)')
ylabel('P_{bound}')
legend('w_{AA} =1 ','w_{AA} =10^1','w_{AA} = 10^2','w_{AA} = 10^3')
StandardFigure(gcf,gca)

% Save the figures
FigPath = 'E:\YangJoon\LivemRNA\Data\Dropbox\Garcia Lab\Figures\Opposing Gradients\Data\Modeling_2A_ParameterExploration';
saveas(gcf,[FigPath,filesep,'P_bound_','p=',num2str(P(1)),'_w_ap=',num2str(w_ap(3)),'_K_a=',num2str(Ka(2)) , '.tif']); 
saveas(gcf,[FigPath,filesep,'P_bound_','p=',num2str(P(1)),'_w_ap=',num2str(w_ap(3)),'_K_a=',num2str(Ka(2)) , '.pdf']); 