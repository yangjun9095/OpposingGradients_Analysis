function OpposingGradient_Modeling_StaticInputs
%% Start from a simple case of six activator sites and one repressor site
% Here, what I calculate is the fold-change (either of rate or mRNA)
% Define the parameters that I need

% A (A/K_A) : Activator (Bcd)
% A can be defined as A^6, to take into account for the cooperativity, like
% Hill model.
% R (R/K_R) : Repressor (Runt)
% Here, Bcd and Runt has AP dependence (gradient), so now I will just take
% one timepoint, and plug those into the model to get some predictions
% (instantaneous rate?)
% for Bcd and Runt, let's start with somewhat ideal shape of gradient,
% expoenential, and made up gradient for Runt. I need to change these more
% realistically later, and also scaling!

%% 1) wR (interaction between Activator and Repressor in terms of binding)
clear all
X = linspace(0,1,41); 
BcdScale = 10;
RuntScale = 10;

% for Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.
A = BcdScale * exp(-5*X); 
A6 =  A.^6;

R = RuntScale * (-1* exp((X-0.5).^2) + 1.3)/0.3;


% NA : number of activator binding sites, 1
% NR : number of repressor binding sites, 1
% wR : interaction between activator and repressor, 0 for competition, and
% 1 for no interaction in terms of the binding.
wR = [0 0.5 1];
% This distinction needs to be revisited since we don't actually kwnow how
% Runt represses Txn.

% r : rate when activator is binding
% rR : RNAP loading rate decrease, in case of repressor binding 
r = 10;
rR = 5;

% Partition function : 
% Rate_1 : Rate of RNAP loading for 1 repressor binding site
%Rate_1 = zeros(41,length(wR));
for i=1:length(wR)
    Rate_1(:,:,i) = (r*A6 + (r-rR)*A6.*R*wR(i)) ./ (1 + A6 + R + A6.*R*wR(i)) ; 
end
% Rate_0 : Rate of RNAP loading for 1 repressor binding site
Rate_0 = A6*r ./ (1+A6) ;

% I'm calculating this since this is now unitless, so easy to compare between data and prediction.
Ratio = Rate_1./Rate_0; 

%% plot Rate_1 and Rate_0 for wR
hold on
    plot(X,A)
    plot(X,R*0.5)
    plot(X,Rate_0)
    plot(X,Rate_1(:,:,1))
    plot(X,Rate_1(:,:,2))
    plot(X,Rate_1(:,:,3))
% for i=1:length(wR)
%     plot(X,Ratio(:,:,i))
% end
title('RNAP loading rate over AP - w_R')
xlabel('AP')
ylabel('RNAP loading Rate (AU)')
legend('Activator','Repressor','Rate_0','Rate_1-wR=0','Rate_1-wR=0.5','Rate_1-wR=1')
%,num2str(wR(1)),num2str(wR(2)),num2str(wR(3)))
standardizeFigure_YJK(gca,legend,[])

%% 2) r/rR ratio
clear all

X = linspace(0,1,41); 
BcdScale = 10;
RuntScale = 10;

% for Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.
A =  BcdScale *exp(-5*X); 
A6 =  A.^6;

R = RuntScale * (-1* exp((X-0.5).^2) + 1.3)/0.3;


% NA : number of activator binding sites, 1
% NR : number of repressor binding sites, 1
% wR : interaction between activator and repressor, 0 for competition, and
% 1 for no interaction in terms of the binding.
wR = 1;%[0 0.5 1];
% This distinction needs to be revisited since we don't actually kwnow how
% Runt represses Txn.

% r : rate when activator is binding
% rR : RNAP loading rate decrease, in case of repressor binding 
r = 10;
%rR = 5;
% Think about r/rR
%rRatio = logspace(0,2,10);
rRatio = [1 2 5 10 100 1000];
rR = r./rRatio;
rRatio = 1./rRatio;
% Partition function : 
% Rate_1 : Rate of RNAP loading for 1 repressor binding site
%Rate_1 = zeros(41,length(wR));
for i=1:length(rRatio)
    Rate_1(:,:,i) = (r*A6 + (r-rR(i)).*A6.*R*wR) ./ (1 + A6 + R + A6.*R*wR) ; 
end
% Rate_0 : Rate of RNAP loading for 1 repressor binding site
Rate_0 = A6*r ./ (1+A6) ;

% I'm calculating this since this is now unitless, so easy to compare between data and prediction.
Ratio = Rate_1./Rate_0; 

%% Plot Rate_1 for different r/rR
hold on
% plot(X,A)
% plot(X,R)
plot(X,Rate_0)
for i=1:length(rR)
    %plot(0:0.025:1,Rate_1(:,:,i)./Rate_0(:,:))
    plot(0:0.025:1,Rate_1(:,:,i))
    pause
end
title('Rate_1 for different r_{R} /r')
xlabel('AP')
ylabel('RNAP loading rate (AU)')
legend('Rate_0',num2str(rRatio(1)),num2str(rRatio(2)),num2str(rRatio(3)),...
        num2str(rRatio(4)),num2str(rRatio(5)),num2str(rRatio(6)))
standardizeFigure_YJK(gca,legend,[])

%% 3) Scaling of R/K_R (basically K_R)
clear all
X = linspace(0,1,41); 
BcdScale = 10;
RuntScale = [10^(-2) 10^(-1) 1 10 10^2 10^3]; %10^4 10^5 10^6 10^7 10^8 10^9 10^10];

% for Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.
A = BcdScale *exp(-5*X);
A6 =  A.^6;

for i=1:length(RuntScale)
    R(i,:) = RuntScale(i) * (-1* exp((X-0.5).^2) + 1.3)/0.3;
end

% NA : number of activator binding sites, 1
% NR : number of repressor binding sites, 1
% wR : interaction between activator and repressor, 0 for competition, and
% 1 for no interaction in terms of the binding.
wR = 1; %[0 0.5 1];
% This distinction needs to be revisited since we don't actually kwnow how
% Runt represses Txn.

% r : rate when activator is binding
% rR : RNAP loading rate decrease, in case of repressor binding 
r = 10;
rR = 5;

% Partition function : 
% Rate_1 : Rate of RNAP loading for 1 repressor binding site
%Rate_1 = zeros(41,length(wR));
% Rate_0 : Rate of RNAP loading for 1 repressor binding site
Rate_0 = A6*r ./ (1 + A6) ;

for i=1:length(RuntScale)
    Rate_1(i,:) = (r*A6 + (r-rR)*A6.*R(i,:)*wR) ./ (1 + A6 + R(i,:) + A6.*R(i,:)*wR) ; 
end

% I'm calculating this since this is now unitless, so easy to compare between data and prediction.
Ratio = Rate_1./Rate_0; 

%% Plot Rate_1 and Rate_0 for different K_R ratio
hold on
plot(X,A)
plot(X,R(4,:))
plot(X,Rate_0)
for i=1:length(RuntScale)
    plot(X,Rate_1(i,:))
    pause
end

title('RNAP loading rate over AP - K_R')
xlabel('AP')
ylabel('RNAP loading rate(AU)')
legend('Activator','Repressor','Rate_0','Rate_1-R=10^{-2}','Rate_1-R=10^{-1}',...
       'Rate_1-R=1','Rate_1-R=10','Rate_1-R=10^2','Rate_1-R=10^3')%,...
       %'Rate_1-R=10^4','Rate_1-R=10^5')
standardizeFigure_YJK(gca,legend,[])

%% 4) Scaling of A/K_A (basically K_A)
clear all
X = linspace(0,1,41); 
BcdScale = [10^(-3) 10^(-2) 10^(-1) 1 10 10^2 10^3 10^4 10^5 10^6];
RuntScale = 10;

% for Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.
for i=1:length(BcdScale)
    A = exp(-5*X);
    A6(i,:) =  BcdScale(i) * A.^6;
end

R= RuntScale * (-1* exp((X-0.5).^2) + 1.3)/0.3;

% NA : number of activator binding sites, 6
% NR : number of repressor binding sites, 1
% wR : interaction between activator and repressor, 0 for competition, and
% 1 for no interaction in terms of the binding.
wR = 1; %[0 0.5 1];
% This distinction needs to be revisited since we don't actually kwnow how
% Runt represses Txn.

% r : rate when activator is binding
% rR : RNAP loading rate decrease, in case of repressor binding 
r = 10;
rR = 5;

% Partition function : 
% Rate_1 : Rate of RNAP loading for 1 repressor binding site
%Rate_1 = zeros(41,length(wR));
% Rate_0 : Rate of RNAP loading for 1 repressor binding site
Rate_0 = A6*r ./ (1 + A6) ;

for i=1:length(BcdScale)
    Rate_1(i,:) = (r*A6(i,:) + (r-rR)*A6(i,:).*R*wR) ./ (1 + A6(i,:) + R + A6(i,:).*R*wR) ; 
end

% I'm calculating this since this is now unitless, so easy to compare between data and prediction.
Ratio = Rate_1./Rate_0; 

%% Plot Rate_1 and Rate_0 for different K_A ratio
hold on
plot(X,A)
plot(X,R*0.5)
plot(X,Rate_0(10,:))
for i=1:length(BcdScale)
    plot(X,Rate_1(i,:))
    pause
end

title('RNAP loading rate over AP - K_R')
xlabel('AP')
ylabel('RNAP loading rate(AU)')
legend('Activator','Repressor','Rate_0','Rate_1-A=10^{-3}','Rate_1-R=10^{-2}',...
       'Rate_1-R=10^{-1}','Rate_1-R=1','Rate_1-R=10','Rate_1-R=10^2',...
       'Rate_1-R=10^3','Rate_1-R=10^4','Rate_1-R=10^5','Rate_1-R=10^6')%,...
       %'Rate_1-R=10^4','Rate_1-R=10^5')
standardizeFigure_YJK(gca,legend,[])

%% r2 prediction for simple vs additive repression
clear all
wR = 1;
r = 10;
rR = 5;

X = linspace(0,1,41); 
BcdScale = 10;
RuntScale = 10;

% for Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.
A = BcdScale *exp(-5*X);
A6 =  A.^6;

R = RuntScale * (-1* exp((X-0.5).^2) + 1.3)/0.3;

% Partition function : 
% Rate_1 : Rate of RNAP loading for 1 repressor binding site
%Rate_1 = zeros(41,length(wR));

Rate_1 = (r*A6 + (r-rR).*A6.*R*wR) ./ (1 + A6 + R + A6.*R*wR) ; 
Rate_2_simple = (r*A6 + (r-rR).*A6.*R*wR + (r-rR).*A6.*R.^2*wR) ./ (1 + A6 + R + R.^2 + A6.*R*wR +  A6.*R.^2*wR) ; 
Rate_2_additive = (r*A6 + (r-rR).*A6.*R*wR + (r-2*rR).*A6.*R.^2*wR) ./ (1 + A6 + R + R.^2 + A6.*R*wR +  A6.*R.^2*wR) ; 


% Rate_0 : Rate of RNAP loading for 1 repressor binding site
Rate_0 = A6*r ./ (1+A6) ;

% I'm calculating this since this is now unitless, so easy to compare between data and prediction.
% This is actually a fold-change, so we can compare this value with the
% data.
Ratio_01 = Rate_1./Rate_0; 
Ratio_02_simple = Rate_2_simple./Rate_0;
Ratio_02_additive = Rate_2_additive./Rate_0;


%% Plot simple vs additive
hold on
plot(X,Rate_0)
pause
plot(X,Rate_1)
pause
plot(X,Rate_2_simple)
pause
plot(X,Rate_2_additive)

title('RNAP loading rate')
xlabel('AP')
ylabel('RNAP loading rate (AU)')
legend('Rate_0','Rate_1','Rate_2_ simple','Rate_2_ additive')
standardizeFigure_YJK(gca,legend,[])

%% r3 : simple vs additive
clear all

wR = 1;
r = 10;
rR = 5;

X = linspace(0,1,41); 
BcdScale = 10;
RuntScale = 10;

% for Bcd, I start with this, with the length constant of 0.2, which is the value from Feng Liu's paper.
A = BcdScale *exp(-5*X);
A6 =  A.^6;

R = RuntScale * (-1* exp((X-0.5).^2) + 1.3)/0.3;

% Partition function : 
Rate_1 = (r*A6 + (r-rR).*A6.*R*wR) ./ (1 + A6 + R + A6.*R*wR) ; 
Rate_2_simple = (r*A6 + (r-rR).*A6.*R*wR + (r-rR).*A6.*R.^2*wR^2) ./...
                    (1 + A6 + R + R.^2 + A6.*R*wR +  A6.*R.^2*wR^2) ; 
Rate_2_additive = (r*A6 + (r-rR).*A6.*R*wR + (r-2*rR).*A6.*R.^2*wR^2) ./...
                    (1 + A6 + R + R.^2 + A6.*R*wR +  A6.*R.^2*wR^2) ; 

% Rate_3 : Rate of RNAP loading for 3 repressor binding sites
Rate_3_simple = (r*A6 + (r-rR).*A6.*R*wR + (r-rR).*A6.*R.^2*wR + (r-rR).*A6.*R.^3*wR) ./...
                    (1 + A6 + R + R.^2 + R.^3 + A6.*R*wR +  A6.*R.^2*wR^2 +  A6.*R.^3*wR^3) ; 
Rate_3_additive = (r*A6 + (r-rR).*A6.*R*wR + (r-2*rR).*A6.*R.^2*wR + (r-3*rR).*A6.*R.^3*wR) ./...
                    (1 + A6 + R + R.^2 + R.^3 + A6.*R*wR +  A6.*R.^2*wR^2 +  A6.*R.^3*wR^3); 


% Rate_0 : Rate of RNAP loading for 1 repressor binding site
Rate_0 = A6*r ./ (1+A6) ;

% I'm calculating this since this is now unitless, so easy to compare between data and prediction.
% This is actually a fold-change, so we can compare this value with the
% data.
Ratio_03_simple = Rate_3_simple./Rate_0;
Ratio_03_additive = Rate_3_additive./Rate_0;

%% Plot for r3
hold on
plot(X,Rate_0)
pause
plot(X,Rate_1)
pause
plot(X,Rate_2_simple)
pause
plot(X,Rate_2_additive)
pause
plot(X,Rate_3_simple)
pause
plot(X,Rate_3_additive)

title('RNAP loading rate')
xlabel('AP')
ylabel('RNAP loading rate (AU)')
legend('Rate_0','Rate_1','Rate_2_ simple','Rate_2_ additive',...
        'Rate_3_ simple','Rate_3_ additive')
standardizeFigure_YJK(gca,legend,[])
     
%% Calculate the fold-change for r1/r0 (from data)to decide which wR to use

%% Think about 3D fitting

end