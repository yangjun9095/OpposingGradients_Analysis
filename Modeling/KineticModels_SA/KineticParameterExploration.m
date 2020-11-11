%initial conditions, state1=1, state2 = state3=0


% k1=0.0001;
% k2=0.0000001;
t = 1000; % a long long time
ConcA = logspace(0,5,200); %activator concentrations
Kd = 7000; %affinity of the activator for its binding site
% d1 = 0.03; %how much the activator accelerates step 1
% d2 = 0.03; %how much the activator accelerates step 2
Occupancies = (ConcA./Kd)./(1+(ConcA./Kd));

KdParams = logspace(0,4.5,100);
NParams = linspace(2,6,5);
k1Params = logspace(-9,3,15);
k2Params = logspace(-9,3,15);
k3Params = logspace(-9,3,15);
d1Params = logspace(-9,3,15);
d2Params = logspace(-9,3,15);
d3Params = logspace(-9,3,15);


% two step kinetic model

TwoStepKineticHalfMaxes = [];
TwoStepKineticSteepnesses = [];
count = 1;
for k1 = k1Params
    for k2 = k2Params
        for d1 = d1Params
            for d2 = d2Params
                D1 = d1.*Occupancies;
                D2 = d2.*Occupancies;
                kineticY=twoIrreversibleSteps(k1,k2,t,D1,D2);
                %find the concentration at the half max (it's actually the 95th percentile)
                [~,maxYPos] = min(abs(kineticY-prctile(kineticY,95)/2));
                HalfMax = ConcA(maxYPos);
                %calculate the derivative at the half max
                firstDer = gradient(kineticY,ConcA);
                steepness = firstDer(maxYPos);
                %plot(HalfMax,steepness,'ko');
                TwoStepKineticHalfMaxes = [TwoStepKineticHalfMaxes HalfMax];
                TwoStepKineticSteepnesses = [TwoStepKineticSteepnesses steepness];
            end
        end
    end
    count = count+1
end

% three step knetic model
ThreeStepKineticHalfMaxes = [];
ThreeStepKineticSteepnesses = [];
count = 1;
for k1 = k1Params
    for k2 = k2Params
        for k3 = k3Params
            for d1 = d1Params
                for d2 = d2Params
                    for d3 = d3Params
                        D1 = d1.*Occupancies;
                        D2 = d2.*Occupancies;
                        D3 = d3.*Occupancies;
                        kineticY=threeIrreversibleSteps(k1,k2,k3,t,D1,D2,D3);
                        %find the concentration at the half max (it's actually the 95th percentile)
                        [~,maxYPos] = min(abs(kineticY-prctile(kineticY,95)/2));
                        HalfMax = ConcA(maxYPos);
                        %calculate the derivative at the half max
                        firstDer = gradient(kineticY,ConcA);
                        steepness = firstDer(maxYPos);
                        %plot(HalfMax,steepness,'ko');
                        ThreeStepKineticHalfMaxes = [ThreeStepKineticHalfMaxes HalfMax];
                        ThreeStepKineticSteepnesses = [ThreeStepKineticSteepnesses steepness];
                    end
                end
            end
        end
        count = count+1
    end   
end


% simple hill function

HillHalfMaxes=[];
HillSteepnesses=[];
for k = KdParams
    HillY = (ConcA./k)./(1+(ConcA./k));
    [~,maxYPos] = min(abs(HillY-prctile(HillY,95)/2));
    HalfMax = ConcA(maxYPos);
    firstDer = gradient(HillY,ConcA);
    steepness = firstDer(maxYPos);
    HillHalfMaxes = [HillHalfMaxes HalfMax];
    HillSteepnesses = [HillSteepnesses steepness];
end

HillNHalfMaxes=[];
HillNSteepnesses=[];


% hill function with cooperativity

for n = NParams
    for k = KdParams
        HillY = (ConcA./k).^n./(1+(ConcA./k).^n);
        [~,maxYPos] = min(abs(HillY-prctile(HillY,95)/2));
        HalfMax = ConcA(maxYPos);
        firstDer = gradient(HillY,ConcA);
        steepness = firstDer(maxYPos);
        HillNHalfMaxes = [HillNHalfMaxes HalfMax];
        HillNSteepnesses = [HillNSteepnesses steepness];
    end
end


% show results
figure
hold on
scatter(TwoStepKineticHalfMaxes,TwoStepKineticSteepnesses,'r')
scatter(ThreeStepKineticHalfMaxes,ThreeStepKineticSteepnesses,'k')
plot(HillNHalfMaxes,HillNSteepnesses,'g-')
plot(HillHalfMaxes,HillSteepnesses,'b-')
ylabel('steepness at half max');xlabel('half max [A]')
set(gca,'XScale','log');set(gca,'YScale','log')

