% Simulation of data for STRATAA methods paper

%% Incidence of typhoid fever over 2 years in those with risk factor =1 or =0
% R_TF=5 in those with rf=1 vs 0
TF1=[ones(1250,1); zeros(25000-1250,1)];
TF0=[ones(750,1); zeros(75000-750,1)];

%% Probability of seeking healthcare for a fever
pHC1_low=binornd(1,.1,25000,1);
pHC0_low=binornd(1,.2,75000,1); % Average probability of 0.18    now 0.1741

pHC1_high=binornd(1,.5,25000,1);
pHC0_high=binornd(1,.7,75000,1); % Average probability of 0.65   now 0.6499

%% Probability of having blood drawn for culturing
% For low probability, p_B=0.4 and phi_B=0.68
pB1_low=[binornd(1,.68,1250,1); binornd(1,.4,25000-1250,1)];  %0.4070
pB0_low=[binornd(1,.68,750,1); binornd(1,.4,75000-750,1)];

% For high probability, p_B=0.9 and phi_B=0.95
pB1_high=[binornd(1,.95,1250,1); binornd(1,.9,25000-1250,1)];
pB0_high=[binornd(1,.95,750,1); binornd(1,.9,75000-750,1)];  %0.9005

%% Probability of testing positive for TF
% Two scenarios: low prob of antibiotic usage (=0.2) or high prob of prior
% antibiotic usage (=0.8)
abx_low=binornd(1,.2,100000,1);
abx_high=binornd(1,.8,100000,1);

% Blood culture sensitivity for those w/out (=0) and with (=1) prior Abx usage
sens_low=[.6; .4]; % mean sensitivity = 0.56
sens_high=[.75; .5]; % mean sensitivity = 0.55

pPOS1_low=zeros(25000,1);
for i=1:1250
pPOS1_low(i,1)=binornd(1,sens_low(abx_low(i)+1));
end
pPOS0_low=zeros(75000,1);
for i=1:750
pPOS0_low(i,1)=binornd(1,sens_low(abx_low(i+25000)+1));
end

pPOS1_high=zeros(25000,1);
for i=1:1250
pPOS1_high(i,1)=binornd(1,sens_high(abx_high(i)+1));
end
pPOS0_high=zeros(75000,1);
for i=1:750
pPOS0_high(i,1)=binornd(1,sens_high(abx_high(i+25000)+1));   %0.0111
end

%% Define four different scenarios for prob of HC seeking, testing, and prior Abx usage
scenarios={'low','high','low'; 'low','high','high'; 'high','low','low'; 'high','low','high'};

% Output variables for each scenario
riskfactor=[ones(25000,4); zeros(75000,4)];
fever=binornd(1,.1,100000,4); % assuming incidence of fever (other than TF) is 5% per year
fever=max(fever,[TF1*ones(1,4); TF0*ones(1,4)]); % incidence of other fever + TF
soughtcare=fever.*[pHC1_low pHC1_low pHC1_high pHC1_high; pHC0_low pHC0_low pHC0_high pHC0_high]; % =1 if had fever and sought care, =0 otherwise
bctest=soughtcare.*[pB1_high pB1_high pB1_low pB1_low; pB0_high pB0_high pB0_low pB0_low]; 
sum(bctest)
TFpositive=bctest.*[pPOS1_low pPOS1_high pPOS1_low pPOS1_high; pPOS0_low pPOS0_high pPOS0_low pPOS0_high];
sum(TFpositive)

%% Replace values that should be missing with NaN
for i=1:100000
    for j=1:4
        if soughtcare(i,j)==0
            bctest(i,j)=NaN; % only have data on those who actually sought care
            TFpositive(i,j)=NaN; % only have data on those who sought care AND were tested
        end
        if bctest(i,j)==0
            TFpositive(i,j)=NaN;
        end
    end
end


"TFpositive" "bctest"     "fever"      "riskfactor" "soughtcare"
