%% Creating a model to predict deltaI and deltaR 

% Defining parameters for simulation
dispersion = [0.439, 0.4068, -10.996]; % coefficients for: E = E0 + ak^2 + bk^2
abc = kconstants; % physical parameters of copper
a0 = abc.a; % space between CO atoms
nhex = 5; % radius of hex, number of COs from center diagonally
nspec = 451;
E = linspace(-0.4, 0.5, nspec); % energies
sf =10; % scale factor, spacing between each CO (sf*a0) % 0.964 is hack for fixing dispersion

% Defining the geometry of artificial graphene lattice
vp = khex(nhex, sf*a0,1); % position of COs
vspec = [0,sf*a0/sqrt(3); sf*a0/2,0]; % position of measurements (top and bond sides)

% Creating Training
training_size = 2500;
training = zeros(training_size, nspec*2 + 2 ); % initialize matrix

rng('default'); % seed of randomness
delta = rand([training_size,2]);
% delta I should be from 0 to 1
% delta R should be between -pi/2 to 0
delta(:,1) = (delta(:,1)-1)*pi/2; % transfering numbers to correct range for delta R

% Inserting delta values in first two columns of training data
training(:,1:2)=delta;

for i = 1:training_size
    
    d = delta(i,1)+sqrt(-1)*delta(i,2);
    training(i,3:end) = reshape(kspec(vp, vspec, E, d, dispersion),1,2*nspec);  
end

%% Saving the training data
save('/Users/emory/Documents/GitHub/DS_ResearchProject_ND/Training_Data/ES_AG_Spec_data.mat', 'training');
csvwrite('/Users/emory/Documents/GitHub/DS_ResearchProject_ND/Training_Data/ES_AG_Spec_data.csv', training);
%% Point Specs on Hex lattices
% Measurement with 10a spacing. % save for predicting data
load 'Spec10a.mat'
h0=ksmooth(mean(didv0,2),5);
h10b=ksmooth(mean(didvb10,2),5); h10t=ksmooth(mean(didvt10,2),5);
h10br=h10b./h0; h10tr=h10t./h0;
figure; plot(v,[h10br h10tr]);

% sim the hexagonal lattice specs
a = 2.55; nhex = 5; E = linspace(-0.4, 0.5, 451)'; 
sf =10; vp = khex(nhex, sf*a,1); vspec=[0,sf*a/sqrt(3); sf*a/2,0]; 
simh10=zeros(size(E,1),size(vspec,1));
for ni=1:size(vspec,1)
    simh10(:,ni) = kspec(vp, vspec(ni,:), E,(-.15+0.05*sqrt(-1)));
end

%% Create Data for Experimental Predicting 
topSide_pnts = interp1(v/1000, h10tr, E); 
bondSide_pnts = interp1(v/1000, h10br, E); 
all_cols = [topSide_pnts', bondSide_pnts'];
%% Saving the training data
save('/Users/emory/Documents/GitHub/DS_ResearchProject_ND/Training_Data/ES_AG_Exp_data.mat', 'all_cols');
csvwrite('/Users/emory/Documents/GitHub/DS_ResearchProject_ND/Training_Data/ES_AG_Exp_data.csv', all_cols);

%% Plotting 
figure; plot(v,h10tr,E*1000,simh10(:,1));
figure; plot(v,h10br,E*1000,simh10(:,2));
