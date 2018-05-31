%% Hexagon Spec Data

fn1 = '/Users/lauracollins/Desktop/Data/2016-08-01 Hexaxgonal Corral/';

cd(fn1)

dv5 = 5*10^-3;
mpp = 0.1;

% Load the data
a = kdat('BiasSpec001.dat');
a1 = kdat('BiasSpec002.dat');
a2 = kdat('BiasSpec003.dat');

%Average forward and backward
spec1 = 0.5*(a.Data(:,3)+a.Data(:,5)); 
spec2 = 0.5*(a1.Data(:,3)+a1.Data(:,5));  %only 201 points
spec3 = 0.5*(a2.Data(:,3)+a2.Data(:,5));

%Loading the bare Cu spec
nv = 501;
b1a = k3ds('GridSpec001.3ds');
dV = 5*10^-3;
b1=reshape(b1a.LIX,[100,nv]);
dataCu=mean(b1)'*10^9;

%Ignoring spec2, since it has a different number of points
spec1 = spec1*10^9;
spec3 = spec3*10^9;

%Normalize by bare Cu
spec1 = spec1./dataCu;
spec3 = spec3./dataCu;

bias = a.Data(:,1);
nv = 201;
bias3 = linspace(-0.4, 0.5, nv);
expSpec = 0.5*(spec1 + spec3);

%First need spec points above E = -0.4 only
bias_exp2 = bias(bias>=-0.4);
expSpec2 = expSpec(bias>= -0.4);

%Now need to downsample
expSpec3 = interp1(bias_exp2, expSpec2, bias3);



close all
% figure;
% plot(bias, spec3+1)
% hold on
% plot(bias, spec1)
%  
% figure;
% findpeaks(spec1, 'MinPeakProminence', mpp)
% 
% figure;
% findpeaks(spec3, 'MinPeakProminence', mpp)

%Checking that the original and downsampled specs look the same
figure; 
plot(bias_exp2, expSpec2);
hold on
plot(bias3, expSpec3);

csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonExperimentalData053118_specPoints.csv', expSpec3);



%% Want to compare the experimental specs to simulated
% First peaks are a bit small, and there may be additional peaks in the
% experiment according to Ken

close all

hex1 = hexagon(1);
hex2 = hexagon2(1);
mpp = 0.1;
bias_limit = -0.4;

x = linspace(-90,90,256);


topofile = ksxm('Topo021.sxm');

figure;
plot(x,topofile.Zf(:,128)*10^9)



figure; 
scatter(hex1(:,1), hex1(:,2))
axis image

hold on
% scatter(hex2(:,1), hex2(:,2))

%Load the points from hexagon.ngef -- were parsed using python notebook 
hex_exp = csvread('/Users/lauracollins/Desktop/DS_ResearchProject_ND/NewHexagonNGEFPoints.csv', 1,1);

%This is in nm, convert to Angstroms
hex_exp = hex_exp*10;

scatter(hex_exp(:,1), hex_exp(:,2))

% scan angle from topograph is 128, rotating back hex_exp to see what the
% spacing is. should line up with hex1 then

theta = -15*pi/180; 
% Angle Found by using an arctan function on the change in y, x, between two 
% points

R = [cos(theta), sin(theta); -sin(theta), cos(theta)];

hex_exp2 = hex_exp*R;

scatter(hex_exp2(:,1), hex_exp2(:,2))

% Want to check that we use the same size hexagon as we built

figure;
imagesc(x,x,topofile.Zf)
axis image 
hold on
scatter(hex1(:,1), hex1(:,2))
scatter(hex_exp2(:,1), hex_exp2(:,2))

vpCO = hex_exp2;

bias2 = linspace(-0.5, 0.5, 501);
bias2 = bias2(bias2>= bias_limit);
expSpecLimited = expSpec(bias>=bias_limit);

bias3 = linspace(-0.4, 0.5, 201);

delta = -0.2 + 0.2*sqrt(-1);
delta2 = -0.5 + 0.2*sqrt(-1);
dispersion1 = [0.439, 0.4068, -10.996];


test1 = kspec(vpCO, [0,0], bias3, delta, dispersion1);
test2 = kspec(vpCO, [0,0], bias3, delta2, dispersion1);

figure; 
plot(bias3, test1+1);
hold on
plot(bias, spec1);

figure;
findpeaks(test1,bias3,'MinPeakProminence', 0.04,'NPeaks', 7)
hold on
findpeaks(expSpecLimited,bias2, 'MinPeakProminence', 0.1, 'NPeaks', 7);
findpeaks(test2, bias3, 'MinPeakProminence', 0.04, 'NPeaks', 7);

% It would be difficult to confirm that at different phases the findpeaks function
% is finding the same peaks, so I think I'm going to use Ken's suggestion
% from before where we use a spec with a lower resolution as the features
% instead of information about the peaks. It'll also help to try this out
% since we'd have to do that for graphene anyways

%% Generating specs for a range of deltas

training_size = 3000;
training1 = cell(training_size,2);

rng('default'); 
vars = rand([2,1,training_size]);
%delta I should be from 0 to 1
%delta R should be between -pi/2 to 0
vars(2,1,:) = (vars(2,1,:)-1)*pi/2;
%x_sim = linspace(0,140,nspec);
%k_sim = kv2k(E,dispersion1);

nv = 201;
bias3 = linspace(-0.4, 0.5, nv);
vspec = [0,0];

trainingA = zeros(training_size, (2+nv));

for i = 1:training_size
    
    
    deltaI = vars(1,1,i);
    deltaR = vars(2,1,i);
    delta = deltaR+sqrt(-1)*deltaI;
    
    training1{i,1} = [deltaI, deltaR];
    training1{i,2} = kspec(vpCO, vspec, bias3, delta, dispersion1);

    trainingA(i,:) = [deltaI deltaR training1{i,2}'];

    i

end
    


%% Saving the training data

save('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonTrainingData052818_specPoints.mat', 'trainingA');
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonTrainingData052818_specPoints.csv', trainingA);
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonBias.csv', bias3);




%% Comparing Simulation with Predicted Phase to Experimental 

pred_phase = -0.03 + sqrt(-1)*0.175; 
Emory_pred_phase = -0.15 + sqrt(-1)*0.05;

sim_pred_spec = kspec(vpCO, vspec, bias3, pred_phase, dispersion1);
sim_pred_spec1 = kspec(vpCO, vspec, bias3, Emory_pred_phase, dispersion1);

figure; 
plot(bias3, expSpec3);

hold on
plot(bias3, sim_pred_spec1, '.r');

plot(bias3, sim_pred_spec);

diff1 = sum((expSpec3-sim_pred_spec1').^2);
diff2 = sum((expSpec3 - sim_pred_spec').^2);













