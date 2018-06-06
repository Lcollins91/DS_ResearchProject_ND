%% Hexagon Spec Data
% Want to make training data that uses a scale factor to effectively alter the
% dispersion simulations

%Also, want to make the training data have the same resolution as the
%experimental data -- any downsampling should happen in python


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
nv2 = 451;
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
bias4 = linspace(-0.4, 0.5, nv2);
expSpec = 0.5*(spec1 + spec3);

%First need spec points above E = -0.4 only
bias_exp2 = bias(bias>=-0.4);
bias_exp3 = linspace(-0.25, 0.25, nv);
expSpec2 = expSpec(bias>= -0.4);

% %Now need to downsample
% expSpec3 = interp1(bias_exp2, expSpec2, bias3);
% expSpec4 = interp1(bias_exp2, expSpec2, bias_exp3);



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
% hold on
% plot(bias_exp3, expSpec4);

% csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonExperimentalData053118_v2.csv', expSpec4);
% csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonExperimentalData053118_specPoints.csv', expSpec3);

%All spec values at energies above -400mV without downsampling
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonExperimentalData060618_v3.csv', expSpec2');

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
scatter(hex_exp2(:,1), hex_exp2(:,2),  'LineWidth', 20)

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

%% Generating specs for a range of deltas, and scale factors

training_size = 12000;
%training_size = 2;
training1 = cell(training_size,2);

rng('default'); 
vars = rand([3,1,training_size]);
%delta I should be from 0 to 1
%delta R should be between -pi/2 to 0
vars(2,1,:) = (vars(2,1,:)-1)*pi/2;

%scale_factor should be between 0.9 and 1.1
vars(3,1,:) = vars(3,1,:)./5+0.9;

%x_sim = linspace(0,140,nspec);
%k_sim = kv2k(E,dispersion1);

% nv = 201;
% bias3 = linspace(-0.4, 0.5, nv);
% bias4 = linspace(-0.25, 0.25, nv);
nv = 451;
bias_sim = linspace(-0.4, 0.5, nv);
dispersion1 = [0.439, 0.4068, -10.996];

vspec = [0,0];
abc = kconstants;
a0 = abc.a;

trainingA = zeros(training_size, (3+nv));

for i = 1:training_size
    
    
    deltaI = vars(1,1,i);
    deltaR = vars(2,1,i);
    
    delta = deltaR+sqrt(-1)*deltaI;
    
    scale_factor = vars(3,1,i);
    
    vpCO_temp = hexagon_v2(scale_factor*a0);
    
    training1{i,1} = [deltaI, deltaR, scale_factor];
    training1{i,2} = kspec(vpCO_temp, vspec, bias_sim, delta, dispersion1);

    trainingA(i,:) = [deltaI deltaR scale_factor training1{i,2}'];

    i

end
    


%% Saving the training data

save('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonTrainingData060518_v4.mat', 'trainingA');
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonTrainingData060518_v4.csv', trainingA);
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonBias_v4.csv', bias_sim);


%% Generate a second set of simulated training data with the predicted scale factor from the first set

training_size = 12000;
%training_size = 2;
training1 = cell(training_size,2);

rng('default'); 
vars = rand([2,1,training_size]);
%delta I should be from 0 to 1
%delta R should be between -pi/2 to 0
vars(2,1,:) = (vars(2,1,:)-1)*pi/2;

% %scale_factor should be between 0.9 and 1.1
% vars(3,1,:) = vars(3,1,:)./5+0.9;

scale_factor = 0.964;

%x_sim = linspace(0,140,nspec);
%k_sim = kv2k(E,dispersion1);

% nv = 201;
% bias3 = linspace(-0.4, 0.5, nv);
% bias4 = linspace(-0.25, 0.25, nv);
nv = 451;
bias_sim = linspace(-0.4, 0.5, nv);
dispersion1 = [0.439, 0.4068, -10.996];

vspec = [0,0];
abc = kconstants;
a0 = abc.a;

trainingA = zeros(training_size, (2+nv));

for i = 1:training_size
    
    
    deltaI = vars(1,1,i);
    deltaR = vars(2,1,i);
    
    delta = deltaR+sqrt(-1)*deltaI;
    
    %scale_factor = vars(3,1,i);
    
    vpCO_temp = hexagon_v2(scale_factor*a0);
    
    training1{i,1} = [deltaI, deltaR];
    training1{i,2} = kspec(vpCO_temp, vspec, bias_sim, delta, dispersion1);

    trainingA(i,:) = [deltaI deltaR training1{i,2}'];

    i

end


csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonTrainingData060618_v5.csv', trainingA);


%% Comparing Simulation with Predicted Phase to Experimental 

vpCO = hexagon_v2(a0);
new_a = 2.42;
vpCO_v2 = hexagon_v2(new_a);
vpCO_v3 = hexagon_v2(2.3);

scale_factor_ML = 0.964;
vpCO_ML = hexagon_v2(scale_factor_ML*a0);
dispersion1 = [0.439, 0.4068, -10.996];


ML_phase_sf = -0.105 + sqrt(-1)*0.175;

Laura_pred_phase = -0.03 + sqrt(-1)*0.175; 
Emory_pred_phase = -0.15 + sqrt(-1)*0.05;
pred_phase_2 = 0.395*sqrt(-1) -0.104;

sim_pred_spec = kspec(vpCO, vspec, bias3, Laura_pred_phase,dispersion1);
sim_pred_spec1 = kspec(vpCO, vspec, bias3, Emory_pred_phase, dispersion1);
sim_pred_spec2 = kspec(vpCO, vspec, bias_exp3, pred_phase_2, dispersion1);


sim_pred_spec_newA = kspec(vpCO_v2, vspec, bias3, Laura_pred_phase, dispersion1);
sim_pred_spec1_newA = kspec(vpCO_v2, vspec, bias3, Emory_pred_phase, dispersion1);
sim_pred_spec2_newA = kspec(vpCO_v2, vspec, bias_exp3, pred_phase_2, dispersion1); 

sim_pred_spec_newA_v3 = kspec(vpCO_v3, vspec, bias3, Laura_pred_phase, dispersion1);
sim_pred_spec1_newA_v3 = kspec(vpCO_v3, vspec, bias3, Emory_pred_phase, dispersion1);
sim_pred_spec2_newA_v3 = kspec(vpCO_v3, vspec, bias_exp3, pred_phase_2, dispersion1); 

sim_pred_spec_MLsf_MLP = kspec(vpCO_ML, vspec, bias3, ML_phase_sf, dispersion1);
sim_pred_spec_MLsf_LP = kspec(vpCO_ML, vspec, bias3, Laura_pred_phase, dispersion1);
sim_pred_spec_MLsf_EP = kspec(vpCO_ML, vspec, bias3, Emory_pred_phase, dispersion1);

figure; 
subplot(1,2,1)
plot(bias3, expSpec3, 'b','LineWidth', 2);

hold on
plot(bias3, sim_pred_spec1, '.m', 'LineWidth', 2);

plot(bias3, sim_pred_spec,'r', 'LineWidth', 2);
plot(bias_exp3, sim_pred_spec2, 'y', 'LineWidth', 2)

subplot(1,2,2)
plot(bias3, expSpec3, 'b','LineWidth', 2);

hold on
plot(bias3, sim_pred_spec1_newA, '.m', 'LineWidth', 2);

plot(bias3, sim_pred_spec_newA,'r', 'LineWidth', 2);
plot(bias_exp3, sim_pred_spec2_newA, 'y', 'LineWidth', 2)

diff1 = sum((expSpec3-sim_pred_spec1').^2)
diff2 = sum((expSpec3 - sim_pred_spec').^2)
diff3 = sum((expSpec4 - sim_pred_spec2').^2)


figure; 
subplot(2,3,1)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec1, '.r', 'LineWidth', 2);
legend('Experimental', 'Emory phase, a = 2.54')

subplot(2,3,2)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec1_newA, '.r', 'LineWidth', 2);
legend('Experimental', 'Emory phase, a = 2.42')

subplot(2,3,3)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec1_newA_v3, '.r', 'LineWidth', 2);
legend('Experimental', 'Emory phase, a = 2.3')

subplot(2,3,4)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec, '.r', 'LineWidth', 2);
legend('Experimental', 'Laura phase, a = 2.54')

subplot(2,3,5)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec_newA, '.r', 'LineWidth', 2);
legend('Experimental', 'Laura phase, a = 2.42')

subplot(2,3,6)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec_newA_v3, '.r', 'LineWidth', 2);
legend('Experimental', 'Laura phase, a = 2.3')

figure;
subplot(1,3,1)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec_MLsf_MLP, '.r', 'LineWidth', 2);
legend('Experimental', 'new ML phase, delta = -0.105 + i*0.175, a = 2.456')
axis square

subplot(1,3,2)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec_MLsf_LP, '.r', 'LineWidth', 2);
legend('Experimental', 'Laura phase, a = 2.456')
axis square

subplot(1,3,3)
plot(bias3, expSpec3, 'b','LineWidth', 2);
hold on
plot(bias3, sim_pred_spec_MLsf_EP, '.r', 'LineWidth', 2);
legend('Experimental', 'Emory phase, a = 2.456')
axis square




% subplot(3,2,5)
% plot(bias3, expSpec3, 'b','LineWidth', 2);
% hold on
% plot(bias3, sim_pred_spec2, '.r', 'LineWidth', 2);
% 
% subplot(3,2,6)
% plot(bias3, expSpec3, 'b','LineWidth', 2);
% hold on
% plot(bias3, sim_pred_spec2_newA, '.r', 'LineWidth', 2);















