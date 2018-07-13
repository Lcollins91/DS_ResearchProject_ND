%% Hexagon Spec Data
% Want to make training data that uses a scale factor to effectively alter the
% dispersion simulations

% In this training data we will save the peak information for the four
% largest peaks in each simulation.


fn1 = '/Users/lauracollins/Documents/Data/2016-08-01 Hexaxgonal Corral/';

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
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/HexagonExperimentalData061318_v4.csv', expSpec2');

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
hex_exp = csvread('/Users/lauracollins/Desktop/DS_ResearchProject_ND/Training_Data/Hexagon/NewHexagonNGEFPoints.csv', 1,1);

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
delta2 = -0.15 + 0.05*sqrt(-1);
delta_emory = -0.109421667 + 0.0604 * sqrt(-1);
delta_predicted = -0.11675 + sqrt(-1)*0.045095;

dispersion1 = [0.439, 0.4068, -10.996];
sf = 0.945;
dispersion2 = [0.439, 0.4068*(sf^2), -10.996*(sf^4)];


test1 = kspec(vpCO, [0,0], bias3, delta, dispersion1);
test2 = kspec(vpCO, [0,0], bias3, delta2, dispersion1);
test3 = kspec(vpCO, [0,0], bias3, delta, dispersion2);
test4 = kspec(vpCO, [0,0], bias3, delta2, dispersion2);

test_predicted = kspec(vpCO, [0,0], bias2, delta_predicted, dispersion2);
test_pred_2 = kspec(vpCO, [0,0], bias2, delta2, dispersion2);
test_emory = kspec(vpCO, [0,0], bias2, delta_emory, dispersion2);

figure; 
% plot(bias3, test1);
hold on
plot(bias3, test2);
% plot(bias3, test3);
plot(bias3, test4);
plot(bias, spec1);

legend('Unscaled Dispersion', 'Scaled Dispersion', 'Experiment')

figure; 
plot(bias_exp2, expSpec2);
hold on
plot(bias2, test_predicted);
plot(bias2, test_pred_2);
plot(bias2, test_emory);
legend('Experimental', 'Predicted phase and scaled dispersion', 'Previous Predicted Phase and scaled dispersion', 'Emory phase')


r2_pred_1 = sum((test_predicted-expSpec2).^2)
r2_pred_2 = sum((test_pred_2 - expSpec2).^2)
r2_pred_emory = sum((test_emory - expSpec2).^2)

rmse_pred_1 = sqrt(mean((test_predicted-expSpec2).^2))
rmse_pred_2 = sqrt(mean((test_pred_2 - expSpec2).^2))
rmse_pred_emory = sqrt(mean((test_emory - expSpec2).^2))
% Finding the four largest peaks in each. 

% figure;
% findpeaks(test1,bias3,'sortstr', 'descend', 'npeaks', 4)
% hold on
% findpeaks(expSpecLimited,bias2,'sortstr', 'descend', 'npeaks', 4);
% findpeaks(test2, bias3, 'sortstr', 'descend', 'npeaks', 4);
%% Simulating training data with spec points, and training data with peak info

sf = 0.945;
dispersion2 = [0.439, 0.4068*(sf^2), -10.996*(sf^4)];

nv = 451;
bias_sim = linspace(-0.4, 0.5, nv);

vspec = [0,0];
abc = kconstants;
a0 = abc.a;

training_size = 12000;

rng('default'); 
vars = rand([2,1,training_size]);
%delta I should be from 0 to 1
%delta R should be between -pi/2 to 0
vars(2,1,:) = (vars(2,1,:)-1)*pi/2;
np = 4;
trainingA = zeros(training_size, (2+nv));
trainingB = zeros(training_size, (2+4*np));

f = waitbar(0,'Simulating Data');

for i = 1:training_size
    
    
    deltaI = vars(1,1,i);
    deltaR = vars(2,1,i);
    
    delta = deltaR+sqrt(-1)*deltaI;
    
    %scale_factor = vars(3,1,i);
    
    vpCO_temp = hexagon_v2(a0);
    
    training1{i,1} = [deltaI, deltaR];
    training1{i,2} = kspec(vpCO_temp, vspec, bias_sim, delta, dispersion2);

    [pks, locs, w, proms] = findpeaks(training1{i,2}, bias_sim, 'sortstr', 'descend', 'npeaks', np);
    
    testA = [locs', pks, locs', w', proms];
    
    testA = sortrows(testA);
    
    trainingA(i,:) = [deltaI deltaR training1{i,2}'];
    trainingB(i,:) = [deltaI deltaR testA(:,2)', testA(:,3)', testA(:,4)', testA(:,5)'];
    
    waitbar(i/training_size, f)
    

end

close(f)


csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/Training_Data/Hexagon/HexagonTrainingData062218_v8_specPoints.csv', trainingA);
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/Training_Data/Hexagon/HexagonTrainingData062218_v8_peakinfo.csv', trainingB);



%% Getting Experimental spec into same format


fn1 = '/Users/lauracollins/Documents/Data/2016-08-01 Hexaxgonal Corral/';

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


[pks, locs, w, proms] = findpeaks(expSpecLimited, bias_sim, 'sortstr', 'descend', 'npeaks', np);
testA = [locs', pks, locs', w', proms];
testA = sortrows(testA); 

expSpec5 = [testA(:,2)', testA(:,3)', testA(:,4)', testA(:,5)'];


expSpec2 = expSpec2';

csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/Training_Data/Hexagon/HexagonExperimentalData062218_v8_specPoints.csv', expSpec2);
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/Training_Data/Hexagon/HexagonExperimentalData062218_v8_peakinfo.csv', expSpec5);


%% Simulating in a small window around average predicted phase

emory_phase = -0.109421667 + 0.0604 * sqrt(-1);
laura_phase = -0.11675556859948584 + 0.04509557223105219 * sqrt(-1);

avg_phase = (emory_phase + laura_phase)/2;

eps = 0.03; 

di = imag(avg_phase);
dr = real(avg_phase);

test_size = 10000;
test_phases = rand(test_size,2);
test_phases = test_phases * (2*eps);
test_phases(:,1) = test_phases(:,1) + dr - eps;
test_phases(:,2) = test_phases(:,2) + di - eps;

test_rmse = zeros(test_size, 1);
test_R2 = zeros(test_size, 1);

f = waitbar(0,'Simulating data');
for i = 1:test_size
    
    test_delta = test_phases(i,1) + test_phases(i,2)*sqrt(-1);
    
    testspec = kspec(vpCO, [0,0], bias2, test_delta, dispersion2);
    
    test_rmse(i) = sqrt(mean((testspec - expSpec2).^2));
    test_R2(i) = sum((testspec - expSpec2).^2);
    
    waitbar(i/test_size,f);
    
end

%% Simulating in a slightly larger window around average predicted phase

emory_phase = -0.109421667 + 0.0604 * sqrt(-1);
laura_phase = -0.11675556859948584 + 0.04509557223105219 * sqrt(-1);

avg_phase = (emory_phase + laura_phase)/2;

eps = 0.03;

di_min = 0.01;
di_max = 0.1;
dr_min = -0.2;
dr_max = -0.01;

dr_range = dr_max - dr_min;
di_range = di_max - di_min;

dr_avg = dr_min + dr_range/2;
di_avg = di_min + di_range/2;

di = imag(avg_phase);
dr = real(avg_phase);

test_size = 10000;
test_phases = rand(test_size,2);

test_phases(:,1) = test_phases(:,1)*dr_range + dr_min;
test_phases(:,2) = test_phases(:,2)*di_range + di_min;

n_delta = 100;
test_dr = linspace(dr_min, dr_max, n_delta);
test_di = linspace(di_min, di_max, n_delta);

[DR, DI] = meshgrid(test_dr, test_di);

test_size = n_delta*n_delta;

% test_phases = test_phases * (2*eps);
% test_phases(:,1) = test_phases(:,1) + dr - eps;
% test_phases(:,2) = test_phases(:,2) + di - eps;

test_rmse = zeros(n_delta, n_delta);
test_rmse_2 = zeros(n_delta, n_delta);
test_R2 = zeros(n_delta, n_delta);
test_R2_2 = zeros(n_delta, n_delta);

skip = 46;

f = waitbar(0,'Simulating data');
k = 0;
for i = 1:n_delta
    for j = 1:n_delta
        
        test_delta = DR(i,j) + DI(i,j)*sqrt(-1);

        testspec = kspec(vpCO, [0,0], bias2, test_delta, dispersion2);

        test_rmse(i,j) = sqrt(mean((testspec - expSpec2').^2));
        test_R2(i,j) = sum((testspec - expSpec2').^2);

        %skipping the first peak or so to see what comparison scores that gives
        test_rmse_2(i,j) = sqrt(mean((testspec(skip:end) - expSpec2(skip:end)').^2));
        test_R2_2(i,j) = sum((testspec(skip:end)-expSpec2(skip:end)').^2);
        k = k + 1;
        waitbar(k/test_size,f);
    end
end
%% Image maps

[min_val1,idx1] = min(test_rmse(:));
[row1,col1] = ind2sub(size(test_rmse),idx1);

[min_val2, idx2] = min(test_rmse_2(:));
[row2, col2] = ind2sub(size(test_rmse_2),idx2);

pred_phase_1 = test_dr(col1)+test_di(row1)*sqrt(-1);
pred_phase_2 = test_dr(col2) + test_di(row2)*sqrt(-1);

figure; 
subplot(1,2,1)
imagesc(test_dr, test_di, test_rmse)
hold on
scatter(test_dr(col1), test_di(row1))
kcm('gold')

subplot(1,2,2)
imagesc(test_dr, test_di, test_rmse_2)
hold on
scatter(test_dr(col2), test_di(row2))
kcm('gold')


%%

close all
[y1, ind1] = min(test_rmse);
[y1_1, ind1_1] = min(test_rmse_2);
[y2, ind2] = min(test_R2);
[y2_1, ind2_1] = min(test_R2_2);

figure; 
subplot(2,2,1)
plot(test_phases(:,1), test_rmse,'.',  'LineStyle', 'None');
hold on
scatter(test_phases(ind1,1), test_rmse(ind1))

subplot(2,2,2)
plot(test_phases(:,2), test_rmse,'.',  'LineStyle', 'None')
hold on
scatter(test_phases(ind1,2), test_rmse(ind1))

subplot(2,2,3)
plot(test_phases(:,1), test_rmse_2,'.',  'LineStyle', 'None');
hold on
scatter(test_phases(ind1_1,1), test_rmse_2(ind1_1))

subplot(2,2,4)
plot(test_phases(:,2), test_rmse_2,'.',  'LineStyle', 'None');
hold on
scatter(test_phases(ind1_1,2), test_rmse_2(ind1_1))

figure;
subplot(2,2,1)
plot(test_phases(:,1), test_R2, '.', 'LineStyle', 'None')
hold on
scatter(test_phases(ind2,1), test_R2(ind2))

subplot(2,2,2)
plot(test_phases(:,2), test_R2, '.', 'LineStyle', 'None')
hold on
scatter(test_phases(ind2,2), test_R2(ind2))

subplot(2,2,3)
plot(test_phases(:,1), test_R2_2,'.',  'LineStyle', 'None');
hold on
scatter(test_phases(ind2_1,1), test_R2_2(ind2_1))

subplot(2,2,4)
plot(test_phases(:,2), test_R2_2,'.',  'LineStyle', 'None');
hold on
scatter(test_phases(ind2_1,2), test_R2_2(ind2_1))


figure; 

pred_delta = test_phases(ind1,1) + test_phases(ind1,2)*sqrt(-1);
pred_spec = kspec(vpCO, [0,0], bias2, pred_delta, dispersion2);
pred_delta_2 = test_phases(ind1_1, 1) + test_phases(ind1_1,2)*sqrt(-1);
pred_spec_2 = kspec(vpCO, [0,0], bias2, pred_delta_2, dispersion2);


plot(bias2, pred_spec)
hold on
plot(bias_exp2, expSpec2)
plot(bias2, pred_spec_2)


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















