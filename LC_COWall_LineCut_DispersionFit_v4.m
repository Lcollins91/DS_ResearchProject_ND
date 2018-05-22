%% CO Line - fitting to get the dispersion
% Saved 04 - 24 - 2018 
% Summarized in Evernote Note: Undergrads/ Scattering Theory Phase fit/ UPDATE - Dispersion Fit for CO Wall Linecuts


close all

%File path for Laura's laptop
fn1 = '/Users/lauracollins/Desktop/Data/2016-08-29 Copper - Fibonacci/';

%File path for network drive location
%fn1 = 'Z:/gomeslab/Data/2016/2016-08-29 Copper - Fibonacci/';

cd(fn1)

% Load the data
a = k3ds('GridSpec001.3ds');
OLMs = a.LIXf;

% Conductance Maps
nE = size(a.LIXf,3);
Ei=str2double(a.Params.bias_spectroscopy_sweep_start__v_)*10^3;
Ef=str2double(a.Params.bias_spectroscopy_sweep_end__v_)*10^3;
E=linspace(Ei,Ef,nE);

dV=str2double(a.Params.lock_in_amplitude);

maps=cell(nE,1);
for i=1:nE
    maps{i}=ksmooth(kinterp(flipud(a.LIXf(:,:,i))/dV*10^9,1024,256),4);
    maps{i}=maps{i}/mean2(maps{i}(824:1024,:));
%     figure; imagesc(maps{i}); axis image off; kcm('sunrise'); caxis([0.7,1.3]);
end

% 
linecut = zeros(nE, 1024);
linecut2 = zeros(nE, 904);
offset = 0;
amplitude = 1;

figure;
for i = 1:nE
    cut1 = mean(maps{i},2);
    linecut(i,:) = ksmooth(amplitude*(cut1+(i-1)*offset),3);
    linecut2(i,:) = linecut(i, 121:end);
    plot(linecut2(i,:)); 
    hold on;
end

hold off

%% Method 1 - Want to fit a cosine function so we can use the wavelength to determine k at each energy
close all

% should center each linecut around y=0, with the first peak after the wall 
% at x=0. 
linecut3 = zeros(nE, 904);
linecutx = zeros(nE, 904);
linecut4 = cell(nE,1);
for i = 1:nE
    
    linecut3(i,:) = linecut2(i,:) - mean(linecut2(i,:));
    [peaks,locs,~, prom] = findpeaks(linecut2(i,:));
    
    if prom(2)>prom(1) %without this sometimes (4-7) it identified a peak from the wall as the first peak. This makes sure it's a normal sized peak at x=0
        linecutx(i,:) = -locs(2):(903-locs(2));
    else
        linecutx(i,:) = -locs(1):(903-locs(1));
    end
    figure; 
    plot(linecutx(i,:), linecut3(i,:))
    xlim([0,900])
    linecut4{i} = linecut3(i,linecutx(i,:)>=0);
end



%% Method 1 - Fitting a sine wave to the linecut

close all
wavelength = zeros(1,nE);
k_estimated = zeros(1,nE);


for i = 1:nE
    Y = linecut4{i};
    X = (0:(length(Y)-1))*(160/1024); %Angstroms
%     figure; 
%     scatter(X,Y)
%     hold on
%     
    [f,g] = fit(X',Y','sin1');
%     plot(f)
    k_estimated(i) = f.b1; 
   
end
%
E = -0.4:0.1:0.4;
k_estimated = 0.5*k_estimated; %Since it should be fitting to sin(2kr)
k_estimated2 = k_estimated(2:end); %excluding k for E = -0.4
k_original = kv2k(E);

figure;
scatter(k_estimated,E);
hold on
scatter(k_original,E);
ylim([-0.5,0.5])

% Momentum fits from David's Mathematica notebook, rescaled to be in
% angstroms
k_david = (1/10)*[1.24707, 1.62095, 1.94075, 2.22925, 2.50278, 2.75491, 2.96558, 3.19892];

scatter(k_david, E(2:end));

abc = kconstants; 
E0 = -abc.E0;
E2 = E - E0; %explicitly use known E0 as constant by offsetting E - for a few of the fits only

%Explicitly set E0
ft1 = fittype({'x^2', 'x^4'});% here x is k

%Allow E0 to be fit as well
ft2 = fittype({'1', 'x^2', 'x^4'});

%% Method 1 - Case 1
%Here we fit k's from all energies, with E0 fixed
[f_case1,g_case1] = fit(k_estimated', E2', ft1);

%Results for all k's, E0 fixed
alpha_case1 = f_case1.b;
mstar_case1 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case1.a);

%% Method 1 - Case 2
%Here we fit k's from all energies, with E0 free
[f_case2, g_case2] = fit(k_estimated', E', ft2);

%Results for all k's, E0 free
E0_case2 = f_case2.a;
alpha_case2 = f_case2.c;
mstar_case2 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case2.b);
%% Method 1 - Case 3
%Here we fit k's from E>=-0.3, with E0 fixed
[f_case3, g_case3] = fit(k_estimated2', E2(2:end)', ft1);

%Results for k's with E >= -0.3, E0 fixed
alpha_case3 = f_case3.b;
mstar_case3 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case3.a);

%% Method 1 - Case 4
%Here we fit k's from E>= -0.3, with E0 free
[f_case4, g_case4] = fit(k_estimated2', E(2:end)', ft2);

%Results for k's with E >= -0.3, E0 free
E0_case4 = f_case4.a;
alpha_case4 = f_case4.c;
mstar_case4 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case4.b);
%% David's Method - Case 5
%Here we fit David's values for k from the Bessel function fit for E
%>=-0.3, with E0 fixed
[f_case5, g_case5] = fit((k_david)', E2(2:end)', ft1);

%Results for David's k's (E >= -0.3), E0 fixed
alpha_case5 = f_case5.b; 
mstar_case5 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case5.a);

%% David's Method - Case 6
%Here we fit David's values for k from the Bessel function fit for E
%>=-0.3, with E0 free
[f_case6, g_case6] = fit((k_david)', E(2:end)', ft2);

%Results for David's k's (E >= -0.3), E0 free
E0_case6 = f_case6.a;
alpha_case6 = f_case6.c;
mstar_case6 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case6.b);

%% Test Case
%This method does yield 0 for b and m* of 0.38 for the original kv2k dispersion, so
%it is accurate at least in that regard. 
%[f_test,g_test] = fit(k_original', E2', ft1);

% alpha_test = f_test.b; 
% mstar_test = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_test.a);

%% Method 1 & David's Method - Results Summary

%In the end, we use the average of the results from Cases 4 and 6, where
%the slight differences in the momentum fits arise due to the difference in
%the equation used to fit them (sin(2kr) vs. Bessel)

E0_result = (E0_case4 + E0_case6)/2; %E0 = -0.439 +- 0.001
alpha_result = (alpha_case4 + alpha_case6)/2; % alpha = 10.996 +- 0.504
mstar_result = (mstar_case4 + mstar_case6)/2; % mstar = 0.4068 +- 0.001



%% Method 1 - Summary Plot
figure;
scatter(k_estimated, E2);

hold on
scatter(k_original, E2);
plot(f_case1)

scatter(k_david, E2(2:end))
plot(f_case5)


%% Method 2

%Here I will fit the data to a combination of sine and cosine, and set x =
%0 to be the position of the wall
close all 

%Use the topograph to find the position of the wall. 
topo = a.Topo; 
topo = kinterp(ltstrip(flipud(topo),3),1024,256);

topocut = mean(topo,2);
[~,wall_position] = min(topocut);


linecut_wall = linecut(:,wall_position:end);

%x_from_wall = (0:(size(linecut_wall,2)-1))*160/1024; %convert to Angstroms

linecut_wall_fit = zeros(9,size(linecut_wall,2));
linecut_wall_fit_final = cell(9);

% This way, we assume k is the same as that determined in Laura's method
% #1, and set it explicitly for each energy
ft3 = fittype( 'a*sin(2*k*x)/sqrt(k*x) + b*cos(2*k*x)/sqrt(k*x)', 'independent', 'x','problem', 'k');
options1 = fitoptions(ft3);
lower_bounds = [(1/(2*pi))*(-1-sqrt(2))/sqrt(2) -1/(2*pi*sqrt(2))];
upper_bounds = [(1/(2*pi))*(1-sqrt(2))/sqrt(2) 1/(2*pi*sqrt(2))];



ft4 = fittype('fourier1');
options = fitoptions('fourier1');
options.Lower = [0 -Inf -Inf -Inf];
options.Upper = [0 Inf Inf Inf]; %removing the added constant as a fit parameter
options;


ft5 = fittype('(exp(-2*di)*cos(2*dr)-1)*sin(2*k*x)/(2*pi*sqrt(k*x))+(exp(-2*di)*sin(2*dr))*cos(2*k*x)/(2*pi*sqrt(k*x))', 'independent', 'x', 'problem', 'k');
%ft6 = fittype('(exp(-2*di)*cos(2*dr)-1)*sin(2*k*x)/(2*pi*(k*x)^0.1)+(exp(-2*di)*sin(2*dr))*cos(2*k*x)/(2*pi*(k*x)^0.1)', 'independent', 'x', 'problem', 'k');
%
ft_sim = fittype('A*sin(k*x+phi)/sqrt(k*x)', 'problem', 'k', 'independent', 'x');

ft3a = fittype( 'a*exp(-p1)*sin(2*k*x) + b*exp(-p1)*cos(2*k*x)', 'independent', 'x','problem', 'k');


ind_fit_start = zeros(9,1);
A = zeros(9,1);
B = zeros(9,1);
k_method2 = k_estimated;
k_method3 = kv2k(E,[0.439, 0.4068, -10.996]); %Using new results for dispersion
k_method4 = zeros(9,1);
A2 = zeros(9,1);
B2 = zeros(9,1);
A31 = zeros(9,1);
B31 = zeros(9,1);
deltaI_5 = zeros(9,1);
deltaR_5 = zeros(9,1);
deltaI_5_1 = zeros(9,1);
deltaR_5_1 = zeros(9,1);
simtest = zeros(9,2);

for i = 1:nE
    x_from_wall = (0:(size(linecut_wall,2)-1))*160/1024;
    linecut_wall_fit(i,:) = linecut_wall(i,:) - mean(linecut_wall(i,:));
    [peaks,locs,~, prom] = findpeaks(linecut_wall(i,:));
    
    if peaks(1)>(peaks(2)+0.2) | prom(2)>prom(1) %without this sometimes (4-7) it identified a peak from the wall as the first peak. This makes sure it's a normal sized peak at x=0
        ind_fit_start(i) = locs(2);
    else
        ind_fit_start(i) = locs(1);
    end
%     figure; 
%     subplot(1,2,1)
%     plot(x_from_wall, linecut_wall_fit(i,:))
    
    linecut_wall_fit_final{i} = linecut_wall_fit(i,ind_fit_start(i):end);
    linecut_wall_fit_final{i} = linecut_wall_fit_final{i}-mean(linecut_wall_fit_final{i});
%     subplot(1,2,2)
    figure; scatter(x_from_wall(ind_fit_start(i):end), linecut_wall_fit_final{i})
    
    
    [f_method4,g_method4] = fit(x_from_wall(ind_fit_start(i):end)',linecut_wall_fit_final{i}', ft4, options);
    
    k_method4(i) = 0.5*f_method4.w;
    
    [f_method3,g_method3] = fit(x_from_wall(ind_fit_start(i):end)',linecut_wall_fit_final{i}', ft3, 'problem', k_method3(i), 'StartPoint', [0, 0], 'Lower', lower_bounds, 'Upper', upper_bounds);
    [f_method3_1, g_method3_1] = fit(x_from_wall(ind_fit_start(i):end)',linecut_wall_fit_final{i}', ft3, 'problem', k_method3(i), 'StartPoint', [0, 0]);
    
    [f_method5, g_method5] = fit(x_from_wall(ind_fit_start(i):end)',linecut_wall_fit_final{i}', ft5, 'problem', k_method3(i), 'StartPoint', [0, 0], 'Lower', [0.001, -Inf]);
    [f_method5_1, g_method5_1] = fit(x_from_wall(ind_fit_start(i):end)',linecut_wall_fit_final{i}', ft5, 'problem', k_method3(i), 'StartPoint', [0, 0]);
    
    [f_simtest1, g_simtest1] = fit(x_from_wall(ind_fit_start(i):end)',linecut_wall_fit_final{i}', ft_sim, 'problem', k_method3(i), 'StartPoint', [0, 0], 'Lower', [0.001, -Inf]);
    [f_simtest2, g_simtest2] = fit(x_from_wall(ind_fit_start(i):end)',linecut_wall_fit_final{i}', ft3a, 'problem', k_method3(i), 'StartPoint', [0, 0,0], 'Lower', lower_bounds, 'Upper', upper_bounds);
    
    simtest(i,1) = f_simtest1.A;
    simtest(i,2) = f_simtest1.phi;
    
    hold on
    plot(f_method3,'g')
    plot(f_simtest1, 'm')
    %plot(f_method4,'k')
    plot(f_simtest2, 'r')
    %plot(f_method5)
    plot(f_method5_1,'k')
    ylim([-0.1, 0.1])
    A(i) = f_method3.a;
    B(i) = f_method3.b;
    
    A2(i) = f_method4.a1;
    B2(i) = f_method4.b1; 
    
    tempA = f_method3_1.a;
    tempB = f_method3_1.b;
    
    A31(i) = f_method3_1.a;
    B31(i) = f_method3_1.b;
    
    deltaI_5(i) = f_method5.di;
    deltaR_5(i) = f_method5.dr;
    
    deltaI_5_1(i) = f_method5_1.di;
    deltaR_5_1(i) = f_method5_1.dr;
    
    tempdeltaI1 = log(((2*pi*A(i)+1).^2+(2*pi*B(i)).^2))./(-4);
    tempdeltaI2 = log(((2*pi*tempA+1).^2+(2*pi*tempB).^2))./(-4);
    tempdeltaI3 = log(((2*pi*A2(i)+1).^2+(2*pi*B2(i)).^2))./(-4);
    
    [g_method3.rsquare tempdeltaI1; g_method3_1.rsquare tempdeltaI2; g_method4.rsquare tempdeltaI3; g_method5.rsquare f_method5.di; g_method5_1.rsquare f_method5_1.di]
    
    %g_method3.rsquare
    
end


delta_R = atan((2*pi*B)./(2*pi*A+1))/2;
delta_R2 = atan((2*pi*B2)./(2*pi*A2+1))/2;

delta_I = log(((2*pi*A+1).^2+(2*pi*B).^2))./(-4);
delta_I2 = log(((2*pi*A2+1).^2+(2*pi*B2).^2))./(-4);

avg_delta_R = mean(delta_R);
avg_delta_I = mean(delta_I);

avg_delta_R2 = mean(delta_R2);
avg_delta_I2 = mean(delta_I2);

%Testing to see if the delta's calculated in the various methods above
%match the data when used in simulation
%I'm only going to check the ones that haven't been restricted with limits

delta_R_3 = atan((2*pi*B)./(2*pi*A+1))/2;
delta_R_3_1 = atan((2*pi*B31)./(2*pi*A31+1))/2;

delta_I_3 = log(((2*pi*A+1).^2+(2*pi*B).^2))./(-4);
delta_I_3_1 = log(((2*pi*A31+1).^2+(2*pi*B31).^2))./(-4);

dispersion1 = [0.439, 0.4068, -10.996];
abc = kconstants; 
a0 = abc.a;
n = 30;
mapsize = 160; %Angstroms, to match original data

vpCOwall = [-n/2*2*a0:2*a0:n/2*2*a0; (-mapsize/2+20)*ones(1,n+1)]';
nspec = 301;
specPoints = [zeros(1,nspec); linspace(0,140,nspec)]';
x_sim = linspace(0,140,nspec);
%k_sim = kv2k(E,dispersion1);
i0 = sqrt(-1);


for i = 1:9
    disp_fit_3 = kspec(vpCOwall, specPoints, E(i), (delta_R_3(i)+i0*delta_I_3(i)), dispersion1);
    disp_fit_3_1 = kspec(vpCOwall, specPoints, E(i), (delta_R_3_1(i)+i0*delta_I_3_1(i)), dispersion1);
    disp_fit_5 = kspec(vpCOwall, specPoints, E(i), (deltaR_5(i)+i0*deltaI_5(i)), dispersion1);
    disp_fit_5_1 = kspec(vpCOwall, specPoints, E(i), (deltaR_5_1(i)+i0*deltaI_5_1(i)), dispersion1);
    figure;
    scatter(x_from_wall(ind_fit_start(i):end), linecut_wall_fit_final{i})
    hold on
    plot(linspace(0,140,nspec), disp_fit_3-1)
    
    
    
end

%% Creating a model to predict deltaI and deltaR based on amplitude and phase
close all
dispersion1 = [0.439, 0.4068, -10.996];

ft_sim = fittype('A*sin(k*x+phi)', 'problem', 'k', 'independent', 'x');
%ft3a = fittype( 'a*exp(-p1)*sin(2*k*x)/sqrt(k*x) + b*exp(-p1)*cos(2*k*x)/sqrt(k*x)', 'independent', 'x','problem', 'k');
ft3a = fittype( 'a*exp(-p1)*sin(2*k*x) + b*exp(-p1)*cos(2*k*x)', 'independent', 'x','problem', 'k');


abc = kconstants; 
a0 = abc.a;
n = 30;
mapsize = 160; %Angstroms, to match original data

vpCOwall = [-n/2*2*a0:2*a0:n/2*2*a0; (-mapsize/2+20)*ones(1,n+1)]';
nspec = 301;
specPoints = [zeros(1,nspec); linspace(-60,80,nspec)]';

E = -0.4:0.1:0.4;

training_size = 3;
training1 = cell(training_size,2);

rng('default'); 
vars = rand([2,1,training_size]);
%delta I should be from 0 to 1
%delta R should be between -pi/2 to 0
vars(2,1,:) = (vars(2,1,:)-1)*pi/2;
x_sim = linspace(0,140,nspec);
k_sim = kv2k(E,dispersion1);

nkeep = 5;
training2 = zeros(training_size*size(E,2), (6+4*nkeep));



for i = 1:training_size
    
    
    deltaI = vars(1,1,i);
    deltaR = vars(2,1,i);
    delta = deltaR+sqrt(-1)*deltaI;
    
    training1{i,1} = [deltaI, deltaR];
    training1{i,2} = kspec(vpCOwall, specPoints, E, delta, dispersion1)-1;
    for j = 1:9
        
        [fit_sim, good_sim] = fit(x_sim', training1{i,2}(j,:)', ft_sim, 'problem', k_sim(j),'StartPoint', [0,0]);
        [fit_sim_3, good_sim_3] = fit(x_sim', training1{i,2}(j,:)', ft3a, 'problem', k_sim(j),'StartPoint', [0,0,0],  'Lower', lower_bounds, 'Upper', upper_bounds);
        [pks, locs, w, p] = findpeaks( training1{i,2}(j,:));
        
        
        if length(pks) >= nkeep
            
            peaks_temp = pks(1:nkeep);
            locs_temp = locs(1:nkeep);
            width_temp = w(1:nkeep);
            prom_temp = p(1:nkeep);
        else
            peaks_temp = [pks zeros(1,nkeep-length(pks))];
            locs_temp = [locs zeros(1,nkeep-length(pks))];
            width_temp = [w zeros(1,nkeep-length(pks))];
            prom_temp = [p zeros(1,nkeep-length(pks))];
            
        end
            
        
        
        training2(((i-1)*9+j),:) = [fit_sim.A fit_sim.phi deltaI deltaR E(j) k_sim(j) peaks_temp locs_temp width_temp prom_temp];
        
        figure; scatter(x_sim', training1{i,2}(j,:))
        hold on
        plot(fit_sim_3,'g')
        plot(fit_sim, 'm')

        
    end
    
    i
    
end


%% Saving the training data

save('/Users/lauracollins/Desktop/LineCutTrainingData050318.mat', 'training2');
csvwrite('/Users/lauracollins/Desktop/LineCutTrainingData050318.csv', training2);

%% Making the models

training3 = table(training2(:,1), training2(:,2), training2(:,3), training2(:,4));
training3.Properties.VariableNames = {'A', 'phi', 'deltaI', 'deltaR'};

fitlm([training3.A, training3.phi], training3.deltaI)
fitlm([training3.A, training3.phi], training3.deltaR)

[beta, sigma] = mvregress([training3.A, training3.phi], [training3.deltaI, training3.deltaR])


