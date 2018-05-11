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
    linecut(i,:) = ksmooth(amplitude*(cut1/cut1(100)+(i-1)*offset),3);
    linecut2(i,:) = linecut(i, 121:end);
    plot(linecut2(i,:)); 
    hold on;
end

hold off

%% Want to fit a cosine function so we can use the wavelength to determine k at each energy
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


%% Fitting a sine wave to the linecut

close all
wavelength = zeros(1,nE);
k_estimated = zeros(1,nE);
k_estimated2 = zeros(1,nE-1);

for i = 1:nE
    Y = linecut4{i};
    X = (0:(length(Y)-1))*(160/1024); %Angstroms
%     figure; 
%     scatter(X,Y)
%     hold on
%     
    [f,g] = fit(X',Y','sin1');
    [fa, ga] = fit(X(2:end)', Y(2:end)', 'sin1');
%     plot(f)
    k_estimated(i) = f.b1; 
    if i>1
        k_estimated2(i-1) = 0.5*fa.b1;
    end
end
%
E = -0.4:0.1:0.4;
k_estimated = 0.5*k_estimated; %Since it should be fitting to sin(2kr)
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
E0 = abc.E0;
E2 = E + E0; %explicitly use known E0 as constant by offsetting E - for a few of the fits only

%Explicitly set E0
ft1 = fittype({'x^2', 'x^4'});% here x is k

%Allow E0 to be fit as well
ft2 = fittype({'1', 'x^2', 'x^4'});

%% Case 1
%Here we fit k's from all energies, with E0 fixed
[f_case1,g_case1] = fit(k_estimated', E2', ft1);

%Results for all k's, E0 fixed
alpha_case1 = f_case1.b; 
mstar_case1 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case1.a);

%% Case 2
%Here we fit k's from all energies, with E0 free
[f_case2, g_case2] = fit(k_estimated', E', ft2);

%Results for all k's, E0 free
E0_case2 = f_case2.a;
alpha_case2 = f_case2.c;
mstar_case2 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case2.b);
%% Case 3
%Here we fit k's from E>=-0.3, with E0 fixed
[f_case3, g_case3] = fit(k_estimated2', E2(2:end)', ft1);

%Results for k's with E >= -0.3, E0 fixed
alpha_case3 = f_case3.b;
mstar_case3 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case3.a);

%% Case 4
%Here we fit k's from E>= -0.3, with E0 free
[f_case4, g_case4] = fit(k_estimated2', E(2:end)', ft2);

%Results for k's with E >= -0.3, E0 free
E0_case4 = f_case4.a;
alpha_case4 = f_case4.c;
mstar_case4 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case4.b);
%% Case 5
%Here we fit David's values for k from the Bessel function fit for E
%>=-0.3, with E0 fixed
[f_case5, g_case5] = fit((k_david)', E2(2:end)', ft1);

%Results for David's k's (E >= -0.3), E0 fixed
alpha_case5 = f_case5.b; 
mstar_case5 = (abc.hbar^2*abc.ec*10^20)/(2*abc.me*f_case5.a);

%% Case 6
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

%% Results Summary

%In the end, we use the average of the results from Cases 4 and 6, where
%the slight differences in the momentum fits arise due to the difference in
%the equation used to fit them (sin(2kr) vs. Bessel)

E0_result = (E0_case4 + E0_case6)/2; %E0 = -0.439 +- 0.001
alpha_result = (alpha_case4 + alpha_case6)/2; % alpha = 10.996 +- 0.504
mstar_result = (mstar_case4 + mstar_case6)/2; % mstar = 0.4068 +- 0.001



%% Summary Plot
figure;
scatter(k_estimated, E2);

hold on
scatter(k_original, E2);
plot(f_case1)

scatter(k_david, E2(2:end))
plot(f_case5)





