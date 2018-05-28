%% Hexagon Spec Data

fn1 = '/Users/lauracollins/Desktop/Data/2016-08-01 Hexaxgonal Corral/';

cd(fn1)

dv5 = 5*10^-3;
mpp = 0.08;

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

close all
figure;
plot(bias, spec3+1)
hold on
plot(bias, spec1)
 
figure;
findpeaks(spec1, 'MinPeakProminence', mpp)

figure;
findpeaks(spec3, 'MinPeakProminence', mpp)

%% Want to compare the experimental specs to simulated
% First peaks are a bit small, and there may be additional peaks in the
% experiment according to Ken

close all

hex1 = hexagon(1);
hex2 = hexagon2(1);

x = linspace(-90,90,256);

%Want to check that we use the same size hexagon as we built
topofile = ksxm('Topo021.sxm');
figure;
imagesc(x,x,topofile.Zf)
axis image 
hold on
scatter(hex1(:,1), hex1(:,2))


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

%scan angle from topograph is 128, rotating back hex_exp to see what the
%spacing is. should line up with hex1 then

% theta = -30*pi/180;
% R = [cos(theta), sin(theta); -sin(theta), cos(theta)];
% 
% hex_exp2 = hex_exp*R;
% 
% scatter(hex_exp2(:,1), hex_exp2(:,2))
























