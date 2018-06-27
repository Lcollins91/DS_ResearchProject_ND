%%
cd('/Users/emory/Documents/GitHub/DS_ResearchProject_ND')

%% Point Specs on Hex lattices
% Measurement with 10a spacing.
load 'Spec10a.mat'
h0=ksmooth(mean(didv0,2),1);
h10b=ksmooth(mean(didvb10,2),1); h10t=ksmooth(mean(didvt10,2),1);
h10br=h10b./h0; h10tr=h10t./h0;
%figure; plot(v,[h10br h10tr]);

% sim the hexagonal lattice specs
sf=0.96;
disp1 = [0.439, 0.4068, -10.996];
disp2 = [0.439, 0.4068*sf^2, -10.996*sf^4];
a = 2.5477; nhex = 5; vsim = linspace(-0.4, 0.5, 451)'; 
sf =10; vp = khex(nhex, sf*a,1); vspec=[0,sf*a/sqrt(3); sf*a/2,0]; 
simh1 = kspec(vp, vspec, vsim,(-0.15+0.05*sqrt(-1)),disp1);
simh2 = kspec(vp, vspec, vsim,(-0.15+0.05*sqrt(-1)),disp2);

figure; plot(v,h10tr,vsim*1000,simh1(:,1),vsim*1000,simh2(:,1));
figure; plot(v,h10br,vsim*1000,simh1(:,2),vsim*1000,simh2(:,2));


%% Optimizing the Phase
% We want to minimize the difference between simulation and data when
% changing the phase. 

% Loading data:
load 'Spec10a.mat'
h0=ksmooth(mean(didv0,2),1);
h10b=ksmooth(mean(didvb10,2),1); h10t=ksmooth(mean(didvt10,2),1);
h10br=h10b./h0; h10tr=h10t./h0;

% Simulation parameters
a = 2.5477; 
nhex = 5; 
E = linspace(-0.4, 0.5, 451)'; 
sf = 0.965;
disp = [0.439, 0.4068*sf^2, -10.996*sf^4];
vp = khex(nhex, 10*a,1);
vspec=[0,10*a/sqrt(3)]; % top site
%vspec=[10*a/2,0];% bond site

% Interp data to match simulation:
data = interp1(v/1000,h10tr,E); % choose between top and bond data

% Range of value of the phase
deltar = -0.2:0.005:-0.01;
deltai = 0.01:0.005:0.1;
ndr=length(deltar);
ndi=length(deltai);
error = zeros(ndi,ndr);
sim=cell(ndi,ndr);
% calculate error for all phase space
for ni=1:ndi
    for nr=1:ndr
        sim{ni,nr} = kspec(vp, vspec, E,(deltar(nr)+deltai(ni)*sqrt(-1)),disp);
        error(ni,nr)=norm(data-sim{ni,nr});
    end
    ndi-ni
end

% Finding the best phase;
[~,minxy]=min(error(:));
[mini, minr]=ind2sub(size(error),minxy);
mindr=deltar(minr); mindi=deltai(mini); 

%% plot error and best spec.
figure;
imagesc(deltar, deltai, error); axis image xy; kcm('gold');
line(mindr, mindi,'marker','x','color','g','markersize',5);
xlabel('deltaR'); ylabel('deltaI'); 
title(['best phase is ' num2str(mindr+mindi*sqrt(-1))]);

figure; 
plot(v,h10tr,E*1000,sim{mini,minr});
legend('data',num2str(mindr+mindi*sqrt(-1)));


