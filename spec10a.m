%%
cd('/Users/ken/Documents/Notre Dame/Research/Projects/Scattering Optimization')

%% Point Specs on Hex lattices
% Measurement with 10a spacing.
load 'Spec10a.mat'
h0=ksmooth(mean(didv0,2),5);
h10b=ksmooth(mean(didvb10,2),5); h10t=ksmooth(mean(didvt10,2),5);
h10br=h10b./h0; h10tr=h10t./h0;
figure; plot(v,[h10br h10tr]);

% sim the hexagonal lattice specs
a = 2.55; nhex = 5; vsim = linspace(-0.4, 0.5, 451)'; 
sf =10; vp = khex(nhex, sf*a,1); vspec=[0,sf*a/sqrt(3); sf*a/2,0]; 
simh10=zeros(size(vsim,1),size(vspec,1));
for ni=1:size(vspec,1)
    simh10(:,ni) = kspec(vp, vspec(ni,:), vsim,(-.15+0.05*sqrt(-1)));
end

figure; plot(v,h10tr,vsim*1000,simh10(:,1));
figure; plot(v,h10br,vsim*1000,simh10(:,2));
