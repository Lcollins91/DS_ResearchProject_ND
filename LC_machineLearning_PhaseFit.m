%% Trying to use machine learning things to optimize dispersion and phase

%dispersion = either [E0, m*, alpha] or [] for usual quadratic fit


% Simulation Settings - CO wall
constants = kconstants;
a = constants.a;
n = 30;
mapsize = 160; %Angstroms, to match original data
vpCOwall = [-n/2*2*a:2*a:n/2*2*a; (-mapsize/2+20)*ones(1,n+1)]';
%vpCOwall2 = [0:2*a:n*2*a; 20*ones(1,n+1)]';



%I'm going to do one and set up a feature engineering thing for that

%Should see if it's faster to use kmap or kspec -- they seem about the same
%time wise and give the same results. (not surprising there)

E = -0.4:0.1:0.4;
delta = 0.2*(-1 +sqrt(-1));
nspec = 101;
specPoints = [zeros(1,nspec); linspace(-80,80,nspec)]';
tic
sim1 = kspec(vpCOwall, specPoints, E, delta);
toc

tic
sim2 = cell(9,1);
for i = 1:9
    sim2{i} = kmap(vpCOwall, mapsize, E(i), 101, delta);
end
toc

figure; plot(linspace(-80,80,nspec), sim1(1,:)); hold on; plot(linspace(-80,80,nspec),sim2{1}(:,50));


[peaks, locs, ws, ps] = findpeaks(sim1(1,:),linspace(-80,80,nspec));


