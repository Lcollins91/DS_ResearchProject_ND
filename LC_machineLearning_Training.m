%% Want to run a bunch of simulations with an assortment of m*, alpha and 
% deltas to use as training data

%Could just save the line cut specs for now in a cell and generate the
%features later. 


abc = kconstants; 
E0 = abc.E0;
nt = 50000; %size of training set
mapsize = 160; %Angstroms, to match original data
n = 30;

vpCOwall = [-n/2*2*a:2*a:n/2*2*a; (-mapsize/2+20)*ones(1,n+1)]';

E = -0.4:0.1:0.4;
%delta = 0.2*(-1 +sqrt(-1));
nspec = 301;
specPoints = [zeros(1,nspec); linspace(-80,80,nspec)]';

%sim1 = kspec(vpCOwall, specPoints, E, delta);


%First we need to generate the m*, alpha and deltas randomly in appropriate
%ranges. 
rng('default'); 

vars = rand([4,1,nt]);
%training = cell(nt,2);

%Rescale the randomly generated variables to be reasonable numbers. 
%m* is probably fine between 0, 1
%alpha probably is too?
%delta should be between -1, 1

vars(3:4,:,:) = vars(3:4,:,:)*2-1;

tic
for i = 1:46875
    
    delta = vars(3,:,i)+sqrt(-1)*vars(4,:,i);
    dispersion = [E0, vars(1,:,i), vars(2,:,i)];
    training{i,1} = vars(:,:,i);
    training{i,2} = kspec(vpCOwall, specPoints, E, delta, dispersion);
    
    i
end
toc
