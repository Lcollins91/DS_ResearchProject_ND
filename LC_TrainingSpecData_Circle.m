%% Generating Training Data for Specs with Multiple peaks



%Need to make the structure we will take the specs in first
abc = kconstants;
a0 = abc.a;


xCenter = 0;
yCenter = 0;
theta = 0 : 0.05 : 2*pi;
radius = 15*2*a0;
x = radius * cos(theta) + xCenter;
y = radius * sin(theta) + yCenter;

vpCO = [x',y'];

figure; scatter(vpCO(:,1), vpCO(:,2))

vspec = [0,0];

dispersion1 = [0.439, 0.4068, -10.996];

E = linspace(-0.4, 0.4, 401);

training_size = 3000;
training1 = cell(training_size,2);

rng('default'); 
vars = rand([2,1,training_size]);
%delta I should be from 0 to 1
%delta R should be between -pi/2 to 0
vars(2,1,:) = (vars(2,1,:)-1)*pi/2;
%x_sim = linspace(0,140,nspec);
%k_sim = kv2k(E,dispersion1);

nkeep = 5;
trainingA = zeros(training_size, (5+4*nkeep));



for i = 1:training_size
    
    
    deltaI = vars(1,1,i);
    deltaR = vars(2,1,i);
    delta = deltaR+sqrt(-1)*deltaI;
    
    training1{i,1} = [deltaI, deltaR];
    training1{i,2} = kspec(vpCO, vspec, E, delta, dispersion1);
    %for j = 1:9
        
        %[fit_sim, good_sim] = fit(x_sim', training1{i,2}(j,:)', ft_sim, 'problem', k_sim(j),'StartPoint', [0,0]);
        
    [pks, locs, w, p] = findpeaks(training1{i,2},E');
    
    avgPeaks_temp = mean(pks);
    avgWidth_temp = mean(w);
    avgProm_temp = mean(p);
    
    if length(pks) >= nkeep

        peaks_temp = pks(1:nkeep)';
        locs_temp = locs(1:nkeep)';
        width_temp = w(1:nkeep)';
        prom_temp = p(1:nkeep)';
    else
        peaks_temp = [pks zeros(1,nkeep-length(pks))];
        locs_temp = [locs zeros(1,nkeep-length(pks))];
        width_temp = [w zeros(1,nkeep-length(pks))];
        prom_temp = [p zeros(1,nkeep-length(pks))];

    end



    trainingA(i,:) = [deltaI deltaR avgPeaks_temp avgWidth_temp avgProm_temp peaks_temp locs_temp width_temp prom_temp];

    i

end
    














%% Saving the training data

save('/Users/lauracollins/Desktop/DS_ResearchProject_ND/LineCutTrainingData051418.mat', 'trainingA');
csvwrite('/Users/lauracollins/Desktop/DS_ResearchProject_ND/LineCutTrainingData051418.csv', trainingA);