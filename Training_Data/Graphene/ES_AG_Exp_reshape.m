%% Open and reshape Graphene Experimental Data
GrapheneExp = load('ES_AG_Exp_data'); 
GrapheneExp = GrapheneExp(:,:).all_cols;
GrapheneExp = GrapheneExp';

%% Save 
save('/Users/emory/Documents/GitHub/DS_ResearchProject_ND/Training_Data/Graphene/ES_AG_Exp_data_reshaped.mat', 'GrapheneExp');
csvwrite('/Users/emory/Documents/GitHub/DS_ResearchProject_ND/Training_Data/Graphene/ES_AG_Exp_data_reshaped.csv', GrapheneExp);