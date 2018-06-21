%% Open and reshape Graphene Experimental Data
GrapheneExp = load('ES_AG_Exp_data'); 
GrapheneExp = GrapheneExp(:,:).all_cols;
GrapheneExp = GrapheneExp';