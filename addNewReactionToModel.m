function [ model ] = addNewReactionToModel( in_model, rxnId, rxnName, lb, ub, eccodes, subSystem, substrates, products, gene, protUse )
% add another clone of a reaction for a different enzyme
model = in_model;
rxn_prot = "prot_"+string(model.enzymes(string(model.enzGenes) == gene));
rxn_prot_i = find(string(model.mets) == rxn_prot);
substrates_i = find(ismember(string(model.mets), substrates));
products_i = find(ismember(string(model.mets), products));
% add the new rxn
model.rxns(end+1) = {char(rxnId)};
model.rxnNames{end+1} = rxnName;
model.lb(end+1) = lb;
model.ub(end+1) = ub;
model.rxnConfidenceScores(end+1) = 0;
model.rev(end+1) = 0;
model.c(end+1) = 0;
model.subSystems{end+1} = subSystem;
model.rxnNotes{end+1} = 'freshly added';
model.rxnGeneMat(end+1,:) = 0;
model.eccodes{end+1,:} = eccodes;

model.grRules(end+1) = {gene};

model.S(:,end+1) = 0;
% take prot usage value from cloned reaction, take rel kCat into account
model.S(rxn_prot_i,end) = -protUse;
model.S(substrates_i,end) = -1;
model.S(products_i,end) = 1;

