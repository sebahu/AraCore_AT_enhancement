function [ model ] = addNewReactionToNonEcModel( in_model, rxnId, rxnName, lb, ub, eccodes, subSystem, substrates, products, gene )
% add another clone of a reaction for a different enzyme
model = in_model;
substrates_i = find(ismember(string(model.mets), string(substrates)));
products_i = find(ismember(string(model.mets), string(products)));
% add the new rxn
model.rxns(end+1) = {char(rxnId)};
model.rxnNames{end+1} = rxnName;
model.lb(end+1) = lb;
model.ub(end+1) = ub;
model.rxnConfidenceScores(end+1) = 0;
model.c(end+1) = 0;
model.subSystems{end+1} = subSystem;
model.rxnNotes{end+1} = 'freshly added';

model.grRules(end+1) = {gene};

model.S(:,end+1) = 0;
% take prot usage value from cloned reaction, take rel kCat into account
model.S(substrates_i,end) = -1;
model.S(products_i,end) = 1;

