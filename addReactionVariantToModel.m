function [ model, rxnId ] = addReactionVariantToModel( in_model, number, gene, rxn2clone, arm_rxn, rel_kCat )
% add another clone of a reaction for a different enzyme
model = in_model;
rxn2cloneInx = find(string(model.rxns) == rxn2clone);
rxnId =extractBefore(rxn2clone, "No")+'No'+string(number);
rxnName =extractBefore(model.rxnNames(rxn2cloneInx), "No")+"No"+string(number)+")";
rxn_prot = "prot_"+string(model.enzymes(string(model.enzGenes) == gene));
rxn_prot_i = find(string(model.mets) == rxn_prot);
rxnMets = model.mets(model.S(:,rxn2cloneInx)~=0);
rxn2clone_prot = string(rxnMets(contains(string(rxnMets),'prot_')));
rxn2clone_prot_i = find(string(model.mets) == rxn2clone_prot);

% add the new rxn
model.rxns(end+1) = {char(rxnId)};
model.rxnNames{end+1} = rxnName;
model.lb(end+1) = model.lb(rxn2cloneInx);
model.ub(end+1) = model.ub(rxn2cloneInx);
model.rxnConfidenceScores(end+1) = model.rxnConfidenceScores(rxn2cloneInx);
model.rev(end+1) = model.rev(rxn2cloneInx);
model.c(end+1) = model.c(rxn2cloneInx);
model.subSystems(end+1) = model.subSystems(rxn2cloneInx);
model.rxnNotes(end+1) = model.rxnNotes(rxn2cloneInx);
model.rxnGeneMat(end+1,:) = model.rxnGeneMat(rxn2cloneInx,:);
model.eccodes(end+1,:) = model.eccodes(rxn2cloneInx,:);

model.grRules(end+1) = {gene};

model.S(:,end+1) = model.S(:,rxn2cloneInx);
% take prot usage value from cloned reaction, take rel kCat into account
model.S(rxn_prot_i,end) = model.S(rxn2clone_prot_i,rxn2cloneInx)/rel_kCat;
model.S(rxn2clone_prot_i,end) = 0;

% put reaction behind cloned one?

% add GR to arm_rxn
model.grRules(string(model.rxns) == arm_rxn) = {string(model.grRules(string(model.rxns) == arm_rxn)) + " or " + gene};
