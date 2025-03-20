function [ model ] = replaceEnzymeForReaction( in_model, gene, rxn2adapt, rel_kCat, arm_rxn )
% add another clone of a reaction for a different enzyme
model = in_model;
rxn2adaptInx = find(string(model.rxns) == rxn2adapt);
rxn_prot = "prot_"+string(model.enzymes(string(model.enzGenes) == gene));
rxn_prot_i = find(string(model.mets) == rxn_prot);
rxnMets = model.mets(model.S(:,rxn2adaptInx)~=0);
rxn2replace_prot = string(rxnMets(contains(string(rxnMets),'prot_')));
gene2replace = extractAfter(rxn2replace_prot,"prot_");
rxn2replace_prot_i = find(string(model.mets) == rxn2replace_prot);

% add the new rxn
model.grRules{rxn2adaptInx} = replace(model.grRules{rxn2adaptInx},gene2replace,gene);

% take prot usage value from cloned reaction, take rel kCat into account
model.S(rxn_prot_i,rxn2adaptInx) = model.S(rxn2replace_prot_i,rxn2adaptInx)/rel_kCat;
model.S(rxn2replace_prot_i,rxn2adaptInx) = 0;

% put reaction behind cloned one?

% add GR to arm_rxn
if ~isempty(arm_rxn)
    model.grRules{string(model.rxns) == arm_rxn} = replace(model.grRules{string(model.rxns) == arm_rxn},gene2replace,gene);
end
