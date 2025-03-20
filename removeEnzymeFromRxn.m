function [ model ] = removeEnzymeFromRxn(model, enzymeGene, rxn_pattern)
model = changeEnzymeCostForRxn(model, enzymeGene, rxn_pattern, 1000, false);
% TODO also change grRules
end