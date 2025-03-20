function [ model ] = changeEnzymeCostForRxn(model, enzymeGene, rxn_pattern, cost, relative)
rxn_prot_old = "prot_"+string(model.enzymes(string(model.enzGenes) == enzymeGene));
rxn_prot_old_i = find(string(model.mets) == rxn_prot_old);
pattern_gene_rxns_i=find(contains(string(model.rxns),rxn_pattern) & model.S(rxn_prot_old_i,:)' ~= 0);
% set prot_use very high (value in S for prot_ is set very negative)
if relative
    model.S(rxn_prot_old_i,pattern_gene_rxns_i) = model.S(rxn_prot_old_i,pattern_gene_rxns_i)*cost;
else
    model.S(rxn_prot_old_i,pattern_gene_rxns_i) = cost;
end

end