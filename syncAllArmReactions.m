function [ model ] = syncAllArmReactions( in_model )
% add another clone of a reaction for a different enzyme
model = in_model;
all_arm_rxns = string(model.rxns(startsWith(string(model.rxns), "arm_")));
for arm_i = 1:length(all_arm_rxns)
   arm_rxn_i = find(string(model.rxns) == all_arm_rxns(arm_i));
   rxn_baseId = extractAfter(all_arm_rxns(arm_i),4)+"No";
   all_sub_rxns_i = find(startsWith(string(model.rxns), rxn_baseId));
   prev_protein_usage = 0;
   prev_rxn = "";
   for sub_i = 1:length(all_sub_rxns_i)
        rxnMets = model.mets(find(model.S(:,all_sub_rxns_i(sub_i))~=0));
        prot = string(rxnMets(contains(string(rxnMets),'prot_')));
        if length(prot) == 1
            prot_i = find(string(model.mets) == prot);
            protein_usage = model.S(prot_i,all_sub_rxns_i(sub_i));
            rxn = model.rxns(all_sub_rxns_i(sub_i));
            if protein_usage ~= 0
                if prev_protein_usage ~= 0
                   factor = protein_usage/prev_protein_usage;                   
                   model = synchronizeReactions(model, prev_rxn, rxn, factor);
                end
                prev_protein_usage = protein_usage;
                prev_rxn = rxn;               
            end
        end
    end    
end

