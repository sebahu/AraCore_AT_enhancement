function [ model, arm_rxnId ] = changeToArmReaction( in_model, rxn2split )
% GPAT_hNo1
% split original reaction into arm_ reaction which has original substrates,
% but not the prot_substrate, and an intermediate pmet_<rxn name> as
% product. rxn2split should end in "No1"
model = in_model;
rxn2splitInx = find(string(model.rxns) == rxn2split);
rxn_baseId = extractBefore(rxn2split,"No");
arm_rxnId ="arm_"+rxn_baseId;
rxnName =extractBefore(model.rxnNames(rxn2splitInx), "No")+"arm)";
rxn_intermediate = "pmet_"+rxn_baseId;
rxnMets_i = find(model.S(:,rxn2splitInx)~=0);
rxnMets = model.mets(rxnMets_i);
rxn2split_prot = string(rxnMets(contains(rxnMets,'prot_')));
rxn2clone_prot_i = find(string(model.mets) == rxn2split_prot);
rxnSubstrates_i = rxnMets_i(model.S(rxnMets_i,rxn2splitInx)<0);
rxnProducts_i = setdiff(rxnMets_i, rxnSubstrates_i);
rxnSubstrates_i = setdiff(rxnSubstrates_i,rxn2clone_prot_i);

rxnComp = model.metComps(rxnSubstrates_i(1));
% add intermediate met
model.S(end+1,:) = 0;
model.metFormulas{end+1} = "";
model.metCharges(end+1) = NaN;
model.mets{end+1} = rxn_intermediate;
model.b(end+1) = 0;
model.metNames{end+1} = rxn_intermediate;
model.metComps(end+1) = rxnComp;
model.metMiriams{end+1} = [];
model.metisinchikeyID{end +1} = "";

%clone rxn for the arm_ rxn



% add the new rxn
model.rxns{end+1} = char(arm_rxnId);
model.rxnNames{end+1} = rxnName;
model.lb(end+1) = model.lb(rxn2splitInx);
model.ub(end+1) = model.ub(rxn2splitInx);
model.rxnConfidenceScores(end+1) = model.rxnConfidenceScores(rxn2splitInx);
model.rev(end+1) = model.rev(rxn2splitInx);
model.c(end+1) = model.c(rxn2splitInx);
model.subSystems(end+1) = model.subSystems(rxn2splitInx);
model.rxnNotes(end+1) = model.rxnNotes(rxn2splitInx);
model.rxnGeneMat(end+1,:) = model.rxnGeneMat(rxn2splitInx,:);
model.eccodes(end+1,:) = model.eccodes(rxn2splitInx,:);

model.grRules(end+1) = model.grRules(rxn2splitInx);

model.S(:,end+1) = model.S(:,rxn2splitInx);
% change S to make the split perfect, original reaction now has the
% intermediate (newly added, last metabolite) and the protein as substrates
model.S(rxnSubstrates_i,rxn2splitInx) = 0;
model.S(end,rxn2splitInx) = -1;
% arm_ rxn has original substrates and intermediate as product
model.S(rxnProducts_i,end) = 0;
model.S(rxn2clone_prot_i,end) = 0;
model.S(end,end) = 1;

% put reaction behind cloned one?


