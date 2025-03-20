function [ model ] = addNewEnzymeToModel( in_model, enzyme, gene, miriam, MW, protCost, sequence, comp )
% add another clone of a reaction for a different enzyme
model = in_model;
prot_pool_i = find(string(model.mets)=="prot_pool");

% add the new enzyme
model.enzGenes{end+1} = gene;
model.enzNames{end+1} = enzyme;
model.enzymes{end+1} = miriam.value;
model.geneMiriams{end+1} = miriam;
model.genes{end+1} = gene;
model.geneShortNames{end+1} = gene;
model.MWs(end+1) = MW;
model.sequences{end+1} = sequence;
model.concs(end+1)=NaN;
model.pathways{end +1} = 'TODO';

model.mets{end + 1} = strcat('prot_',miriam.value);
model.metNames{end + 1} = strcat('prot_',miriam.value);
model.metFormulas{end +1} = '';
model.metCharges(end +1) = NaN;
model.metComps(end +1) = comp;
model.metisinchikeyID{end +1} = '';
model.metMiriams{end +1} = [];
model.b(end +1) = 0;


model.rxns{end+1} = strcat('draw_prot_',miriam.value);
model.rxnNames{end+1} =  strcat('draw_prot_',miriam.value);
model.lb(end+1) = 0;
model.ub(end+1) = Inf;
model.rxnConfidenceScores(end+1) = NaN;
model.rev(end+1) = 0;
model.c(end+1) = 0;
model.subSystems{end+1} = '';
model.rxnNotes{end+1} = '';
model.rxnGeneMat(end+1,end +1) = 0;
model.eccodes{end+1,:} = '';

model.grRules(end+1) = {gene};

model.S(end+1,end+1) = 1;
model.S(prot_pool_i, end) = -protCost;

end

