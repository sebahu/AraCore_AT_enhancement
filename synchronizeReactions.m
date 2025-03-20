function [ model ] = synchronizeReactions( in_model, rxn1, rxn2, ratio1_2 )
% GPAT_hNo1
% create a new sync metabolite, which is produced in rxn1 and comsumed in
% rxn2, with a given ratio
model = in_model;
rxn1_i = find(string(model.rxns) == rxn1);
rxn2_i = find(string(model.rxns) == rxn2);

% add intermediate met
model.S(end+1,:) = 0;
model.metFormulas{end+1} = "";
model.metCharges(end+1) = NaN;
model.mets{end+1} = char("sync_"+string(rxn1)+"_"+string(rxn2));
model.b(end+1) = 0;
model.metNames{end+1} = char("sync_"+string(rxn1)+"_"+string(rxn2));
model.metComps(end+1) = 1; % does it matter?
model.metMiriams{end+1} = [];
model.metisinchikeyID{end +1} = "";

% use it in the two reactions
model.S(end,rxn1_i) = -1;
model.S(end,rxn2_i) = ratio1_2;


