function [ opt_obj, fluxes, exitflag, model ] = testKnockOut( in_model, koProteinMets )
model = in_model;
koMet_i = find(ismember(string(model.mets), koProteinMets));
model.ub(sum(model.S(koMet_i,:)~=0)>0)=0;
model.lb(sum(model.S(koMet_i,:)~=0)>0)=0;
[ opt_obj, fluxes, exitflag ] = FBA(model);
