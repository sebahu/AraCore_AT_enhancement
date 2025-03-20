function [mapping_matrix, atom_met_inx] = create_mapping_matrix(model, ...
    atom_names, atom_map_rxns, atom_map_mapping, last_fluxes, export_rxn_ids)
atom_name_prefix_length = 2;

atom_cpds = extractBefore(atom_names, ":");
atom_met_inx = zeros(size(atom_names));
for met_i = 1:length(model.mets)
    met = string(model.mets(met_i));
    if ismember(met, atom_cpds)
        atom_met_inx(atom_cpds == met) = met_i;
    end
end

mapping_matrix = zeros(length(atom_names),length(atom_names)+1);

for map_i = 1:length(atom_map_rxns) %284:284 %
    rxn = atom_map_rxns(map_i);
    sourceAtom = extractAfter(extractBefore(atom_map_mapping(map_i),"="),atom_name_prefix_length);
    destAtom = extractAfter(extractAfter(atom_map_mapping(map_i),"="),atom_name_prefix_length);
    source_i = (atom_names == sourceAtom);
    dest_i = (atom_names == destAtom);
    if(sum(source_i)+sum(dest_i) ==0)
        continue
    end
    rxn_i = find(model.rxns == rxn);
    switch_source = false;
    flux_direction = 1;
    S_value = model.S(atom_met_inx(source_i), rxn_i);
    if  S_value > 0
        switch_source = true;
    else
        S_value = -S_value;
    end
    % Workaround: we have two different cases, where S_value is different
    % from 1: a) two (or more) molecules of the same metabolite are involved. Then,
    % there should be as many mappings, thus we need to use 1 here b) a
    % fractional amount is given, as in biomass reaction - we take the
    % fraction - workaround solution: if it is a whole number, take 1
    if S_value == round(S_value)
        S_value=1;
    end
    if last_fluxes(rxn_i) < 0
        switch_source = ~switch_source;
        flux_direction = -1;
    end
    if switch_source
        tmp = source_i;
        source_i = dest_i;
        dest_i = tmp;
    end
    mapping_matrix(source_i,dest_i) = mapping_matrix(source_i,dest_i)+S_value*flux_direction*last_fluxes(rxn_i);
end

for rxn_i = 1:length(export_rxn_ids)
    export_rxn_id = export_rxn_ids(rxn_i);
    met_flux = model.S(:,export_rxn_id)*last_fluxes(export_rxn_id);
    export_mappings = atom_map_mapping(atom_map_rxns == model.rxns(export_rxn_id));
    atoms_with_mapping = string(zeros(2*length(export_mappings),1));
    for map_i = 1:length(export_mappings)
        atoms_with_mapping(2*map_i-1) = extractAfter(extractBefore(export_mappings(map_i),"="),atom_name_prefix_length);
        atoms_with_mapping(2*map_i) = extractAfter(extractAfter(export_mappings(map_i),"="),atom_name_prefix_length);
    end
    atom_flux = met_flux(atom_met_inx);
    atom_flux(ismember(atom_names, atoms_with_mapping))=0;    
    mapping_matrix(atom_flux < 0,end) = mapping_matrix(atom_flux < 0,end) - atom_flux(atom_flux < 0);
end


end