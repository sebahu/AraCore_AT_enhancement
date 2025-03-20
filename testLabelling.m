function [ history_ratio ] = testLabelling( in_model, fluxes )
atom_name_prefix_length = 2;
atom_N_id_table = readtable('all_atoms.N.sorted.txt', 'ReadVariableNames', false, 'Delimiter', ' ');
atom_names = extractAfter(string(atom_N_id_table.Var1),atom_name_prefix_length);

AtomTransitionRDT_table = readtable("all_mapping.N.sorted.txt","Delimiter"," ", 'ReadVariableNames', false);
atom_map_rxns = string(AtomTransitionRDT_table.Var1);
atom_map_mapping = string(AtomTransitionRDT_table.Var2);

export_rxn_ids = find(startsWith(string(in_model.rxns), ["Ex_", "Bio_"]));
N15_feed_atoms = 532; % NH4[c]:N#1
N15_feed_rxn_ids = 439; % Im_NH4

[mapping_matrix, atom_met_inx] = create_mapping_matrix(in_model, atom_names, atom_map_rxns, atom_map_mapping, fluxes, export_rxn_ids);
atomChange = sum(mapping_matrix,2) - sum(mapping_matrix(:,1:end-1))';

met_pools = 0.02*ones(length(in_model.mets),1); %default: 20nmol/gdw, for umol/gdw*h fluxes
total_atoms = met_pools(atom_met_inx);
simDurationHours = 4;
logsPerHour = 10;

% filter out unconnected atoms
extended_mapping_matrix = mapping_matrix;
extended_mapping_matrix(end+1,:)=0;
connected_atoms = false(length(atom_names)+1,1);
new_connected_atoms = connected_atoms;
new_connected_atoms(end) = true;
while sum(connected_atoms ~= new_connected_atoms)>0
    connected_atoms = new_connected_atoms;
    new_connected_atoms = connected_atoms | sum(extended_mapping_matrix(:,connected_atoms),2);    
end
unconnected_atoms = ~connected_atoms(1:end-1);
connected_mapping_matrix = mapping_matrix;
connected_mapping_matrix(unconnected_atoms,:)=0;
connected_mapping_matrix(:,unconnected_atoms)=0;

% find minimal step
history_length = simDurationHours * logsPerHour;
total_out_rate = sum(connected_mapping_matrix,2)./total_atoms;
max_flux_rate = max(total_out_rate);

steps_per_log = 200*ceil(max_flux_rate/100);

[history_ratio, ~, ~, ~] = ...
    simLabelingEuler(fluxes, history_length, connected_mapping_matrix, ...
    N15_feed_rxn_ids, N15_feed_atoms, total_atoms, steps_per_log);

