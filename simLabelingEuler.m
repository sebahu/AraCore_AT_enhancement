function [history_ratio, history_unlabeled, history_labeled, mapping_matrix, high_flux_atoms] = simLabelingEuler(sim_fluxes, ...
    history_length, mapping_matrix, label_feed_rxn_ids, label_feed_atoms, total_atoms, flux_part)
% Performs simulation of an experiment with stable isotope labeling
%
% USAGE:
%
%    [history_ratio, history_unlabeled, history_labeled, mapping_matrix, high_flux_atoms] = simlabeling(sim_fluxes, history_length, mapping_matrix, label_feed_rxn_ids, label_feed_atoms, total_atoms, flux_part)
%
% INPUT:
%    sim_fluxes:       fluxes of the reactions in model to be simulated, in
%                      model.rxns. unit is per (time interval of one
%                      history step)
%    history_length:   number of simulated history steps (* sim_fluxes is
%                      the total flux), for which the labeling values are
%                      stored and returned as result of the function
%    atom_map_mapping: matrix of atom transitions, contains the atom
%                      transitions from atom pool "row" to atom pool
%                      "column"
%    label_feed_rxn_ids: reactions which feed labeled atoms into the network
%    label_feed_atoms: atom pools of input reactions, which are labeled
%    total_atoms:      size of atom pools, for the atom pools referenced in
%                      the matrix
%    flux_part:        how many mini-steps shall be simulated per history step
%

unlabeled_atoms = total_atoms;
labeled_atoms = zeros(size(unlabeled_atoms));

history_unlabeled = zeros(size(unlabeled_atoms,1),history_length+1);
history_labeled = zeros(size(unlabeled_atoms,1),history_length+1);
history_unlabeled(:,1) = total_atoms;
history_labeled(:,1) = labeled_atoms;

%minimal_flux = 1E-10;
%sanitize fluxes
%sim_fluxes(abs(sim_fluxes)<minimal_flux)=0;
base_rxn_flux = sim_fluxes/flux_part;


%precalculating some values, which are reused every step
base_mapping_matrix = sparse(mapping_matrix/flux_part);
total_out_flux = sum(base_mapping_matrix,2);
total_in_flux = sum(base_mapping_matrix(:,1:end-1),1)';
high_flux_atoms = total_out_flux > total_atoms/2;

labeled_atom_import = zeros(size(unlabeled_atoms));
for feed_rxn_i = 1:length(label_feed_rxn_ids)
    dest=label_feed_atoms(feed_rxn_i);
    % relying on selection of feed rxns to have pos. flux
    labeled_atom_import(dest) = labeled_atom_import(dest) + base_rxn_flux(label_feed_rxn_ids(feed_rxn_i));
    if base_rxn_flux(label_feed_rxn_ids(feed_rxn_i)) < 0
        "negative feed! RXN #" + label_feed_rxn_ids(feed_rxn_i)
    end
end


% history_length steps with flux_part mini steps of 1/flux_part poolsize
for history_step = 1:history_length
    for mini_step = 1:flux_part
        labeled_ratios = labeled_atoms./total_atoms;

        % first calculate all outgoing atoms
        labeled_atoms_out = total_out_flux .* labeled_ratios;
        unlabeled_atoms_out = total_out_flux - labeled_atoms_out;
        labeled_atom_gain = base_mapping_matrix(:,1:end-1)'*labeled_ratios;

        unlabeled_atom_gain = total_in_flux - labeled_atom_gain;
        labeled_atom_gain = labeled_atom_gain + labeled_atom_import;
        labeled_atom_change = labeled_atom_gain - labeled_atoms_out;            

        if sum(abs(labeled_atom_change(high_flux_atoms))>0)
            labeled_atom_change(high_flux_atoms) = ((labeled_atom_gain(high_flux_atoms) .* total_atoms(high_flux_atoms)) ...
                ./ total_out_flux(high_flux_atoms)) - labeled_atoms(high_flux_atoms);
        end

        unlabeled_atom_change = unlabeled_atom_gain - unlabeled_atoms_out;
        %unchanged_atoms = abs(labeled_atom_change)<1e-16;
        %labeled_atom_change(unchanged_atoms) = 0;
        %unlabeled_atom_change(unchanged_atoms) = 0;

        unlabeled_atoms = unlabeled_atoms + unlabeled_atom_change;
        labeled_atoms = labeled_atoms + labeled_atom_change;
    end

    history_unlabeled(:,history_step+1) = unlabeled_atoms;
    history_labeled(:,history_step+1) = labeled_atoms;
end

history_ratio = history_labeled ./ total_atoms;

end