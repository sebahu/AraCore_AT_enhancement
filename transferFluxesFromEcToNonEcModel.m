function [ fluxes ] = transferFluxesFromEcToNonEcModel( ecModel, nonEcModel, in_fluxes )
    fluxes = zeros(size(nonEcModel.rxns));
    for rxn_i = 1:length(nonEcModel.rxns)
        destRxn = nonEcModel.rxns(rxn_i);
        % is this rxn in ecModel?
        srcRxn_i = find(string(ecModel.rxns) == string(destRxn));
        if ~ isempty(srcRxn_i)
            fluxes(rxn_i) = in_fluxes(srcRxn_i);
            srcRxn_i = find(string(ecModel.rxns) == string(destRxn)+"_REV");
            if ~ isempty(srcRxn_i)
                fluxes(rxn_i) = fluxes(rxn_i)-in_fluxes(srcRxn_i);
            end            
            continue;
        end
        srcRxn_i = find(string(ecModel.rxns) == "arm_"+string(destRxn));
        if ~ isempty(srcRxn_i)
            fluxes(rxn_i) = in_fluxes(srcRxn_i);
            srcRxn_i = find(string(ecModel.rxns) == "arm_"+string(destRxn)+"_REV");
            if ~ isempty(srcRxn_i)
                fluxes(rxn_i) = fluxes(rxn_i)-in_fluxes(srcRxn_i);
            else
                srcRxn_i = find(string(ecModel.rxns) == string(destRxn)+"_REVNo1");
                if ~ isempty(srcRxn_i)
                    fluxes(rxn_i) = fluxes(rxn_i)-in_fluxes(srcRxn_i);
                end
            end
            continue;
        end
        srcRxn_i = find(string(ecModel.rxns) == string(destRxn)+"No1");
        if ~ isempty(srcRxn_i)
            fluxes(rxn_i) = in_fluxes(srcRxn_i);
            srcRxn_i = find(string(ecModel.rxns) == string(destRxn)+"_REVNo1");
            if ~ isempty(srcRxn_i)
                fluxes(rxn_i) = fluxes(rxn_i)-in_fluxes(srcRxn_i);
            end
            continue;
        end
        
        fluxes(rxn_i) = -1111;
    end
end


