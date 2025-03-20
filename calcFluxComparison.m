function [ flux9_fixed0_minE, flux9_like_fixed0_minE ] = calcFluxComparison( ecModel_expanded_allSynced, minimizeE_expanded )
% forced flux threw new reactions AlaTA_cNo1, arm_AlaHyPyrAT_m, AlaHyPyrAT_hNo1, AlaHyPyrAT_pNo1, 
% arm_AlaGlyoxyAT_m, arm_AlaGlyoxyAT_p, AlaGlyoxyAT_hNo1, Ala4MOPAT_hNo1
ecModel_expanded_allSynced_fixed0 = ecModel_expanded_allSynced;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="AlaTA_cNo1")=0.003;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="arm_AlaHyPyrAT_m")=0.006;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="AlaHyPyrAT_hNo1")=0.003;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="AlaHyPyrAT_pNo1")=0.003;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="arm_AlaGlyoxyAT_m")=0.009;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="arm_AlaGlyoxyAT_p")=0.006;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="AlaGlyoxyAT_hNo1")=0.003;
ecModel_expanded_allSynced_fixed0.lb(string(ecModel_expanded_allSynced.rxns)=="Ala4MOPAT_hNo1")=0.0015;
flux9_fixed0_minE = pFBA( ecModel_expanded_allSynced_fixed0, minimizeE_expanded );

% blocked new reactions, but keeping the other fluxes at least at half
% their previous value, to have a similarity
ecModel_expanded_allSynced_like_fixed0 = ecModel_expanded_allSynced_fixed0;
ecModel_expanded_allSynced_like_fixed0.lb = flux9_fixed0_minE / 2;
ecModel_expanded_allSynced_like_fixed0.lb(startsWith(string(ecModel_expanded_allSynced_like_fixed0.rxns),"draw"))=0;
ecModel_expanded_allSynced_like_fixed0.lb(2517:end) = 0;
ecModel_expanded_allSynced_like_fixed0.ub(2517:end) = 0;
flux9_like_fixed0_minE = pFBA( ecModel_expanded_allSynced_like_fixed0, minimizeE_expanded );
end