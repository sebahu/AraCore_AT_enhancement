clearvars
addpath(fullfile('~/apps/gurobi903','linux64', 'matlab'));
addpath('~/apps/cobratoolbox/');
initCobraToolbox(0);
changeCobraSolver('gurobi', 'all');

load('ecAraCore_batch_raw.mat')

% First Part: creating a new version of the model, including all the newly
% found side reactions


% AlaTA_h and reversed get assigned to AT1G80600 , even if the found
% activity was weak, with a corresponding high enzyme cost
additional_AlaTA_h_enzyme = string(ecModel_batch.enzymes(string(ecModel_batch.enzGenes)=="AT1G80600"));
ecModel_batch.S(string(ecModel_batch.mets)=="prot_"+additional_AlaTA_h_enzyme, [1,266]) = [1.8E-4,7E-4];

% starting with the model with additional enzymes
ecModel_enhanced = ecModel_batch;

% add new enzymes/genes: BCAT7 AT1G50090; TAA1 AT1G70560; AGT3 AT2G38400;
% ALD1 AT2G13810
ecModel_enhanced = addNewEnzymeToModel( ecModel_enhanced, 'BCAT7', 'AT1G50090', struct("name", {'uniprot'}, "value", {'Q9LPM8'}), 79.99, 164.3, ...
strcat('MAPSVHPSSSPLFTSKADEKYANVKWDELGFALVPTDYMYVAKCKQGESFSTGEIVPYGDISISPCAGILNYGQGLFEGLKAYRTEDGRITLFRPDQNAIRMQTGADRLCMTPPSPEQFVEAVKQT', ...
    'VLANNKWVPPPGKGALYIRPLLIGTGAVLGVASAPEYTFLIYTSPVGNYHKASSGLNLKVDHNHRRAHFGGTGGVKSCTNYSPVVKSLIEAKSSGFSDVLFLDAATGKNIEEVSTCNIFILKGNIVSTP', ...
    'PTSGTILPGITRKSICELARDIGYEVQERDLSVDELLEAEEVFCTGTAVVIKAVETVTFHDKRVKYRTGEEAFSTKLHLILTNIQMGVVEDKKGWMMEIDHLVGTDSFPDET'), 1);
ecModel_enhanced = addNewEnzymeToModel( ecModel_enhanced, 'TAA1', 'AT1G70560', struct("name", {'uniprot'}, "value", {'Q9S7N2'}), 89.60, 201.34, ...
strcat('MVKLENSRKPEKISNKNIPMSDFVVNLDHGDPTAYEEYWRKMGDRCTVTIRGCDLMSYFSDMTNLCWFLEPELEDAIKDLHGVVGNAATEDRYIVVGTGSTQLCQAAVHALSSLARSQPVSVVAAA', ...
    'PFYSTYVEETTYVRSGMYKWEGDAWGFDKKGPYIELVTSPNNPDGTIRETVVNRPDDDEAKVIHDFAYYWPHYTPITRRQDHDIMLFTFSKITGHAGSRIGWALVKDKEVAKKMVEYIIVNSIGVSKES', ...
    'QVRTAKILNVLKETCKSESESENFFKYGREMMKNRWEKLREVVKESDAFTLPKYPEAFCNYFGKSLESYPAFAWLGTKEETDLVSELRRHKVMSRAGERCGSDKKHVRVSMLSREDVFNVFLERLANMKLIKSIDL'), 1);
ecModel_enhanced = addNewEnzymeToModel( ecModel_enhanced, 'AGT3', 'AT2G38400', struct("name", {'uniprot'}, "value", {'Q94AL9'}), 103.82, 213.4, ...
strcat('MQRFAAKRSVQNISVSLWRRCISSTSQAATASVKDSDEFQARLPPFAYTPPPYTGPSADVILSKRKEFLSPSMFCLYRKPLNIVDGKMQYLFDESGRRYLDAFAGIAVVNCGHCHPDVVEPVINQI', ...
    'KRLQHPTVLYLNHAIADFSEALASKLPGDLKVVFFTNSGTEANELALMMAKLYTGCQDIVAVRNGYHGNAAATMGATGQSMWKFNVVQNSVHHALNPDPYRGVFGSDGEKYAKDLQDLIQYGTTGHIAG', ...
    'FICEAIQGVGGIVELAPGYLSAAYDTVKKAGGLFIADEVQSGFARTGNFWGFEAHNVVPDIVTMAKGIGNGFPLGAVVTTPEIAGVLTRRSYFNTFGGNSVSTTAGLAVLNVIEKEKLQENAAMVGSYL', ...
    'KEKLTQLKEKHEIIGDVRGRGLMLGVELVSDRKLKTPATAETLHIMDQMKELGVLIGKGGYFGNVFRITPPLCFTKDDADFLVEAMDYSMSKM'), 5);
ecModel_enhanced = addNewEnzymeToModel( ecModel_enhanced, 'ALD1', 'AT2G13810', struct("name", {'uniprot'}, "value", {'Q9ZQI7'}), 101.15, 205.6, ...
strcat('MVSLMFFSSASPLCSSPSKIPKASLDFEMKKLGGSTKLVRNVNLEKLKNNYLFPEINRRELEHIEKHPNVQLISLGTGDTTEPIPEQITSHMSNFAHGLSTVEGYRGYGLEQGNKTLRKAIAETFY', ...
    'RDLHVKSNEVFVSDGAQSDISRLQLLLGSNVTIAVQDPTFPAYIDSSVIIGQTGHFHEKTKKYQNVVYMPCGPNNSFFPDLAMTPRTDVIFFCSPNNPTGYVASRKQLHQLVDFAKTNGSIIIFDSAYA', ...
    'AFIEDGSPRSIYEIPGAREVAIEVSSFSKFAGFTGVRLGWSIIPDELLYSNGFPIINDFHRIVTTSFNGASNIAQAGGLACLSSGGLKEIRSVNNYYKENRKILMDTLVSLGLKVYGGVNAPYLWVHFK', ...
    'GSKSWDVFNEILENTHIITVPGSGFGPGGEEYLRISGFGRRDHIVEASKRLQNFFNTRTKHFTYLSSTSNTN'), 6);

% No pbserved rev activity
% for rxn BAATA1 add enzyme BCAT7 = AT1G50090 - new part of AraCore
BAATA1_additionalGenes = ["AT1G50090"];
BAATA1_additional_relKcats = [0.42];
for i = (1:length(BAATA1_additionalGenes))
    ecModel_enhanced = addReactionVariantToModel(ecModel_enhanced, 3+i, BAATA1_additionalGenes(i), "BAATA1_hNo3", "arm_BAATA1_h", BAATA1_additional_relKcats(i));
end

% for rxn BAATA2 add enzyme BCAT7 = AT1G50090 - new part of AraCore
BAATA2_additionalGenes = ["AT1G50090"];
BAATA2_additional_relKcats = [0.19];
for i = (1:length(BAATA2_additionalGenes))
    ecModel_enhanced = addReactionVariantToModel(ecModel_enhanced, 3+i, BAATA2_additionalGenes(i), "BAATA2_hNo3", "arm_BAATA2_h", BAATA2_additional_relKcats(i));
end

% for rxn BAATA3 add enzyme BCAT7 = AT1G50090 - new part of AraCore -
% instead of BCAT3
ecModel_enhanced = replaceEnzymeForReaction(ecModel_enhanced, "AT1G50090", "BAATA3_hNo1", 1.0, "arm_BAATA3_h");

% for AlaTA, remove AT3G08860
ecModel_enhanced = removeEnzymeFromRxn(ecModel_enhanced, "AT3G08860", "AlaTA");

% Now adding completely new reactions

ecModel_expanded = ecModel_enhanced;
% new rxn , GGAT in m, for AT1G17290 (AlaAT1) - but no GLX_m
% ecModel_expanded = addNewReactionToModel( ecModel_expanded, '"GGAT_mNo1"', 'GGAT_m (No1)', 0, 1000, ...
%                                        '', '', ["Glu_m"; "GLX_m"] , ["KG_m"; "Gly_m"], 'AT1G17290', 4E-5 );      
% There is already GGAT_h - no enzyme yet
% new rxns , GGAT in h, for  "AT4G35630" (PSAT1), "AT2G17630" (PSAT2), "AT1G80600" (WIN1);
rxn_prot = "prot_"+string(ecModel_expanded.enzymes(string(ecModel_expanded.enzGenes) == "AT4G35630"));
rxn_prot_i = find(string(ecModel_expanded.mets) == rxn_prot);
rxn_i = find(string(ecModel_expanded.rxns) == "GGAT_h");
ecModel_expanded.rxns(rxn_i) = {'GGAT_hNo1'};
ecModel_expanded.rxnNames{rxn_i} = 'GGAT_h (No1)';
ecModel_expanded.grRules(rxn_i) = {'AT4G35630'};
ecModel_expanded.S(rxn_prot_i,rxn_i) = -2.4E-5;

ecModel_expanded = changeToArmReaction( ecModel_expanded, "GGAT_hNo1");
ecModel_expanded = addReactionVariantToModel(ecModel_expanded, 2, "AT2G17630", "GGAT_hNo1", "arm_GGAT_h", 1);
ecModel_expanded = addReactionVariantToModel(ecModel_expanded, 3, "AT1G80600", "GGAT_hNo1", "arm_GGAT_h", 0.33);

ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'AlaTA_cNo1', 'AlaTA_c (No1)', 0, 1000, ...
                                        '', '', ["Glu_c"; "Pyr_c"] , ["KG_c"; "Ala_c"], 'AT1G70560', 6.5E-5 );

%Ala_m, Ala_p HPR_m	Pyr	Ser  Ala:hydroxypyruvate AT   AGT1 (per), PYD4 (m), ALD1 (p), AGT3 (m) 
ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'AlaHyPyrAT_mNo1', 'Alanine:hydroxypyruvate AT (mito) (No1)', 0, 1000, ...
                                        '', '', ["Ala_m"; "HPR_m"] , ["Pyr_m"; "Ser_m"], 'AT3G08860', 4E-5 );
ecModel_expanded = changeToArmReaction( ecModel_expanded, 'AlaHyPyrAT_mNo1');
ecModel_expanded = addReactionVariantToModel(ecModel_expanded, 2, "AT2G38400", "AlaHyPyrAT_mNo1", "arm_AlaHyPyrAT_m", 1.25);
ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'AlaHyPyrAT_hNo1', 'Alanine:hydroxypyruvate AT (chloroplast) (No1)', 0, 1000, ...
                                        '', '', ["Ala_h"; "HPR_h"] , ["Pyr_h"; "Ser_h"], 'AT2G13810', 2E-4 );
ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'AlaHyPyrAT_pNo1', 'Alanine:hydroxypyruvate AT (peroxisome) (No1)', 0, 1000, ...
                                        '', '', ["Ala_p"; "HPR_p"] , ["Pyr_p"; "Ser_p"], 'AT2G13360', 6.5E-5 );

%Ala_m, Ala_p GLX_m	Pyr	Gly  Ala:glyoxylate AT
% AT3G08860 AT2G13360 AT3G22200   AT1G23310	AT1G70580	AT1G17290 AT2G38400
ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'AlaGlyoxyAT_mNo1', 'Alanine:glyoxylate AT (mito) (No1)', 0, 1000, ...
                                        "", "", ["Ala_m"; "GLX_m"] , ["Pyr_m"; "Gly_m"], 'AT3G08860', 2.4E-5 );
ecModel_expanded = changeToArmReaction( ecModel_expanded, 'AlaGlyoxyAT_mNo1');
additional_cell21_genes = ["AT3G22200"; "AT1G17290"; "AT2G38400"];
additional_cell21_relkatCs = [1,0.5,1];
for i=1:length(additional_cell21_genes)
    ecModel_expanded = addReactionVariantToModel(ecModel_expanded, 1+i, additional_cell21_genes(i), "AlaGlyoxyAT_mNo1", ...
                                                    "arm_AlaGlyoxyAT_m", additional_cell21_relkatCs(i));
end
ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'AlaGlyoxyAT_pNo1', 'Alanine:glyoxylate AT (peroxisome) (No1)', 0, 1000, ...
                                        "", "", ["Ala_p"; "GLX_p"] , ["Pyr_p"; "Gly_p"], 'AT1G23310', 2.4E-5 );
ecModel_expanded = changeToArmReaction( ecModel_expanded, 'AlaGlyoxyAT_pNo1');
ecModel_expanded = addReactionVariantToModel(ecModel_expanded, 2, "AT1G70580", "AlaGlyoxyAT_pNo1", "arm_AlaGlyoxyAT_p", 1);

ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'AlaGlyoxyAT_hNo1', 'Alanine:glyoxylate AT (chloroplast) (No1)', 0, 1000, ...
                                        "", "", ["Ala_h"; "GLX_h"] , ["Pyr_h"; "Gly_h"], 'AT2G13360', 2.4E-5 );

% Ala:4MOP AT
ecModel_expanded = addNewReactionToModel( ecModel_expanded, 'Ala4MOPAT_hNo1', 'Alanine:4MOP AT  (chloroplast) (No1)', 0, 1000, ...
                                        "", "", ["Ala_h"; "4MOP_h"] , ["Pyr_h"; "Leu_h"], 'AT2G13810', 12.0E-5  );

%AspAT2_c
% remove: {'AT5G19550'} {'AT2G22250'} {'AT5G11520'}
%AspAT2_h
% remove:  {'AT5G19550'} {'AT2G22250'} {'AT5G11520'}
%AspAT2_m
% remove:  {'AT2G30970'}
%AspAT2_p
% remove: {'AT5G11520'}
% total removed: {'AT5G19550'} {'AT2G22250'} {'AT5G11520'} {'AT2G30970'}
ecModel_expanded = removeEnzymeFromRxn(ecModel_expanded, "AT5G19550", "AspAT2");
ecModel_expanded = removeEnzymeFromRxn(ecModel_expanded, "AT2G22250", "AspAT2");
ecModel_expanded = removeEnzymeFromRxn(ecModel_expanded, "AT5G11520", "AspAT2");
ecModel_expanded = removeEnzymeFromRxn(ecModel_expanded, "AT2G30970", "AspAT2");

% total kept: {'AT1G62800'} (cost*10), {'AT1G62960'}  {'AT5G51690'} {'AT4G31990'}
ecModel_expanded = changeEnzymeCostForRxn(ecModel_expanded, "AT1G62800", "AspAT2", 10, true);

%adding synchronization for arm reactions according to their enzyme costs
ecModel_expanded_allSynced = syncAllArmReactions(ecModel_expanded);

% Starting the analysis, 1st max biomass production under the same enzyme
% cost constraints
range = zeros(size(ecModel_expanded.rxns));
range5 = range;
[minFlux, maxFlux, exit_code, opt_bio]=slowFVA(ecModel_batch,0.99);
range(maxFlux > 0) = 1 - minFlux(maxFlux > 0) ./ maxFlux(maxFlux > 0);

[minFlux5, maxFlux5, exit_code, opt_bio_expanded]=slowFVA(ecModel_expanded,0.99);
range5(maxFlux5 > 0) = 1 - minFlux5(maxFlux5 > 0) ./ maxFlux5(maxFlux5 > 0);

% biomass reactions have indexes 232 to 234
disp("Adding the new reactions increased the maximum biomass production by "+ ...
    (sum(maxFlux5(232:234))/sum(maxFlux(232:234))-1)*100 + "%.")

disp("Of the newly added reactions, "+join(string(ecModel_expanded.rxnNames(2507+find(minFlux5(2508:end)>0))),", ")+ ...
    " have a non-zero minimum flux, indicating their importance in the increased biomass prodcution.");
disp("Reactions that have a changed range (more than +- 20%) with the added reactions "+...
    "AND have now a non-zero minimum flux are: "+ ...
join(string(ecModel_batch.rxnNames(abs(1-range5(1:1835)./range(1:1835))>0.2 & minFlux(1:1835)==0 & minFlux5(1:1835)>0)),", ")+ ...
".");
disp("Reactions that have a changed range (more than +- 20%) with the added reactions "+...
    "AND can now have zero flux for optimal biomass, while they had a non-zero minimum flux before, are: "+ ...
join(string(ecModel_batch.rxnNames(abs(1-range5(1:1835)./range(1:1835))>0.2 & minFlux(1:1835)>0 & minFlux5(1:1835)==0)),", ")+ ...
".");

% checking changes in enrichment, creating a defined flux distribution for
% the expanded model, minimizing enzyme cost and total flux sum (to
% regulate reactions without enzyme cost)
minimizeE_expanded = 0.0001*ones(size(ecModel_expanded_allSynced.rxns));
minimizeE_expanded(2507)=1;

[ flux9_fixed0_minE, flux9_like_fixed0_minE ] = calcFluxComparison( ecModel_expanded_allSynced, minimizeE_expanded );

load('AraCore_v2_1.mat');

nonEcModel = model;

% converting the ecModel-flux distribution to a non-EC-model flux
% distribution for enrichment simulation
% new reactions: AlaTA_cNo1, arm_AlaHyPyrAT_m, AlaHyPyrAT_hNo1, AlaHyPyrAT_pNo1, 
% arm_AlaGlyoxyAT_m, arm_AlaGlyoxyAT_p, AlaGlyoxyAT_hNo1, Ala4MOPAT_hNo1
nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel,          'AlaTA_c', 'AlaTA_c',           0, 1000, '', '', ["Glu[c]"; "Pyr[c]"] , ["KG[c]"; "Ala[c]"], 'AT1G70560' );

nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel_expanded, 'AlaHyPyrAT_m', 'Alanine:hydroxypyruvate AT (mito)', 0, 1000, '', '', ["Ala[m]"; "HPR[m]"] , ["Pyr[m]"; "Ser[m]"], 'AT3G08860' );
nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel_expanded, 'AlaHyPyrAT_h', 'Alanine:hydroxypyruvate AT (chloroplast)', 0, 1000, '', '', ["Ala[h]"; "HPR[h]"] , ["Pyr[h]"; "Ser[h]"], 'AT2G13360' );
nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel_expanded, 'AlaHyPyrAT_p', 'Alanine:hydroxypyruvate AT (peroxisome)', 0, 1000, '', '', ["Ala[p]"; "HPR[p]"] , ["Pyr[p]"; "Ser[p]"], 'AT2G13360' );

nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel_expanded, 'AlaGlyoxyAT_m', 'Alanine:glyoxylate AT (mito)', 0, 1000, "", "", ["Ala[m]"; "GLX[m]"] , ["Pyr[m]"; "Gly[m]"], 'AT3G08860' );
nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel_expanded, 'AlaGlyoxyAT_h', 'Alanine:glyoxylate AT (chloroplast)', 0, 1000, "", "", ["Ala[h]"; "GLX[h]"] , ["Pyr[h]"; "Gly[h]"], 'AT2G13360' );
nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel_expanded, 'AlaGlyoxyAT_p', 'Alanine:glyoxylate AT (peroxisome)', 0, 1000, "", "", ["Ala[p]"; "GLX[p]"] , ["Pyr[p]"; "Gly[p]"], 'AT1G23310' );

nonEcModel_expanded = addNewReactionToNonEcModel( nonEcModel_expanded, 'Ala4MOPAT_h', 'Alanine:4MOP AT (peroxisome)', 0, 1000, "", "", ["Ala[h]"; "4MOP[h]"] , ["Pyr[h]"; "Leu[h]"], 'AT2G13810' );

% vary fluxes +- 20% using CHRR by setting lb und ub of fluxes
num_samples = 100;

atom_name_prefix_length = 2;
atom_N_id_table = readtable('all_atoms.N.sorted.txt', 'ReadVariableNames', false, 'Delimiter', ' ');
atom_names = extractAfter(string(atom_N_id_table.Var1),atom_name_prefix_length);
decomp_atom_names = extractBefore(atom_names,"[")+extractAfter(atom_names,"]");
decomp_atom_names_unique = unique(decomp_atom_names);
    
fluxes_nonEc_expanded_fixed0_minE = transferFluxesFromEcToNonEcModel(ecModel_expanded_allSynced, nonEcModel_expanded, flux9_fixed0_minE);
fluxes_nonEc_expanded_like_fixed0_minE = transferFluxesFromEcToNonEcModel(ecModel_expanded_allSynced, nonEcModel_expanded, flux9_like_fixed0_minE);
max_random_variation = 0.2;

fluxes_nonEc_expanded_fixed0_minE_samples = zeros(size(fluxes_nonEc_expanded_fixed0_minE,1),num_samples);
fluxes_nonEc_expanded_fixed0_minE_samples(:,1) = fluxes_nonEc_expanded_fixed0_minE;
fluxes_nonEc_expanded_like_fixed0_minE_samples = zeros(size(fluxes_nonEc_expanded_like_fixed0_minE,1),num_samples);
fluxes_nonEc_expanded_like_fixed0_minE_samples(:,1) = fluxes_nonEc_expanded_like_fixed0_minE;

model_sample = nonEcModel_expanded;
model_sample.lb=(1-max_random_variation)*fluxes_nonEc_expanded_fixed0_minE;
model_sample.ub=(1+max_random_variation)*fluxes_nonEc_expanded_fixed0_minE;
switch_boundarys = model_sample.lb > model_sample.ub;
tmp = model_sample.lb(switch_boundarys);
model_sample.lb(switch_boundarys) = model_sample.ub(switch_boundarys);    
model_sample.ub(switch_boundarys) = tmp;
[P,model_post] = chrrParseModel(model_sample);
model_post.c(:)=1;
[samples, roundedPolytope, minSampledFlux, maxSampledFlux] = chrrExpSampler(model_post, 731, num_samples);
fluxes_nonEc_expanded_fixed0_minE_samples = samples;

model_sample = nonEcModel_expanded;
model_sample.lb=(1-max_random_variation)*fluxes_nonEc_expanded_like_fixed0_minE;
model_sample.ub=(1+max_random_variation)*fluxes_nonEc_expanded_like_fixed0_minE;
switch_boundarys = model_sample.lb > model_sample.ub;
tmp = model_sample.lb(switch_boundarys);
model_sample.lb(switch_boundarys) = model_sample.ub(switch_boundarys);    
model_sample.ub(switch_boundarys) = tmp;
[P,model_post] = chrrParseModel(model_sample);
model_post.c(:)=1;
[samples, roundedPolytope, minSampledFlux, maxSampledFlux] = chrrExpSampler(model_post, 731, num_samples);
fluxes_nonEc_expanded_like_fixed0_minE_samples = samples;

for sample_i = 1:num_samples
    hist_nonEc_expanded_fixed0_minE = testLabelling(nonEcModel_expanded, fluxes_nonEc_expanded_fixed0_minE_samples(:,sample_i));
    hist_nonEc_expanded_like_fixed0_minE = testLabelling(nonEcModel_expanded, fluxes_nonEc_expanded_like_fixed0_minE_samples(:,sample_i));
    
    if sample_i == 1
        decomp_hist_nonEc_expanded_like_fixed0_minE = zeros(length(decomp_atom_names_unique), size(hist_nonEc_expanded_like_fixed0_minE,2),num_samples);
        decomp_hist_nonEc_expanded_fixed0_minE = zeros(length(decomp_atom_names_unique), size(hist_nonEc_expanded_fixed0_minE,2),num_samples);
    end

    for i = 1:length(decomp_atom_names_unique)
        decomp_atom = decomp_atom_names_unique(i);
        comp_atom_i = find(decomp_atom_names == decomp_atom);
        if length(comp_atom_i) > 1
            decomp_hist_nonEc_expanded_like_fixed0_minE(i,:,sample_i) = mean(hist_nonEc_expanded_like_fixed0_minE(comp_atom_i,:));
            decomp_hist_nonEc_expanded_fixed0_minE(i,:,sample_i) = mean(hist_nonEc_expanded_fixed0_minE(comp_atom_i,:));
        else
            decomp_hist_nonEc_expanded_like_fixed0_minE(i,:,sample_i) = hist_nonEc_expanded_like_fixed0_minE(comp_atom_i,:);
            decomp_hist_nonEc_expanded_fixed0_minE(i,:,sample_i) = hist_nonEc_expanded_fixed0_minE(comp_atom_i,:);
        end
    end
end
% A-Glu - 17, AGN - 44, Ala - 68, A-Orn - 19;20, Arg - 73-76, Arg-SCA - 69-72, A-Ser - 21, Asp - 80, Cys - 107, DAP - 108;109, 
% GluP - 176, Gly - 177, His - 184-186, Ile - 202, Leu - 207, NH4 - 255,
% Ser - 303, Val - 328
% maybe just add the other atoms as well, select then afterwards
atom_i = 68;
test_x = squeeze(decomp_hist_nonEc_expanded_fixed0_minE(atom_i,1:21,:))';
test_y = squeeze(decomp_hist_nonEc_expanded_like_fixed0_minE(atom_i,1:21,:))';
[h,p,ci,stats] = ttest2(test_x,test_y);
%atoms_to_plot = [68;107;303;316;17;44;21;80;176;177;202;207;255;328];
%met_names = ["Alanine", "Cysteine", "Serine", "Threonine", "A-Glu", "AGN", "A-Ser", "Asp", "GluP", "Gly", "Ile", "Leu", "NH4", "Val"];
atoms_to_plot = [68;107;303;316]; %;17;44;21;80;176;177;202;207;255;328];
met_names = ["Alanine", "Cysteine", "Serine", "Threonine"]; %, "A-Glu", "AGN", "A-Ser", "Asp", "GluP", "Gly", "Ile", "Leu", "NH4", "Val"];
for i=1:length(atoms_to_plot)
    atom_i = atoms_to_plot(i);
    figure(); 
    for sample_i = 1:num_samples
        plot(0:6:120,decomp_hist_nonEc_expanded_fixed0_minE(atom_i,1:21,sample_i), 'LineWidth',2, 'Color', 'red');
        if sample_i == 1
            hold on; 
        end
        plot(0:6:120,decomp_hist_nonEc_expanded_like_fixed0_minE(atom_i,1:21,sample_i), 'LineWidth',2, 'Color', 'blue');         
    end
    hold off;
    xlabel("Simulation time [min]");
    ylabel("relative 15N enrichment");
    ylim([0,1]);
    xlim([0,120]);
    leg1=legend("with flux in"+newline+"Ala:hydroxypyruvate AT (mito)","without flux in "+newline+"Ala:hydroxypyruvate AT (mito)");
    set(leg1,'Box','off');
    set(gca,'fontsize', 32);
    set(gcf, 'color', [1,1,1]);    
end

% Knock-out analysis
ko_prots_b = "prot_"+string(ecModel_batch.enzymes);

[optBiomass1_b, ko_prots1_b] = testAAKnockOuts(ecModel_batch,ko_prots_b);
[optBiomass2_b, ko_prots2_b] = testAAKnockOuts(ecModel_expanded,ko_prots_b);
[optBiomass4_b, ko_prots4_b] = test2AAKnockOuts(ecModel_batch,ko_prots_b);
[optBiomass5_b, ko_prots5_b] = test2AAKnockOuts(ecModel_expanded,ko_prots_b);
improved12_b = find(optBiomass2_b > (optBiomass1_b + 1e-5));
decreased23_b = find(optBiomass2_b < (optBiomass1_b - 1e-5));
improved45_b = find(optBiomass5_b > (optBiomass4_b + 1e-5));
decreased45_b = find(optBiomass5_b < (optBiomass4_b - 1e-5));

% no single gene knockout affected the optimal growth rate, obtained by
% flux balance analysis:
disp("Number of single knock outs affecting optimal growth: "+string(sum(optBiomass1_b<0.999*opt_bio)));

% 33% of double knockouts inhibited growth with and without inclusion of
% the side reactions.
disp("Number of total double knock outs: "+string(length(ko_prots4_b)));
disp("Number of double knock outs inhibiting growth w/o side reactions: "+ ...
    string(sum(optBiomass4_b<0.000001*opt_bio)));
disp("Of those, inhibiting growth also with side reactions: "+ ...
    (sum(optBiomass5_b<0.000001*opt_bio_expanded & optBiomass4_b<0.000001*opt_bio)));

% prot_Q8RXU4 is Threonine aldolase 1
disp("Total double KO with Threonine aldolase 1: "+string(sum(contains(ko_prots4_b,"prot_Q8RXU4"))));
disp("Blocked w/o side reactions: "+string(sum(optBiomass4_b(contains(ko_prots4_b,"prot_Q8RXU4"))<0.000001*opt_bio)));
disp("Blocked also with side reactions: "+ ...
        string(sum(optBiomass5_b(contains(ko_prots5_b,"prot_Q8RXU4"))<0.000001*opt_bio_expanded ...
                    & optBiomass4_b(contains(ko_prots4_b,"prot_Q8RXU4"))<0.000001*opt_bio)));
disp("Optimal biomass prodcution reduced by 4.9%, no reduction with side reactions:" + ...
    string(sum(optBiomass5_b(contains(ko_prots5_b,"prot_Q8RXU4"))>0.999*opt_bio_expanded ...
                & optBiomass4_b(contains(ko_prots4_b,"prot_Q8RXU4"))<0.951*opt_bio)));

optBiomass5_max = max(optBiomass5_b);
optBiomass4_max = max(optBiomass4_b);
percentChangeBiomass5 = round(100*(optBiomass5_b./optBiomass5_max - 1),2);
percentChangeBiomass4 = round(100*(optBiomass4_b./optBiomass4_max - 1),2);
ko_prot_1 = extractAfter(extractBefore(ko_prots5_b,"___"),"prot_");
ko_prot_2 = extractAfter(extractAfter(ko_prots5_b,"___"),"prot_");
resultTable = table(ko_prot_1, ko_prot_2, percentChangeBiomass4, percentChangeBiomass5, ...
'VariableNames', ["Knock out protein 1"; "Knock out protein 2";
    "percentage change of biomass w/o new activities"; "percentage change of biomass w/ new activities"]);
writetable(resultTable, "double_knockouts.csv");


