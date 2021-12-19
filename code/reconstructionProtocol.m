%% PROTOCOL FOR THE RECONSTRUCTION OF A Papiliotrema laurentii GENOME SCALE MODEL USING THE RAVEN TOOLBOX

% This protocol was prepared based on the reconstruction protocols and scripts developed for
% Rhodotorula toruloides and Hansenula polymorpha, both available on GitHub:
% https://github.com/SysBioChalmers/rhto-GEM and https://github.com/SysBioChalmers/hanpo-GEM.

% Define paths for model reconstruction.
clear; clc;
if ~exist([pwd() '/reconstructionProtocol.m']); error(['Make sure that '...
        'your Current Folder is the one containing the reconstructionProtocol file.']); end
cd ../;  root = [pwd() '/'];
data = [root 'data/'];
code = [root 'code/'];
cd(code)

%% IMPORT TEMPLATE MODEL
% Load R.toruloides template GEM using importModel()
% Source: https://github.com/SysBioChalmers/rhto-GEM
modelRhto = importModel([data 'templateModels/rhto.xml'], true);
modelRhto.id = 'rhto';

% Correct duplicate metabolite name in rhto-GEM v1.3.0.
metIdx=getIndexes(modelRhto,'N-[(R)-4-phosphonopantothenoyl]-L-cysteine','metnames');
modelRhto.metNames(metIdx(1))={'tmp'};
modelRhto=replaceMets(modelRhto,'tmp',modelRhto.metNames(metIdx(2)),true);

% Yarrowia 
modelYli        = importModel([data 'templateModels/iYali.xml'],true);
modelYli.id     = 'yli';
[modelYli.grRules, modelYli.rxnGeneMat]=standardizeGrRules(modelYli);

% Correct duplicate metabolite name in iYali v4.1.2
metIdx=getIndexes(modelYli,{'N-[(R)-4-phosphonopantothenoyl]-L-cysteine',...
    '6-O-{2-O-[(2-aminoethyl)phosphoryl]-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol',...
    '6-O-{alpha-D-mannosyl-(1-6)-2-O-[(2-aminoethyl)phosphoryl]-alpha-D-mannosyl-(1-4)-alpha-D-glucosaminyl}-O-acyl-1-phosphatidyl-1D-myo-inositol'},...
    'metnames');
for i=1:numel(metIdx)
    modelYli.metNames(metIdx{i}(1))={'tmp'};
    modelYli=replaceMets(modelYli,'tmp',modelYli.metNames(metIdx{i}(2)),true);
end

% Change reversibility fields to match boundaries. Prevents problems with
% MENECO.
for i=1:length(modelYli.rxns)
	if modelYli.lb(i)==0 && not(modelYli.ub(i)==0)
        modelYli.rev(i)=0;
    elseif modelYli.lb(i)<=0 && modelYli.ub(i)>=0
        modelYli.rev(i)=1;
	end
end
% Confirm that the model is functional, set objective to growth
modelYli        = setParam(modelYli,'obj','r_4041',1);
solveLP(modelYli)

modelYli.rxns = regexprep(modelYli.rxns,'y00','r_');
modelYli = removeReactions(modelYli,contains(modelYli.rxns,'y10'),true,true,true);

% Create the 'scrap' folder, in which all files created during the reconstruction,
% which are not the final model, will be stored.
% Create the excel file of template for easy inspection during reconstruction steps.
mkdir([root 'scrap'])
%exportToExcelFormat(modelRhto, [root 'scrap/rhto.xlsx']);
%exportToExcelFormat(modelYli, [root 'scrap/yli.xlsx']);
% Prepare for gap-filling with meneco
exportModel(modelRhto, [root 'data/meneco/rhto.xml']);

% Save MATLAB environment.
save([root '/scrap/modelTemplate.mat'], 'model*');
% load([root 'scrap/modelTemplate.mat'],'model*')

%% GENERATE MODEL FROM HOMOLOGY 
% PAPLA protein FASTA IDs were replaced by shorter IDs and are now following the 
% naming pattern PAPLA_01234.

% DETERMINE HOMOLOGY by BLAST 
% BLAST the whole-genome protein FASTA of P. laurentii against the
% R. toruloides protein FASTA.
blastRhto = getBlast('papla',[data '/genomes/Papla_protein.fasta'], ...
            'rhto',[data '/genomes/rhto_np11.faa']);

% BLAST the whole-genome protein FASTA of P. laurentii against the
% Y. lipolytica protein FASTA.        
blastYli=getBlast('papla',[data '/genomes/Papla_protein.fasta'], ...
            'yli',[data '/genomes/yli_clib122.faa']);
%  
% Save intermediate files in 'scrap' folder.
save([root '/scrap/blastStruct.mat'],'blast*');
%load([root '/scrap/blastStruct.mat'])

% Use the blast results to generate the first draft model.
model=getModelFromHomology(modelRhto,blastRhto,'papla',{},1,false,10^-20,150,35);

% Add some reactions as based on homology with Yarrowia lipolytica
modelYli    = getModelFromHomology(modelYli,blastYli,'papla',{},1,false,10^-20,150,35);

% Discard reactions that were already in draft rhto-GEM
modelYli    = removeReactions(modelYli,contains(modelYli.rxns,model.rxns),true,true,true);

% Focus on reactions derived from yeast-GEM
tmp         = removeReactions(modelYli,cellfun(@isempty,regexp(modelYli.rxns,'r_\d{4}$')),true,true,true);

% How the Yarrowia model was constructed, there is a set of new metabolites
% that were introduced by simplifying lipid metabolism. Discard these
% metabolites and associated reactions.
tmp         = removeMets(tmp,contains(tmp.mets,'m'),false,true,true,true);
tmp         = removeReactions(tmp,~ismember(tmp.rxns,modelRhto.rxns));
model       = addRxnsGenesMets(model,modelRhto,tmp.rxns,tmp.grRules,'Identified from homology to Yarrowia lipolytica',2);

% Add Yarrowia specific reactions
tmp                 = removeReactions(modelYli,~contains(modelYli.rxns,'y'),true,true,true);
% Replace old identifiers with new format
oldIdx              = find(contains(tmp.mets,'m'));
tmp.mets(oldIdx)    = generateNewIds(model,'mets','m_',numel(oldIdx));

modelComb           = mergeModels({model,tmp});
model               = contractModel(modelComb);

save([root '/scrap/homologyModel.mat'],'model');
%load([root 'scrap/homologyModel.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
% 
% To inspect the first draft model:
exportToExcelFormat(model,[root '/scrap/r1_paplaGEM.xlsx']);
clear blastRhto blastYli
% 
%% DEFINE BIOMASS COMPOSITION
% Use the biomass pseudoreactions from rhto-GEM as template to modify.

% Find all reactions with 'pseudreaction' in reactio name in rhto-GEM, and
% add these to the draft model.
biomassRxns = modelRhto.rxns(endsWith(modelRhto.rxnNames, 'pseudoreaction'));
model = addRxnsGenesMets(model, modelRhto, biomassRxns);

% Add exchange reactions for media components
mediumComps = {'r_1654', 'r_1672', 'r_1808', 'r_1832', 'r_1861', ...
               'r_1992', 'r_2005', 'r_2060', 'r_2100', 'r_2111'};
model = addRxnsGenesMets(model, modelRhto, mediumComps);

% Add all exchange rxns
% These were not gene annotated, and therefore not added in draft.
% Might not require all exchange rxns, but easier to remove unconnected ones later.
rxns    = getExchangeRxns(modelRhto);
model   = addRxnsGenesMets(model,modelRhto,rxns,false,'Modelling reaction',1);

% Add all non-gene annotated transport reactions
noGeneIdx   = find(cellfun(@isempty,modelRhto.grRules)); % Which rxns have no genes
rxnIdx      = find(getTransportRxns(modelRhto));
rxnIdx      = intersect(rxnIdx,noGeneIdx); % Keep the ones without gene anotation
rxns        = modelRhto.rxns(rxnIdx); % Obtain reaction IDs
model       = addRxnsGenesMets(model,modelRhto,rxns,false,'Modeling reaction required for intercellular transport, gene unknown',1);

% For the lipid curation and gapfilling 
model = addRxnsGenesMets(model, modelRhto,{'r_4062', 'r_4064', 'r_4046'});
model = setParam(model, 'ub', {'r_4062', 'r_4064', 'r_4046'}, 1000);
model = setParam(model, 'lb', {'r_4062', 'r_4064', 'r_4046'}, 0);

% Load biomass information
fid           = fopen([data 'biomass/biomassCuration.csv']);
loadedData    = textscan(fid, '%q %q %q %f','delimiter', ',', 'HeaderLines', 1);
fclose(fid);
BM.name       = loadedData{1};    BM.mets     = loadedData{2};
BM.pseudorxn  = loadedData{3};    BM.coeff    = loadedData{4};

% Nucleotides (DNA)
% Find out which rows contain the relevant information
indexes = find(contains(BM.pseudorxn, 'DNA'));
% Define new stoichiometries
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
% Change reaction
model = changeRxns(model, 'r_4050', equations, 1);

% Ribonucleotides (RNA)
indexes = find(contains(BM.pseudorxn, 'RNA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4049', equations, 1);

% Amino acids (protein)
indexes = find(contains(BM.pseudorxn, 'AA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4047', equations, 1);

% Carbohydrates
indexes = find(contains(BM.pseudorxn, 'carbohydrate'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4048', equations, 1);

% Lipid backbones
indexes = find(contains(BM.pseudorxn, 'backbone'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4063', equations, 1);

% Lipid chains
indexes = find(contains(BM.pseudorxn, 'chain'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4065', equations, 1);

save([root 'scrap/biomass.mat'],'model')
% load([root 'scrap/biomass.mat'])

clear ans biomassRxns BM equations fid indexes loadedData mediumComps noGeneIdx rxnIdx rxns

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to inspect:
exportToExcelFormat(model, [root 'scrap/r2_paplaGEM.xlsx']);

%% CURATION OF LIPID REACTIONS
% P. laurentii has unique fatty acid and lipid class compositions. SLIMEr explicitly models each
% lipid moiety, with unique chain distribution, but to reduce complexity we will only include a
% subset of possible chain distributions. To do this, files with templates reactions will be
% modified to match the desired chain distributions. That distribution followed the S. cerevisiae
% lipid chain distribution (Yeast8 - ) based on the simmilarities about the acil chains found on
% both species. Here, R. toruloides distribution was not used once its model present some
% specificities adjusted from experimental data about the presence of specific chains in certain
% positions, which are not confirmed in P. laurentii.
fid         = fopen([data '/reconstruction/lipidTemplates.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', ...
    '\t', 'HeaderLines', 1);
fclose(fid);

% Reorganize the content so that it can be used by the addLipidReactions
% function.
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.grRules     = loadedData{4};
template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+4}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');

% Remove reactions that match a lipid template reaction (ignoring acyl-chains)
toRemove    = regexprep(template.rxns, 'CHAIN.*', '');
toRemove    = find(startsWith(model.rxnNames, toRemove));
model       = removeReactions(model, toRemove);

% Now use the templates to add the relevant reactions to the model. If a
% reaction already existed in the R. toruloides template model, then it
% will use the same reaction identifier.
cd([code 'lipidMetabolism'])
model = addLipidReactions(template, model, modelRhto);

fid         = fopen([data '/reconstruction/lipidTransport.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', '\t');
fclose(fid);

clear template
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.chains = {};
for k = 1:numCols-2; template.chains(:,k) = loadedData{k+3}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');
toRemove    = regexprep(template.rxns,'CHAIN.*','');
toRemove    = find(startsWith(model.rxnNames,toRemove));
model       = removeReactions(model,toRemove);

model = addLipidReactions(template, model, modelRhto);

% Apply SLIME reactions
% Remove any SLIME reactions that might exist in the draft model.
model = removeReactions(model,contains(model.rxnNames,'SLIME rxn'));

% Load SLIME template reactions and parse it through the addSLIMEreactions
% function to amend the model.
fid             = fopen([data '/reconstruction/SLIMERtemplates.txt']);
firstLine       = fgets(fid);
numCols         = numel(strfind(firstLine,char(9))); % number of \t
loadedData      = textscan(fid,['%q %q %f' repmat(' %q',[1,numCols-2])],'delimiter','\t');
fclose(fid);
template.metName    = loadedData{1};	template.bbID   = loadedData{2};
template.bbMW       = loadedData{3};    template.comps  = loadedData{4};
template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+4}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');

model = addSLIMEreactions(template, model, modelRhto);
cd(code)
% Remove unused genes and metabolites
model = deleteUnusedGenes(model);
model = removeMets(model,all(model.S == 0,2),false,true,true,true);

save([root 'scrap/lipids.mat'],'model')
% load([root 'scrap/lipids.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to inspect:
exportToExcelFormat(model, [root 'scrap/r3_paplaGEM.xlsx']);
exportModel(model,[root 'data/meneco/r3_paplaGEM.xml'])
clear ans fid firstLine k loadedData numCols template toRemove

%% PERFORM GAP-FILLING
% MENECO
% Find targets: any substrate for the pseudoreactions, us the following
% text to reconstruct menecoTargets.sbml.
rxnIdx  = find(contains(model.rxnNames,'pseudoreaction'));
targets = find(any(model.S(:,rxnIdx)<0,2));
 [model.mets(targets), model.metNames(targets)]
 targetSBML=strcat('<species id="M_',model.mets(targets),...
     '" name="',model.metNames(targets),'"/>');

% Identified by MENECO (see meneco.txt for output file).
% A minimum of 34 reactions are required, with different combinations of 36
% reactions. Not to favour one reaction over the other, as we don't
% have any prove at the moment which one is more likely to be present, we
% will add the union of reactions.
fid         = fopen([data '/meneco/menecoRxns.txt']);
menecoRxns  = textscan(fid,'%s'); fclose(fid);
menecoRxns  = menecoRxns{1};

% If these reactions are present, that means that their respective enzymes
% are present. Any other reaction annotated to the same enzymes should also
% be added.
menecoRxns  = getAllRxnsFromGenes(modelRhto,menecoRxns);
model       = addRxnsGenesMets(model,modelRhto,menecoRxns,true,'Identified by MENECO to produce biomass components',1);

model   = setParam(model,'obj','r_2111',1);
sol     = solveLP(model,1)
printFluxes(model, sol.x)

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

save([root 'scrap/meneco.mat'],'model')
%load([root 'scrap/meneco.mat'])

% Export to Excel format for easy inspection
%exportToExcelFormat(model,[root '/scrap/r4_paplaGEM.xlsx']);

% RAVEN fillGaps using glucose as carbon source
% Use biomass production as obj func for gapfilling
model = setParam(model, 'obj', 'r_2111', 1);

% Set biomass production to arbitrary low flux, to force gap-filling to
% produce biomass.
model = setParam(model, 'lb', 'r_2111', 0.01);
% Block 3',5'-AMP exchange, as it will otherwise be suggested by fillGaps
model = setParam(model, 'ub', 'r_1641', 0);

% From the Rhto model, remove all exchange reactions (the
% necessary ones we already added, don't want to add new ones)
modelRhto2 = removeReactions(modelRhto, getExchangeRxns(modelRhto, 'both'), true, true, true);
biomassRxns = modelRhto2.rxns(endsWith(modelRhto2.rxnNames, 'pseudoreaction'));
modelRhto2 = removeReactions(modelRhto2, biomassRxns, true, true, true);

% Run fillGaps function
[~, ~, addedRxns, model] = fillGaps(model, modelRhto2, false, true);

% Verify that model can now grow
sol = solveLP(model, 1)
printFluxes(model, sol.x)

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

save([root 'scrap/glucosefillgaps.mat'],'model')
%load([root 'scrap/glucosefillgaps.mat'])

% Export to Excel format for easy inspection
%exportToExcelFormat(model,[root '/scrap/r4_glcfillgaps_paplaGEM.xlsx']);

% Remove duplicate reaction
model       = removeReactions(model, 'r_4046_rhto');

% RAVEN fillGaps using xylose as carbon source
% remove glucose and set xylose as carbon source
model = setParam(model, 'lb', {'r_1714','r_1718'},[0,-1]);    

% Verify that model can grow
sol = solveLP(model, 1)

% Run fillGaps function
[~, ~, addedRxns, model] = fillGaps(model, modelRhto2, false, true);

% Verify that model can now grow
sol = solveLP(model, 1)
printFluxes(model, sol.x)

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

save([root 'scrap/xylosefillgaps.mat'],'model')
%load([root 'scrap/xylosefillgaps.mat'])

% Export to Excel format for easy inspection
%exportToExcelFormat(model,[root '/scrap/r4_xylfillgaps_paplaGEM.xlsx']);

% Remove duplicate reaction
model       = removeReactions(model, 'r_4046_rhto');

% return glucose as carbon source
model = setParam(model, 'lb', {'r_1714','r_1718'},[-1,0]);

cd([code 'lipidMetabolism'])
model = scaleLipids(model, 'tails');
cd(code)

% Verify the change of fluxes in the model after lipid scaling
sol = solveLP(model, 1)
printFluxes(model, sol.x)

% Reset lower bound of biomass production and allow 3',5'-AMP exchange
model = setParam(model, 'lb', 'r_2111', 0);
model = setParam(model, 'ub', 'r_1641', 1000);

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
%exportToExcelFormat(model,[root '/scrap/r5_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/r5_paplaGEM.xml'])

save([root '/scrap/gapfilling.mat'],'model');
% load([root 'scrap/gapfilling.mat'])
clear addedRxns ans biomassRxns fid menecoRxns modelRhto2 modelRhto3 rxnIdx sol targets targetSBML

%% MANUAL CURATION
% Include some missing essential reactions, e.g. complex IV and reactions
% related to assimilation of alternative carbon sources. Some of them, which were not
% present in the R. toruloides template model, were imported from S. cerevisiae (Yeast8).
fid         = fopen([data '/reconstruction/manualCuration.txt']);
loadedData  = textscan(fid,'%q %q','delimiter','\t'); fclose(fid);
rxns        = loadedData{1};
grRules     = regexprep(loadedData{2},'***','');

model = addRxnsGenesMets(model,modelRhto,rxns,grRules,'Identified from homology, manual curation',2);

% Manual curation importing reactions from modelSce
fid2         = fopen([data '/reconstruction/manualCuration2.txt']);
loadedData2  = textscan(fid2,'%q %q','delimiter','\t'); fclose(fid2);
rxns2        = loadedData2{1};
grRules2     = regexprep(loadedData2{2},'***','');

modelSce = importModel([data 'templateModels/sce.xml'], true);
modelSce.id = 'sce';

model = addRxnsGenesMets(model,modelSce,rxns2,grRules2,'Identified from homology, manual curation',2);

% Export to Excel format for easy inspection
%exportToExcelFormat(model,[root '/scrap/r6_paplaGEM.xlsx']);

% Modify gene associations of reactions annotationed
% with R. toruloides and Y. lipolytica genes.
rhtoRxns = find(contains(model.grRules, 'RHTO'));
model.rxnNames(rhtoRxns);
model.rxns(rhtoRxns);
model.grRules(rhtoRxns);

yliRxns = find(contains(model.grRules, 'yli'));
model.rxnNames(yliRxns);
model.rxns(yliRxns);
model.grRules(yliRxns);

fid         = fopen([data '/reconstruction/updateGrRules.txt']);
loadedData  = textscan(fid,'%q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

rxns = loadedData{1};
grRules = regexprep(loadedData{2},'***','');
model = changeGrRules(model,rxns,grRules);

% Remove any reactions containing 18:3 chain that might exist in the draft model.
model = removeReactions(model,contains(model.rxnNames,'18:3'));

% Clean up model to remove rhto remnants in subSystems and unused metabolites
run([code 'curation/cleanupModel']);

% Set exchange reactions to alternative carbon sources to reversible
model = setParam(model,'rev',{'r_1634','r_1706','r_1710','r_1711','r_1714','r_1715','r_1718','r_1808','r_1878','r_2058'},1);

model = deleteUnusedGenes(model);
% Save workspace
save([root 'scrap/cleanup.mat'],'model')
% load([root 'scrap/cleanup.mat'])

%Export to inspect:
%exportToExcelFormat(model, [root '/scrap/r7_paplaGEM.xlsx']);

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

%% Set NGAM and fit GAEC
% Due to limited bioreactor data, unable to reliably fit NGAM. Instead, use the
% NGAM that was fitted in rhto-GEM to >17 datapoints = 3.3928
model = setParam(model,'lb','r_4046',3.3928);
%nutrient uptake reactions to simulate complex medium conditions

% Fit GAEC based on bioreactor cultivation data gathered in this study
cd([root 'code/curation'])
[model, GAEC] = fitGAEC(model);
disp(['GAEC is set to: ' num2str(GAEC)])
% model = setParam(model,'lb', aminoacidRxns, 0);

%% Add model information.
model.annotation.defaultLB    = -1000; % Default lower bound
model.annotation.defaultUB    = +1000; % Default upper bound
model.annotation.taxonomy     = 'taxonomy/460523';
model.annotation.givenName    = 'Rafaela'; 'Mauricio';
model.annotation.familyName   = 'Ventorim'; 'Ferreira';
model.annotation.email        = 'rafaela.ventorim@ufv.br'; 'mauricio.moura@ufv.br';
model.annotation.organization = 'Universidade Federal de Vicosa';
model.annotation.note         = 'Genome-scale model of Papiliotrema laurentii UFV-1';
model.id                      = 'paplaGEM';
model.name                    = 'Papiliotrema laurentii-GEM';
model.annotation.sourceUrl    = 'https://github.com/SysBioChalmers/papla-GEM';

% Remove unnecessary fields
model = rmfield(model,{'metFrom','rxnFrom','geneFrom','geneShortNames'});
model = deleteUnusedGenes(model);

% Save workspace
save([root 'scrap/finalmodel.mat'],'model')
% load([root 'scrap/finalmodel.mat'])

%Export to inspect:
exportToExcelFormat(model,[root '/scrap/papla-GEM.xlsx']);

%% Store model
cd([root 'code'])
newCommit(model)