%% PROTOCOL FOR THE RECONSTRUCTION OF A Papiliotrema laurentii GENOME SCALE MODEL USING THE RAVEN TOOLBOX
%  AUTHORS: Rafaela Z. Ventorim, Mauricio A. M. Ferreira, Eduard J. Kerkhoven, Wendel B. da Silveira

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

% Create the 'scrap' folder, in which all files created during the reconstruction,
% which are not the final model, will be stored.
% Create the excel file of template for easy inspection during reconstruction steps.
mkdir([root 'scrap'])
exportToExcelFormat(modelRhto, [root 'scrap/rhto.xlsx']);

% Save MATLAB environment.
save([root 'scrap/importModel.mat'])
% load([root 'scrap/importModel.mat'])

%% GENERATE MODEL FROM HOMOLOGY 
% PAPLA protein FASTA IDs were replaced by shorter IDs and are now following the 
% naming pattern PAPLA_01234.

% DETERMINE HOMOLOGY by BLAST 
% BLAST the whole-genome protein FASTA of P. laurentii against the
% R. toruloides protein FASTA.
blastRhto = getBlast('papla',[data '/genomes/Papla_protein.fasta'], ...
            'rhto',[data '/genomes/rhto_np11.faa']);
%  
% Save intermediate files in 'scrap' folder.
save([root '/scrap/blastStruct.mat'],'blast*');
% 
% Use the blast results to generate the first draft model.
model=getModelFromHomology(modelRhto,blastRhto,'papla',{},1,false,10^-20,150,35);
% 
save([root '/scrap/model_r1.mat'],'model');
%load([root 'scrap/model_r1.mat'])
% 
disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
% 
% To inspect the first draft model:
exportToExcelFormat(model,[root '/scrap/r1_paplaGEM.xlsx']);
clear blast mediumComps
% 
%% DEFINE BIOMASS COMPOSITION
% Use the biomass pseudoreactions from rhto-GEM as template to modify.
% 
% Find all reactions with 'pseudreaction' in reactio name in rhto-GEM, and
% add these to the draft model.
biomassRxns = modelRhto.rxns(endsWith(modelRhto.rxnNames, 'pseudoreaction'));
model = addRxnsGenesMets(model, modelRhto, biomassRxns);
% 
% Add exchange reactions for media components
mediumComps = {'r_1654', 'r_1672', 'r_1808', 'r_1832', 'r_1861', ...
               'r_1992', 'r_2005', 'r_2060', 'r_2100', 'r_2111'};
model = addRxnsGenesMets(model, modelRhto, mediumComps);
% 
% Add all exchange rxns
% These were not gene annotated, and therefore not added in draft.
% Might not require all exchange rxns, but easier to remove unconnected ones later.
rxns    = getExchangeRxns(modelRhto);
model   = addRxnsGenesMets(model,modelRhto,rxns,false,'Modelling reaction',1);
% 
% Add all non-gene annotated transport reactions
noGeneIdx   = find(cellfun(@isempty,modelRhto.grRules)); % Which rxns have no genes
rxnIdx      = find(getTransportRxns(modelRhto));
rxnIdx      = intersect(rxnIdx,noGeneIdx); % Keep the ones without gene anotation
rxns        = modelRhto.rxns(rxnIdx); % Obtain reaction IDs
model       = addRxnsGenesMets(model,modelRhto,rxns,false,'Modeling reaction required for intercellular transport, gene unknown',1);
% 
% For the lipid curation and gapfilling 
model = addRxnsGenesMets(model, modelRhto,{'r_4062', 'r_4064', 'r_4046'});
model = setParam(model, 'ub', {'r_4062', 'r_4064', 'r_4046'}, 1000);
model = setParam(model, 'lb', {'r_4062', 'r_4064', 'r_4046'}, 0);
% 
% Load biomass information
fid           = fopen([data 'biomass/biomassCuration.csv']);
loadedData    = textscan(fid, '%q %q %q %f','delimiter', ',', 'HeaderLines', 1);
fclose(fid);
BM.name       = loadedData{1};    BM.mets     = loadedData{2};
BM.pseudorxn  = loadedData{3};    BM.coeff    = loadedData{4};
% 
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

clear indexes equations loadedData fid BM biomassRxns ans

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Extport to inspect:
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
for k = 1:length(loadedData)-4; template.chains(:,k) = loadedData{k+4}; end

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
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', ...
    '\t', 'HeaderLines', 1);
fclose(fid);
clear template
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+3}; end

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
clear template
template.metName    = loadedData{1};	template.bbID   = loadedData{2};
template.bbMW       = loadedData{3};    template.comps  = loadedData{4};
template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+4}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');

model = addSLIMEreactions(template, model, modelRhto);
cd(code)

save([root 'scrap/lipids.mat'],'model')
% load([root 'scrap/lipids.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to inspect:
exportToExcelFormat(model, [root 'scrap/r3_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/r3_paplaGEM.xml'])
clear fid loadedData template k toRemove ans

%% PERFORM GAP-FILLING

% RAVEN fillGaps
% Use biomass production as obj func for gapfilling
model = setParam(model, 'obj', 'r_2111', 1);

% Set biomass production to arbitrary low flux, to force gap-filling to
% produce biomass.
model = setParam(model, 'lb', 'r_2111', 0.01);

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

cd([code 'lipidMetabolism'])
model = scaleLipids(model, 'tails');
cd(code)

% Verify the change of fluxes in the model after lipid scaling
sol = solveLP(model, 1)
printFluxes(model, sol.x)

% Back lower bound of biomass production to zero.
model = setParam(model, 'lb', 'r_2111', 0);
model = deleteUnusedGenes(model);

save([root '/scrap/gapfilling.mat'],'model');
% load([root 'scrap/gapfilling.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/r4_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/r4_paplaGEM.xml'])
clear addedRxns rxns sol biomassRxns

% MENECO
% Find targets: any substrate for the pseudoreactions, us the following
% text to reconstruct menecoTargets.sbml.
rxnIdx  = find(contains(model.rxnNames,'pseudoreaction'));
targets = find(any(model.S(:,rxnIdx)<0,2));
 [model.mets(targets), model.metNames(targets)]
 targetSBML=strcat('<species id="M_',model.mets(targets),...
     '" name="',model.metNames(targets),'"/>');

% Identified by MENECO (see meneco.txt for output file).
% A minimum of 13 reactions are required, with different combinations of 19
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

% Back lower bound of biomass production to zero.
model = setParam(model, 'lb', 'r_2111', 0);
model = deleteUnusedGenes(model);

save([root '/scrap/gapfillingfinal.mat'],'model');
% load([root 'scrap/gapfillingfinal.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/r5_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/r5_paplaGEM.xml'])
clear addedRxns rxns sol biomassRxns

%% MANUAL CURATION
% Include some missing essential reactions, e.g. xylulokinase, complex IV and reactions
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
exportToExcelFormat(model,[root '/scrap/r6_paplaGEM.xlsx']);

% Modify gene associations of gap-filled reactions, annotationed
% with R. toruloides genes.
rhtoRxns = find(contains(model.grRules, 'RHTO'));
model.rxnNames(rhtoRxns);
model.grRules(rhtoRxns);

fid         = fopen([data '/reconstruction/updateGrRules.txt']);
loadedData  = textscan(fid,'%q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

rxns = loadedData{1};
grRules = regexprep(loadedData{2},'***','');
model = changeGrRules(model,rxns,grRules);
%     
% Remove any reactions containing 18:3 chain that might exist in the draft model.
model = removeReactions(model,contains(model.rxnNames,'18:3'));

% Remove duplicate reaction
model       = removeReactions(model, 'r_4046_rhto');

% Clean up model to remove rhto remnants in subSystems and unused metabolites
run([code 'curation/cleanupModel']);

% Set exchange reactions to alternative carbon sources to reversible
model = setParam(model,'rev',{'r_1634','r_1706','r_1710','r_1711','r_1714','r_1715','r_1718','r_1808','r_1878','r_2058'},1);

% Set lipid related reactions ICDH, THFS and all fatty-acid-CoA ligases as irreversible
model = setParam(model,'lb',{'r_0659','r_0446'},0);
model = setParam(model,'lb',contains(model.rxnNames,'fatty-acid--CoA ligase'),0);
model = setParam(model,'ub','r_4046',1000);
model = removeUnusedGenes(model);
% Save workspace
save([root 'scrap/cleanup.mat'])
% load([root 'scrap/cleanup.mat'])

%Export to inspect:
exportToExcelFormat(model, [root '/scrap/r7_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/r7_paplaGEM.xml'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%  
%% VALIDATION
% Parameters for simulating growth in the experimental conditions evaluated
%% Exponential phase before before the establishment of the permanent regime (Glucose 5 g/L)
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));
rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
% 
targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));
% 
% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);
% 
% nutrient uptake reactions to simulate complex medium conditions
aminoacidRxns = {'r_1810'; ... % L-glycine
                 'r_1873'; ... % L-alanine
                 'r_1879'; ... % L-arginine
                 'r_1880'; ... % L-asparagine
                 'r_1881'; ... % L-aspartate
                 'r_1883'; ... % L-cysteine
                 'r_1889'; ... % L-glutamate
                 'r_1891'; ... % L-glutamine
                 'r_1893'; ... % L-histidine
                 'r_1897'; ... % L-isoleucine
                 'r_1899'; ... % L-leucine
                 'r_1900'; ... % L-lysine
                 'r_1902'; ... % L-methionine
                 'r_1903'; ... % L-phenylalanine
                 'r_1904'; ... % L-proline
                 'r_1905'; ... % L-serine
                 'r_1911'; ... % L-threonine
                 'r_1912'; ... % L-tryptophan
                 'r_1913'; ... % L-tyrosine
                 'r_1914'};    % L-valine
%              
model = setParam(model, 'lb', aminoacidRxns, -0.01);
% 
% set the carbon source and unlimited O2 for aerobic growth
model = setParam(model, 'lb', {'r_1714'}, -1.741);    % glucose
model = setParam(model, 'ub', {'r_1714'}, 0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
%
% set biomass pseudoreaction as objective
model = setParam(model, 'lb', {'r_2111'}, 0);   % block biomass uptake
%model = setParam(model, 'ub', {'r_2111'}, 1000);
model = setParam(model, 'obj',{'r_2111'}, 1);  
sol = solveLP(model)
printFluxes(model, sol.x, true);
% 
%% Chemostat simmulation
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));
rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
% 
targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));
% 
% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);
% 
% nutrient uptake reactions to simulate complex medium conditions
aminoacidRxns = {'r_1810'; ... % L-glycine
                 'r_1873'; ... % L-alanine
                 'r_1879'; ... % L-arginine
                 'r_1880'; ... % L-asparagine
                 'r_1881'; ... % L-aspartate
                 'r_1883'; ... % L-cysteine
                 'r_1889'; ... % L-glutamate
                 'r_1891'; ... % L-glutamine
                 'r_1893'; ... % L-histidine
                 'r_1897'; ... % L-isoleucine
                 'r_1899'; ... % L-leucine
                 'r_1900'; ... % L-lysine
                 'r_1902'; ... % L-methionine
                 'r_1903'; ... % L-phenylalanine
                 'r_1904'; ... % L-proline
                 'r_1905'; ... % L-serine
                 'r_1911'; ... % L-threonine
                 'r_1912'; ... % L-tryptophan
                 'r_1913'; ... % L-tyrosine
                 'r_1914'};    % L-valine
%              
model = setParam(model, 'lb', aminoacidRxns, -0.01);
% 
% adjust parameters for chemostat growth
model = setParam(model, 'eq', {'r_2111'}, 0.05);  % fix specific growth rate at the dilution rate value
% 
uptake = find(strcmp(model.rxnNames,'D-glucose exchange')); % remove constraints on substrate uptake
model = setParam(model, 'lb', uptake, -Inf);
model = setParam(model, 'ub', uptake, Inf);
% 
% minimize substrate uptake
model = setParam(model, 'obj',{'r_2111'}, 0);
model = setParam(model, 'obj', uptake, 1);
sol = solveLP(model,1)
printFluxes(model, sol.x, true);
% printFluxes(model, sol.x, false);
% 
%% Batch cultivation under Carbon limitation (Glucose 1.5 g/L)
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));
rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
% 
targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));
% 
% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);
% 
% nutrient uptake reactions to simulate complex medium conditions
aminoacidRxns = {'r_1810'; ... % L-glycine
                 'r_1873'; ... % L-alanine
                 'r_1879'; ... % L-arginine
                 'r_1880'; ... % L-asparagine
                 'r_1881'; ... % L-aspartate
                 'r_1883'; ... % L-cysteine
                 'r_1889'; ... % L-glutamate
                 'r_1891'; ... % L-glutamine
                 'r_1893'; ... % L-histidine
                 'r_1897'; ... % L-isoleucine
                 'r_1899'; ... % L-leucine
                 'r_1900'; ... % L-lysine
                 'r_1902'; ... % L-methionine
                 'r_1903'; ... % L-phenylalanine
                 'r_1904'; ... % L-proline
                 'r_1905'; ... % L-serine
                 'r_1911'; ... % L-threonine
                 'r_1912'; ... % L-tryptophan
                 'r_1913'; ... % L-tyrosine
                 'r_1914'};    % L-valine
%              
model = setParam(model, 'lb', aminoacidRxns, -0.01);
% 
% set the carbon source and unlimited O2 for aerobic growth
model = setParam(model, 'lb', {'r_1714'}, -1.816);    % glucose
model = setParam(model, 'ub', {'r_1714'}, 0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
%
% set biomass pseudoreaction as objective
model = setParam(model, 'lb', {'r_2111'}, 0);   % block biomass uptake
%model = setParam(model, 'ub', {'r_2111'}, 1000);
model = setParam(model, 'obj',{'r_2111'}, 1);  
sol = solveLP(model)
printFluxes(model, sol.x, true);
% 
%% Batch cultivation under Nitrogen limitation - Glucose (30 g/L)
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));
rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
% 
targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));
% 
% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);
% 
% nutrient uptake reactions to simulate complex medium conditions
aminoacidRxns = {'r_1810'; ... % L-glycine
                 'r_1873'; ... % L-alanine
                 'r_1879'; ... % L-arginine
                 'r_1880'; ... % L-asparagine
                 'r_1881'; ... % L-aspartate
                 'r_1883'; ... % L-cysteine
                 'r_1889'; ... % L-glutamate
                 'r_1891'; ... % L-glutamine
                 'r_1893'; ... % L-histidine
                 'r_1897'; ... % L-isoleucine
                 'r_1899'; ... % L-leucine
                 'r_1900'; ... % L-lysine
                 'r_1902'; ... % L-methionine
                 'r_1903'; ... % L-phenylalanine
                 'r_1904'; ... % L-proline
                 'r_1905'; ... % L-serine
                 'r_1911'; ... % L-threonine
                 'r_1912'; ... % L-tryptophan
                 'r_1913'; ... % L-tyrosine
                 'r_1914'};    % L-valine
%              
model = setParam(model, 'lb', aminoacidRxns, -0.01);
% 
% set the carbon source and unlimited O2 for aerobic growth
model = setParam(model, 'lb', {'r_1714'}, -4.158);    % glucose
model = setParam(model, 'ub', {'r_1714'}, 0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
%
% set biomass pseudoreaction as objective
model = setParam(model, 'lb', {'r_2111'}, 0);   % block biomass uptake
%model = setParam(model, 'ub', {'r_2111'}, 1000);
model = setParam(model, 'obj',{'r_2111'}, 1);  
sol = solveLP(model)
printFluxes(model, sol.x, true);
% 
%% Batch cultivation under Nitrogen limitation - Xylose (30 g/L)
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));
rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
% 
targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));
% 
% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);
% 
% nutrient uptake reactions to simulate complex medium conditions
aminoacidRxns = {'r_1810'; ... % L-glycine
                 'r_1873'; ... % L-alanine
                 'r_1879'; ... % L-arginine
                 'r_1880'; ... % L-asparagine
                 'r_1881'; ... % L-aspartate
                 'r_1883'; ... % L-cysteine
                 'r_1889'; ... % L-glutamate
                 'r_1891'; ... % L-glutamine
                 'r_1893'; ... % L-histidine
                 'r_1897'; ... % L-isoleucine
                 'r_1899'; ... % L-leucine
                 'r_1900'; ... % L-lysine
                 'r_1902'; ... % L-methionine
                 'r_1903'; ... % L-phenylalanine
                 'r_1904'; ... % L-proline
                 'r_1905'; ... % L-serine
                 'r_1911'; ... % L-threonine
                 'r_1912'; ... % L-tryptophan
                 'r_1913'; ... % L-tyrosine
                 'r_1914'};    % L-valine
%              
model = setParam(model, 'lb', aminoacidRxns, -0.01);
% 
% set the carbon source and unlimited O2 for aerobic growth
model = setParam(model, 'lb', {'r_1718'}, -7.77);    % xylose
model = setParam(model, 'ub', {'r_1718'}, 0);
model = setParam(model, 'lb', {'r_1714'}, 0);    % glucose
model = setParam(model, 'ub', {'r_1714'}, 0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
%
% set biomass pseudoreaction as objective
model = setParam(model, 'lb', {'r_2111'}, 0);   % block biomass uptake
%model = setParam(model, 'ub', {'r_2111'}, 1000);
model = setParam(model, 'obj',{'r_2111'}, 1);  
sol = solveLP(model)
printFluxes(model, sol.x, true);
% 
%% SIMULATIONS
% Parameters for simulating batch growth in minimal medium with different carbon sources
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));
% 
rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
% 
targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));
% 
% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);
%
% set the carbon source and unlimited O2 for aerobic growth
% Here, we adjuste lb for each carbon source one at a time,
% while the others remain set to zero and perform the simulation
model = setParam(model, 'lb', {'r_1714'}, -10);    % glucose
model = setParam(model, 'ub', {'r_1714'}, 0);
model = setParam(model, 'lb', {'r_1634'}, 0);    % acetate
model = setParam(model, 'ub', {'r_1634'}, 0);
model = setParam(model, 'lb', {'r_1710'}, 0);    % galactose
model = setParam(model, 'ub', {'r_1710'}, 0);
model = setParam(model, 'lb', {'r_2058'}, 0);    % sucrose
model = setParam(model, 'ub', {'r_2058'}, 0);
model = setParam(model, 'lb', {'r_1718'}, 0);    % xylose
model = setParam(model, 'ub', {'r_1718'}, 0);
model = setParam(model, 'lb', {'r_1711'}, 0);    % galacturonate
model = setParam(model, 'ub', {'r_1711'}, 0);
model = setParam(model, 'lb', {'r_1808'}, 0);    % glycerol
model = setParam(model, 'ub', {'r_1808'}, 0);
model = setParam(model, 'lb', {'r_1878'}, 0);    % L-arabinose
model = setParam(model, 'ub', {'r_1878'}, 0);
model = setParam(model, 'lb', {'r_1706'}, 0);    % D-arabinose
model = setParam(model, 'ub', {'r_1706'}, 0);
model = setParam(model, 'lb', {'r_1715'}, 0);    % mannose
model = setParam(model, 'ub', {'r_1715'}, 0);
model = setParam(model, 'lb', {'r_1546'}, 0);    % (R)-lactate
model = setParam(model, 'ub', {'r_1546'}, 0);
model = setParam(model, 'lb', {'r_1551'}, 0);    % (S)-lactate
model = setParam(model, 'ub', {'r_1551'}, 0);
% 
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
%
% set biomass pseudoreaction as objective
model = setParam(model, 'lb', {'r_2111'}, 0);   % block biomass uptake
%model = setParam(model, 'ub', {'r_2111'}, 1000);
model = setParam(model, 'obj',{'r_2111'}, 1);   
sol = solveLP(model,1)
printFluxes(model, sol.x, true);
%printFluxes(model, sol.x, false);
%

%%  OLEAGINOUS FBA
%  Add exchange reactions for triglyceride (16:0/18:1/18:1-TAG).
idx = getIndexes(model, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]', ...
    'triglyceride (1-16:0, 2-18:1, 3-18:1)[lp]'}, 'metscomps');
% Add exchange reactions for products
rxnsToAdd.rxns          = 'exch_TAG';
rxnsToAdd.mets          = model.mets(idx);
rxnsToAdd.stoichCoeffs  = {[-1, -1]}; 
rxnsToAdd.lb            = 0;
model = addRxns(model,rxnsToAdd);

model = setParam(model,'obj','exch_TAG',1);
model = setParam(model,'lb','r_4046',0);
sol=solveLP(model,1)
printFluxes(model,sol.x)

% Make sure that COBRA Toolbox version >3 is installed
initCobraToolbox(false)

model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [-1, 0]);
Glc.grRatio = singleRxnDeletion(model,'FBA');

modelXyl    = setParam(model, 'eq', {'r_1714', 'r_1718'}, [0, -1]);
Xyl.grRatio = singleRxnDeletion(modelXyl,'FBA');


idx = find(Glc.grRatio < 0.90 | Xyl.grRatio < 0.90);

out = [num2cell(Glc.grRatio(idx)*100) num2cell(Xyl.grRatio(idx)*100) ...
    model.rxns(idx) model.rxnNames(idx) ...
    constructEquations(model,idx)];

mkdir([root 'data/results'])
fid = fopen([data 'results/oleaginous.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',["glucose" "xylose" "rxns" "rxnName" "eqn"]);
for j=1:length(idx)
    fprintf(fid,'%d\t%d\t%s\t%s\t%s\n',out{j,:});
end
fclose(fid);
%
%% This script predicts metabolic engineering targets for increased
% production of triglycerides.
% Add exchange reactions for linolenate and triglyceride. For the
% triglyceride, we specifically choose 16:0/18:1/18:1-TAG as target.
idx = getIndexes(model, {'linolenate[c]','triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]'}, 'metscomps');
%
% Add exchange reactions for products
model = addExchangeRxns(model, 'out', idx);
% Keep track of ids of exchange reactions
rxn1 = model.rxns(end-1);
rxn2 = model.rxns(end);
%
% Perform FSEOF for linolenic acid and TAG on glucose
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [-1, 0]);
%
targets{1} = FSEOF(model, 'r_2111', rxn1, 10, 0.9);
targets{2} = FSEOF(model, 'r_2111', rxn2, 10, 0.9);
%
% Perform FSEOF for linolenic acid and TAG on xylose
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [0, -1]);
%
targets{3} = FSEOF(model, 'r_2111', rxn1, 10, 0.9);
targets{4} = FSEOF(model, 'r_2111', rxn2, 10, 0.9);
%
% Summarize results in table
geneAssoc = ~cellfun('isempty',model.grRules);
for i=1:size(targets,2)
    target(:,i)=targets{i}.logical;
end
for i=1:size(targets,2)
    slope(:,i)=targets{i}.slope;
    slope(~target(:,i),i)=nan;
end
%
target  = find(sum(target,2) & geneAssoc);
[~,I]=sort(sum(slope(target,:),2,'omitnan'),'descend');
out     = [num2cell(slope(target(I),:)), model.rxnNames(target(I)), model.grRules(target(I))];
%
fid = fopen([data '/results/fseof_TAG_linolenicAcid.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\n',["glu_LA" "glu_TAG" "xyl_LA" "xyl_TAG" ...
    "rxnName" "grRule"]);
for j=1:length(I)
    fprintf(fid,'%d\t%d\t%d\t%d\t%s\t%s\n',out{j,:});
end
fclose(fid);
%
%% Add model information.
model.annotation.defaultLB    = -1000; % Default lower bound
model.annotation.defaultUB    = +1000; % Default upper bound
model.annotation.taxonomy     = 'taxonomy/460523';
model.annotation.givenName    = 'Rafaela'; 'Maurício';
model.annotation.familyName   = 'Ventorim'; 'Ferreira';
model.annotation.email        = 'rafaela.ventorim@ufv.br'; 'mauricio.moura@ufv.br';
model.annotation.organization = 'Universidade Federal de Vicosa';
model.annotation.note         = 'First draft model';
model.id                      = 'paplaGEM';
model.description             = 'Papiliotrema laurentii-GEM';

% Save workspace
save([root 'scrap/finalmodel.mat'],'model')
% load([root 'scrap/finalmodel.mat'])

%Export to inspect:
exportToExcelFormat(model,[root '/scrap/paplaGEM.xlsx']);
exportModel(model,[root 'scrap/paplaGEM.xml'])
