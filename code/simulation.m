%% SIMULATIONS
clear; clc;
if ~exist([pwd() '/simulation.m']); error(['Make sure that '...
        'your Current Folder is the one containing the simulation file.']); end
cd ../;  root = [pwd() '/'];
data = [root 'data/'];
code = [root 'code/'];
cd(code)
% Load model
model = importModel('../model/papla-GEM.xml');

%% Validate with rates measured in this study

%Load data:
fid         = fopen('../data/biomass/bioreactor_growth.csv');
fluxData    = textscan(fid,'%f32 %f32 %s','Delimiter',',','HeaderLines',1);
cSourceData = fluxData{3};
fluxData    = [fluxData{1} fluxData{2}];
fclose(fid);

% Yeast extract in medium, allow uptake of amino acids
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
modelTmp = setParam(model, 'lb', aminoacidRxns, -0.1);

disp('Growth on glucose:')
for i=1:4
    modelTmp = setParam(modelTmp,'lb','r_1714',-fluxData(i,2));
    sol = solveLP(modelTmp,1);
    fprintf('Growth rate measured: %.4f / simulated: %.4f\n',fluxData(i,1),-sol.f)
end
disp('Growth on xylose:')
modelTmp = setParam(modelTmp,'lb',{'r_1714','r_1718'},[0,-fluxData(5,2)]);
sol = solveLP(modelTmp,1);
fprintf('Growth rate measured: %.4f / simulated: %.4f\n',fluxData(5,1),-sol.f)

%% Simulate different carbon sources
% set the carbon source and unlimited O2 for aerobic growth
cSource = {'glucose','acetate','galactose','sucrose','xylose',...
    'glycerol','L-arabinose','D-arabinose','mannose',...
    '(R)-lactate','(S)-lactate'};
cAtoms = [6,2,6,12,5,3,5,5,6,3,3];
exchRxn = {'r_1714','r_1634','r_1710','r_2058','r_1718',...
    'r_1808','r_1878','r_1706','r_1715','r_1546','r_1551'};
zeroFlx = repmat(0,1,numel(cSource));
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
% set biomass pseudoreaction as objective
model = setParam(model, 'obj',{'r_2111'}, 1);

% Fluxes to track
rxnId = getIndexes(model,{'r_1672','r_1992','r_2111'},'rxns');

% Loop through all carbon sources and print exchange fluxes
for i = 1:numel(cSource)
    fprintf('\nSimulate growth on %s:\n', cSource{i})
    modelTmp = setParam(model,'lb',exchRxn,zeroFlx);
    modelTmp = setParam(modelTmp,'lb',exchRxn{i},-5);
    sol = solveLP(modelTmp);
    cIdx = getIndexes(model,exchRxn{i},'rxns');
    fluxList = [cIdx; rxnId];
    table2(i,:) = abs(sol.x(fluxList)/(cAtoms(i)/6));
end
table2=[cSource', num2cell(table2)];

fid = fopen([data 'results/carbonSources.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',["cSource" "cSource flux" "CO2 flux" "O2 flux" "growth"]);
for j=1:length(cSource)
    fprintf(fid,'%s\t%4.2f\t%4.2f\t%4.2f\t%4.2f\n',table2{j,:});
end
fclose(fid);

%%  OLEAGINOUS FBA
%  Add exchange reactions for triglyceride (16:0/18:1/18:1-TAG).
idx = getIndexes(model, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]', ...
    'triglyceride (1-16:0, 2-18:1, 3-18:1)[lp]'}, 'metcomps');
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
    fprintf(fid,'%3.1f\t%3.1f\t%s\t%s\t%s\n',out{j,:});
end
fclose(fid);

%% This script predicts metabolic engineering targets for increased
% production of triglycerides.
% Add exchange reactions for triglyceride (16:0/18:1/18:1-TAG as target).
idx = getIndexes(model, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]'}, 'metscomps');
%
% Add exchange reactions for products
model = addExchangeRxns(model, 'out', idx);
% Keep track of ids of exchange reactions
rxn = model.rxns(end);
%
% Perform FSEOF for TAG on glucose
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [-1, 0]);
%
targets{1} = FSEOF(model, 'r_2111', rxn, 10, 0.9);
%
% Perform FSEOF for TAG on xylose
model       = setParam(model, 'eq', {'r_1714', 'r_1718'}, [0, -1]);
%
targets{2} = FSEOF(model, 'r_2111', rxn, 10, 0.9);

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
fid = fopen([data '/results/fseof_TAG.tsv'],'w');
fprintf(fid,'%s\t%s\t%s\t%s\n',["glu_TAG" "xyl_TAG" ...
    "rxnName" "grRule"]);
for j=1:length(I)
    fprintf(fid,'%d\t%d\t%s\t%s\n',out{j,:});
end
fclose(fid);