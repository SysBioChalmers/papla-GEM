%% Simulate TAG production in different cultivations, setting growth and 
% sugar uptake with experimental values
clear; clc;
if ~exist([pwd() '/simulation.m']); error(['Make sure that '...
        'your Current Folder is the one containing the simulation file.']); end
cd ../../;  root = [pwd() '/'];
data = [root 'data/'];
code = [root 'code/'];
cd(code)

% Load model
model = importModel('../model/papla-GEM.xml');

%%  Add exchange reactions for triglyceride (16:0/18:1/18:1-TAG).
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
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);

%% Batch YNB-glucose 5g/L
model = setParam(model,'eq','r_1714',-4.926);
model = setParam(model,'eq','r_2111',0.3527);
sol=solveLP(model);
printFluxes(model, sol.x)

%% Batch Complex- glucose 1.5 g/L
model = setParam(model,'eq','r_1714',-1.82);
model = setParam(model,'eq','r_2111',0.2382);
% %nutrient uptake reactions to simulate complex medium conditions
% aminoacidRxns = {'r_1810'; ... % L-glycine
%                  'r_1873'; ... % L-alanine
%                  'r_1879'; ... % L-arginine
%                  'r_1880'; ... % L-asparagine
%                  'r_1881'; ... % L-aspartate
%                  'r_1883'; ... % L-cysteine
%                  'r_1889'; ... % L-glutamate
%                  'r_1891'; ... % L-glutamine
%                  'r_1893'; ... % L-histidine
%                  'r_1897'; ... % L-isoleucine
%                  'r_1899'; ... % L-leucine
%                  'r_1900'; ... % L-lysine
%                  'r_1902'; ... % L-methionine
%                  'r_1903'; ... % L-phenylalanine
%                  'r_1904'; ... % L-proline
%                  'r_1905'; ... % L-serine
%                  'r_1911'; ... % L-threonine
%                  'r_1912'; ... % L-tryptophan
%                  'r_1913'; ... % L-tyrosine
%                  'r_1914'};    % L-valine              
% model = setParam(model, 'lb', aminoacidRxns, -0.1);

sol=solveLP(model);
printFluxes(model, sol.x)

%% Batch Complex-Glucose Carbon:Nitrogen 100:1
model = setParam(model,'eq','r_1714',-4.158);
model = setParam(model,'eq','r_2111',0.3526);
%nutrient uptake reactions to simulate complex medium conditions
% aminoacidRxns = {'r_1810'; ... % L-glycine
%                  'r_1873'; ... % L-alanine
%                  'r_1879'; ... % L-arginine
%                  'r_1880'; ... % L-asparagine
%                  'r_1881'; ... % L-aspartate
%                  'r_1883'; ... % L-cysteine
%                  'r_1889'; ... % L-glutamate
%                  'r_1891'; ... % L-glutamine
%                  'r_1893'; ... % L-histidine
%                  'r_1897'; ... % L-isoleucine
%                  'r_1899'; ... % L-leucine
%                  'r_1900'; ... % L-lysine
%                  'r_1902'; ... % L-methionine
%                  'r_1903'; ... % L-phenylalanine
%                  'r_1904'; ... % L-proline
%                  'r_1905'; ... % L-serine
%                  'r_1911'; ... % L-threonine
%                  'r_1912'; ... % L-tryptophan
%                  'r_1913'; ... % L-tyrosine
%                  'r_1914'};    % L-valine              
% model = setParam(model, 'lb', aminoacidRxns, -0.1);

sol=solveLP(model);
printFluxes(model, sol.x)

%% Batch Complex-Xylose Carbon:Nitrogen 100:1
model = setParam(model,'eq','r_1714',0);
model = setParam(model,'eq','r_1718',-7.77);
model = setParam(model,'eq','r_2111',0.44);
%nutrient uptake reactions to simulate complex medium conditions
% aminoacidRxns = {'r_1810'; ... % L-glycine
%                  'r_1873'; ... % L-alanine
%                  'r_1879'; ... % L-arginine
%                  'r_1880'; ... % L-asparagine
%                  'r_1881'; ... % L-aspartate
%                  'r_1883'; ... % L-cysteine
%                  'r_1889'; ... % L-glutamate
%                  'r_1891'; ... % L-glutamine
%                  'r_1893'; ... % L-histidine
%                  'r_1897'; ... % L-isoleucine
%                  'r_1899'; ... % L-leucine
%                  'r_1900'; ... % L-lysine
%                  'r_1902'; ... % L-methionine
%                  'r_1903'; ... % L-phenylalanine
%                  'r_1904'; ... % L-proline
%                  'r_1905'; ... % L-serine
%                  'r_1911'; ... % L-threonine
%                  'r_1912'; ... % L-tryptophan
%                  'r_1913'; ... % L-tyrosine
%                  'r_1914'};    % L-valine              
% model = setParam(model, 'lb', aminoacidRxns, -0.1);

sol=solveLP(model);
printFluxes(model, sol.x)
