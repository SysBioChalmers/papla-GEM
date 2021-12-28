%This script adjust simulation parameters to evaluate the TAG production
%following eMOMA evaluation as described by Kim et al. (2019)- doi:10.1186/s13068-019-1518-4

clear; clc;
if ~exist([pwd() '/growthPrediction.m']); error(['Make sure that '...
        'your Current Folder is the one containing the simulation file.']); end
cd ../../;  root = [pwd() '/'];
data = [root 'data/'];
code = [root 'code/'];
cd(code)
mkdir([root 'data/results'])
% Load model
model = importModel('../model/papla-GEM.xml');

%% Add TAG exchange reaction
idx = getIndexes(model, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]'}, 'metcomps');
model = addExchangeRxns(model, 'out', idx);
%% Fix NGAM to low value
model = setParam(model,'eq','r_4046',3.3928);

%% Block exchange reactions
% Block various TCA cycle intermediates (see eMOMA paper) except from
% citrate, which is a known overflow metabolite of various oleaginous
% yeasts (was it also measured in HPLC?).
model = setParam(model, 'eq', {'r_1798', 'r_1586', 'r_2056'}, 0); % Fumarate, 2-oxoglutarate, succinate, 
model = setParam(model, 'eq', {'r_1552', 'r_1989', 'r_1815', 'r_1634'}, 0); % Malate, oxaloacetate, glyoxylate, acetate

% Block various other lipids, we know that TAGs accumulate
model = setParam(model, 'eq', 'r_1727', 0); % Decanoate
model = setParam(model, 'eq', 'r_2134', 0); % 13-demethyllanosterol
model = setParam(model, 'eq', 'r_1753', 0); % fecosterol
model = setParam(model, 'eq', 'r_1757', 0); % ergosterol
model = setParam(model, 'eq', 'r_1788', 0); % episterol
model = setParam(model, 'eq', 'r_1915', 0); % Lanosterol
model = setParam(model, 'eq', 'r_2106', 0); % Lanosterol
model = setParam(model, 'eq', 'r_2137', 0); % ergosta-5,7,22,24(28)-tetraen-3beta-ol
model = setParam(model, 'eq', 'r_1994', 0); % Palmitoleate
model = setParam(model, 'eq', 'r_2189', 0); % Oleate
sol   = solveLP(model,1) % Check that model still functions

% Alternatively, block all exchange out fluxes and only allow CO2, water,
% H+, biomass and TAG to be excreted. This ONLY allows TAGs as overflow
% metabolite, which gives a higher TAG production in the MOMA of reference
% condition, but it also gives less opportunity for the production to be
% improved. For instance, maybe in the reference condition there is flux
% towards citrate, but a KO would abolish that flux and thereby increase
% the TAG production. This situation would not occur if citrate secretion
% was not allowed in the reference condition.
%model = setParam(model,'ub',getExchangeRxns(model,'out'),0);
%model = setParam(model,'ub',{'r_1672','r_2100','r_1832','r_2111','EXC_OUT_s_3036'},1000);

%% Get reference flux distribution by FBA under non-limited growth conditions
model       = setParam(model, 'obj',{'r_2111'}, 1); % Growth as objective
model       = setParam(model, 'lb', {'r_1714'}, -5); % Glucose uptake
modelRef    = model;
solRef      = solveLP(modelRef,1);
printFluxes(modelRef, solRef.x);

%% Block nitrogen exchange to mimick nitrogen limiting conditions
modelLim    = setParam(model, 'eq', 'r_1654', 0);
solLim      = solveLP(modelLim,1); % Confirm that no growth is possible
solMOMA     = MOMA(modelRef,modelLim);
printFluxes(modelRef,solMOMA.x)
solMOMA.x(end)

% RefMuKO, LimProdKO, RefMuOE, LimProdOE
out=zeros(numel(model.rxns),5);
f=waitbar(0,'Running eMOMA...');
for i=1:numel(modelRef.rxns)
    waitbar(i/numel(modelRef.rxns),f,sprintf('Running eMOMA... %.1f%%',(i/numel(modelRef.rxns))*100))
    % Knockout
    modelRefKO  = setParam(modelRef,'eq',i,0);
    solRefKO    = solveLP(modelRefKO,1);
    try
        out(i,1)    = -solRefKO.f; % Growth rate in reference conditions
    catch
        out(i,1)    = 0;
    end
    modelLimKO  = setParam(modelRefKO, 'eq', {'r_1654'}, 0);
    solLimKO    = MOMA(modelRefKO, modelLimKO);
    try
        out(i,2)    = solLimKO.x(end);
    catch
        out(i,2)    = 0;
    end
   
    % Overexpression (forcing 2x higher flux, catch non-functional models)
    noGrowth = true;
    OEfactor = 2.1;
    while (noGrowth==true && OEfactor>1)
        OEfactor = OEfactor - 0.1; % To start at 2, in 10% steps
        modelRefOE  = setParam(modelRef,'eq',i,OEfactor*solRef.x(i));
        solRefOE    = solveLP(modelRefOE,1);
        try
            out(i,3)    = -solRefOE.f; % Growth rate in reference conditions
            noGrowth    = false;
            out(i,5)    = OEfactor;
        catch
            out(i,3)    = 0;
            out(i,5)    = OEfactor;
        end
    end
    modelLimOE  = setParam(modelRefOE, 'eq', {'r_1654'}, 0);
    solLimOE    = MOMA(modelRefOE, modelLimOE);
    try
        out(i,4)    = solLimOE.x(end);
    catch
        out(i,4)    = 0;
    end
end
close(f);

%% Filter to growth rate and TAG exchange from reference condition
refMu = -solRef.f; % Reference growth rate
refEx = solMOMA.x(end); % Reference TAG production in N-restriction

filtRes = out;
filtRes(:,[1,3]) = filtRes(:,[1,3])/refMu; % Normalize growth rates to reference
filtRes(:,[2,4]) = filtRes(:,[2,4])/refEx; % Normalize TAG production to reference

% Only keep those cases where N-restriction resulted in TAG production that
% was at least 5% higher than in the non-mutated strain in N-restriction,
% and the growth rate when N-exchange was allowed was at least 90% of the
% growth rate of the non-mutated strain when allowing N-exchange. Keep
% those reactions where the above is true for either the knockout or the
% overexpression (or both).
prodImpr = 1.05; % Minimum 5% production increase during N-restriction
growRed = 0.90 % Minimum 90% of non-mutated growth rate when N-uptake is allowed
exUp = find((filtRes(:,2)>prodImpr & filtRes(:,1)>growRed) | (filtRes(:,4)>prodImpr & filtRes(:,3)>growRed));
exUp = [exUp, filtRes(exUp,:)];
exUp(:,7) = max(exUp(:,[3,5]),[],2); % Add extra column: maximum production for each reaction
[~,I]=sort(exUp(:,7),1,'descend'); % Sort by maximum production
exUp = exUp(I,:); % Sort by maximum production

exUp = [modelRef.rxns(exUp(:,1)),modelRef.rxnNames(exUp(:,1)),num2cell(exUp)];
% Columns of table now refer to:
% 1: Reaction identifier
% 2: Reaction name
% 3: Position in output vector (not really useful anymore, since we already
%    extracted the reaction identifiers and names
% 4: Growth rate when this reaction is knocked out.
% 5: TAG production when this reaction is knocked out and no N-uptake.
% 6: Growth rate when this reaction is overexpressed.
% 7: TAG production when this reaction is overexpressed and no N-uptake.
% 8: Fraction by which the reaction is overexpressed. Two-fold by default,
%    but could be lower if two-fold did not allow growth.
% 9: Maximum of row 5 and 7, to sort the most promising reactions.

fid = fopen([data 'results/eMOMA_noSterolExch.txt'],'w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',["rxnID" "rxnName" "idx" ...
    "GR_KO" "EX_KO" "GR_OE" "EX_OE" "OEfactor" "EXmax"]);
for j=1:length(exUp)
    fprintf(fid, '%s\t%s\t%i\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%i\t%5.4f\n',exUp{j,:});
end
fclose(fid);