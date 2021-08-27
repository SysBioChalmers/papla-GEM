%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAEC = fitGAEC(model)
%
%   Modified from GECKO, under MIT License:
%   https://github.com/SysBioChalmers/SLIMEr/blob/master/models/scaleAbundancesInModel.m
%
% 2021-07-25    Eduard Kerkhoven
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [model,GAEC,err] = fitGAEC(model)

%Load chemostat data:
fid         = fopen('../../data/biomass/bioreactor_growth.csv');
fluxData    = textscan(fid,'%f32 %f32 %s','Delimiter',',','HeaderLines',1);
fluxData    = [fluxData{1} fluxData{2}];
fclose(fid);
fluxData(5,:) = [];

%GAECs to span:
disp('Estimating GAEC:')
GAEC = 1:10:301;

%1st iteration:
GAEC = iteration(model,GAEC,fluxData);

%2nd iteration:
GAEC = iteration(model,GAEC-10:1:GAEC+10,fluxData);

%3rd iteration:
[GAEC, err] = iteration(model,GAEC-1:0.1:GAEC+1,fluxData);

model = setGAEC(model,GAEC);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [GAEC,err] = iteration(model,GAEC,fluxData)

fitting = ones(size(GAEC))*1000;

for i = 1:length(GAEC)
    %Simulate model and calculate fitting:
    mod_data   = abs(simulateChemostat(model,fluxData,GAEC(i)));
    R          = (mod_data - fluxData)./fluxData;
    fitting(i) = sqrt(sum(sum(R.^2)));
    disp(['GAEC = ' num2str(GAEC(i)) ' -> Error = ' num2str(fitting(i))])
end

%Choose best:
[~,best] = min(fitting);

if best == 1 || best == length(GAEC)
    error('GAEC found is sub-optimal: please expand GAEC search bounds.')
else
    GAEC = GAEC(best);
    err = min(fitting);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mod_data = simulateChemostat(model,fluxData,GAEC)

model = setGAEC(model,GAEC);

pos = getIndexes(model',{'r_2111','r_1714'},'rxns');

%Simulate chemostats:
mod_data = zeros(size(fluxData));
for i = 1:4
    %Fix glucose and maximize biomass
    model = setParam(model,'lb',pos(2),-fluxData(i,2));
    model = setParam(model,'ub',pos(2),0);
    model = setParam(model,'obj',pos(1),1);
    sol   = solveLP(model);
    %Store relevant variables:
    mod_data(i,:) = sol.x(pos)';

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = setGAEC(model,GAEC)

xr_pos = getIndexes(model','r_4041','rxns');
for i = 1:length(model.mets)
    S_ix  = model.S(i,xr_pos);
    isGAEC = sum(strcmp({'ATP','ADP','H2O','H+','phosphate'},model.metNames{i})) == 1;
    if S_ix ~= 0 && isGAEC
        model.S(i,xr_pos) = sign(S_ix) * GAEC;
    end
end
end