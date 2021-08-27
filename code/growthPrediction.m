clear; clc;
if ~exist([pwd() '/growthPrediction.m']); error(['Make sure that '...
        'your Current Folder is the one containing the simulation file.']); end
cd ../;  root = [pwd() '/'];
data = [root 'data/'];
code = [root 'code/'];
cd(code)
% Load model
model = importModel('../model/papla-GEM.xml');
%% This script compares experimentally measured with model predicted growth rates
model       = setParam(model,'eq',{'r_1718','r_1714'},0);
model       = setParam(model,'obj','r_2111',1);

% Cultivations were in complex medium containing amino acids
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
model = setParam(model, 'lb', aminoacidRxns, -0.1);

% Load file
fid         = fopen('../data/biomass/bioreactor_growth.csv');
fluxData    = textscan(fid,'%f32 %f32 %s','Delimiter',',','HeaderLines',1);
growth      = fluxData{1};
rate        = fluxData{2};
source      = fluxData{3};
fclose(fid);

clear out
for i = 1:length(growth)
    if strcmp(source(i),'glucose')
        modelTmp = setParam(model,'lb','r_1714',-rate(i));
    elseif strcmp(source(i),'xylose')
        modelTmp = setParam(model,'lb','r_1718',-rate(i));
    end
    sol=solveLP(modelTmp);
    out(i) = -sol.f;
end
out=transpose(out);

% Calculate Pearson R, not useful, as deviation from x=y is more relevant.
R = corr(growth,out);
% Calculate RMSE instead
RMSE = sqrt(mean((growth-out).^2)); 
RMSEglc = sqrt(mean((growth(1:end-1)-out(1:end-1)).^2)); 

% Color according to carbon source
[~,~,ic] = unique(source);
cols = [55,126,184; 228,26,28; 77,175,74; 152,78,163; 255,127,0; 255,255,51; 166,86,40; 247,129,191; 153,153,153];
cols = cols/255;
for i=1:2
    plot(growth(ic == i), out(ic == i), 'o', 'LineWidth', 4,...
        'Color', cols(i,:));
    hold on;
end
plot([0 0.5],[0 0.5],':k')
hold off
title('Model-predicted growth rate')
xlabel('Measured growth, h^-^1')
ylabel('Predicted growth, h^-^1')
legend(unique(source),'Location','eastoutside')%,'NumColumns',2)
text(0.05,0.5,['RMSE: ' num2str(RMSE)],'HorizontalAlignment','left','FontSize',12)
text(0.05,0.45,['RMSE (glucose only): ' num2str(RMSEglc)],'HorizontalAlignment','left','FontSize',12)
set(gca,'FontSize',12) % Creates an axes and sets its FontSize to 18
print([root 'data/results/growthPrediction.pdf'],'-dpdf')

