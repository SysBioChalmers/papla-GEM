clear; clc;
if ~exist([pwd() '/growthPrediction.m']); error(['Make sure that '...
        'your Current Folder is the one containing the simulation file.']); end
cd ../;  root = [pwd() '/'];
data = [root 'data/'];
code = [root 'code/'];
cd(code)
mkdir([root 'data/results'])
% Load model
model = importModel('../model/papla-GEM.xml');
%% This script compares experimentally measured with model predicted growth rates
model       = setParam(model,'eq',{'r_1718','r_1714'},0);
model       = setParam(model,'obj','r_2111',1);

% Load file
fid         = fopen('../data/biomass/bioreactor_growth.csv');
fluxData    = textscan(fid,'%f32 %f32 %s','Delimiter',',','HeaderLines',1);
growth      = fluxData{1};
rate        = fluxData{2};
source      = fluxData{3};
fclose(fid);

clear out
for i = 1:length(growth)
    switch source{i}
        case {'YNBglc','glucose'}
            modelTmp = setParam(model,'lb','r_1714',-rate(i));
        case 'xylose'
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
for i=1:3
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
