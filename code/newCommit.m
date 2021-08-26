function newCommit(model)

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if strcmp(currentBranch,'master')
    error('ERROR: in master branch. For new releases, use newRelease.m')
end

if ~exist('model')
    %Load model:
    model = importModel('../model/papla-GEM.xml');
end

%Save model
exportForGit(model,'papla-GEM','../model',{'txt', 'xml', 'yml'},false,false);
end