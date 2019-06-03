%% This script generates a first draft of a genome-scale model for P. laurentii (papla)
% Prepare some redirects to folders in the repository
root  = regexprep(pwd(),'(.*)\\[^\\]*\\.*','$1');
scripts = [root '/ComplementaryScripts'];
data    = [root '/ComplementaryData'];

% BLAST against R. toruloides genome. Rhto protein fasta obtained from JGI, modified
% to have genes in format RHTO_00001. Papla protein fasta obtained from UFV.
blastRhto = getBlast('papla',[data '/genome/Papla_protein.fasta'], ...
            'rhto',[data '/genome/rhto_np11.faa']);
% Save intermediate files in 'scrap' folder that is not tracked by Git
mkdir([root '/scrap'])
save([root '/scrap/blastStruct.mat'],'blast*');

% Load R. toruloides model version 1.1.1 (UPDATE WITH NEWER VERSION WHEN
% AVAILABLE!), downloaded from
% https://github.com/SysBioChalmers/rhto-GEM/blob/master/ModelFiles/xml/rhto.xml
modelRhto=importModel([data,'/reconstruction/rhto.xml'],true);
% Check that the model is functional:
sol=solveLP(modelRhto)

% Generate draft model, based on homology:
model=getModelFromHomology(modelRhto,blastRhto,'papla',{},1,false,10^-20,150,35);

% Add meta data
model.annotation.defaultLB    = -1000;
model.annotation.defaultUB    = +1000;
model.annotation.taxonomy     = 'taxonomy/5418'; % NCBI taxonomy ID
model.annotation.givenName    = 'Eduard';
model.annotation.familyName   = 'Kerkhoven';
model.annotation.email        = 'eduardk@chalmers.se';
model.annotation.organization = 'Chalmers University of Technology';
model.annotation.note         = 'Papiliotrema laurentii UFV-1';
model.id                      = 'papla';
model.description             = 'Papiliotrema laurentii-GEM_v0.0.1';

save([root '/scrap/model_r1.mat'],'model');
% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/papla_v0.0.1.xlsx']);
% To prepare files for commit to GitHub, run newCommit function. This
% writes the necessary models files (xml, yml, txt) etc.
cd(scripts); newCommit(model); cd('reconstruction');



