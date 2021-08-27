function newRelease(bumpType)
% newRelease
%   Prepares a new release of the papla-GEM model, for direct submission to
%   GitHub. This function should be run from the /code directory.
%
%   bumpType    string specifying the type of release, either 'major',
%               'minor' or 'patch'
%
%   Usage: newRelease(bumpType)

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if ~strcmp(currentBranch,'master')
    error('ERROR: not in master')
end

%Bump version number:
oldModel   = load('../model/papla-GEM.mat');
oldVersion = oldModel.model.name;
oldVersion = oldVersion(strfind(oldVersion,'_v')+2:end);
oldVersion = str2double(strsplit(oldVersion,'_'));
newVersion = oldVersion;
switch bumpType
    case 'major'
        newVersion(1) = newVersion(1) + 1;
        newVersion(2) = 0;
        newVersion(3) = 0;
    case 'minor'
        newVersion(2) = newVersion(2) + 1;
        newVersion(3) = 0;
    case 'patch'
        newVersion(3) = newVersion(3) + 1;
    otherwise
        error('ERROR: invalid input. Use "major", "minor" or "patch"')
end
newVersion = num2str(newVersion,'%d.%d.%d');

%Check if history has been updated:
fid     = fopen('../history.md','r');
history = fscanf(fid,'%s');
fclose(fid);
if ~contains(history,['papla-GEMv' newVersion ':'])
    error('ERROR: update history.md first')
end

%Load model:
model = importModel('../model/papla-GEM.xml');

%Include tag and save model:
model.name = ['papla-GEM_v' strrep(newVersion,'.','_')];
nGenes=num2str(numel(model.genes));
nMets=num2str(numel(model.mets));
nRxns=num2str(numel(model.rxns));

%Save model
exportForGit(model,'papla-GEM','../model',{'mat', 'txt', 'xlsx', 'xml', 'yml'},true,false);

%Update version file:
fid = fopen('../version.txt','wt');
fprintf(fid,newVersion);
fclose(fid);

%Update model stats in README.md
newStats = ['$1' datestr(now,'dd-mmm-yyyy') ' | ' newVersion ' | ' nRxns ' | ' nMets ' | ' nGenes ' |'];
searchStats = '^(\| \_Papiliotrema laurentii_ UFV-1 \| )\d{2}-\D{3}-\d{4} \| \d+\.\d+\.\d+ \| \d+ \| \d+ \| \d+ \|';
fOld = fopen('../README.md','rt');
fNew = fopen('../README.new','w+');
while ~feof(fOld)
    str = fgets(fOld);
    fwrite(fNew,regexprep(str,searchStats,newStats));
end
fclose(fNew);
fclose(fOld);
movefile '../README.new' '../README.md' f
end