function model = getEarlierModel(version,unversion)
% getEarlierModel
%   Obtain an earlier model version from the Git repository. If no output
%   is specified, it will keep the file as _earlierModel.xml in the current
%   directory, otherwise it will return the model as loaded by importModel.
%
%   Input:
%   version     string of either 'main' for latest release, or e.g. 
%               '1.0.2' for a specific release.
%   unversion   logical whether version information should be stripped from
%               the model (opt, default false, only applies if output is
%               specified)
%
%   Output:
%   model       model structure from obtained model (opt)
%
%   Usage: model = getEarlierModel(version,unversion)

nargoutchk(0,1)
if nargin<2
    unversion = false;
end

if strcmp(version,'main')
    status=system('git show main:model/papla-GEM.xml > _earlierModel.xml')
elseif regexp(version,'^\d+\.\d+\.\d+$')
    tagpath = ['refs/tags/' version ':model/papla-GEM.xml'];
    status=system(['git show ' tagpath ' > _earlierModel.xml']);
else
    error('''version'' should be either ''main'' or of the format ''1.0.2''.')
end
switch nargout
    case 0
        disp('Earlier model version is stored as ''_earlierModel.xml'' in the current working directory')   
    case 1
        disp('Loading earlier model version.')
        model=importModel('_earlierModel.xml');
        if unversion==true
            model.description='papla-GEM';
        end
        delete('_earlierModel.xml');
end
