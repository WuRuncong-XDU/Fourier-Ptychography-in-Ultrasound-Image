function AddFullPath(varargin)
if isempty(varargin)
    varargin{1} = './';
end
len = length(varargin);
for i = 1 : len
    if isfolder(varargin{i})
        addpath(varargin{i});
        addsubpath(varargin{i});
    else
        error([varargin{i} ' is not a valid path.']);
    end
end
% fprintf('Add sub folder done !\n');
end

%% Local
function addsubpath(path)
    subpath = dir(path);
    for i = 3 : length(subpath)
        if subpath(i).isdir
            folder = fullfile(subpath(i).folder, subpath(i).name);
            addpath(folder);
            addsubpath(folder);
        end
    end
end