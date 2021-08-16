function restore_session(session_dir, type, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 3
    error('Not enough input argument!');
end

% check input extension
TYPES = {'func', 'anat', 'both'};
if ~ismember(type, TYPES)
    error('Type must be func, anat or both');
end

%% parameters

%% clean anatomical directory if necessary
if strcmp(type, 'anat') || strcmp(type, 'both')
    
    fprintf('Restore anatomical directory.\n');
    
    % get files
    files = dir(fullfile(session_dir, 'MPRAGE'));
    N_files = length(files);
    
    % delete directories
    for i = 3:N_files
        % delete if it is a directory ant not ants related
        if files(i).isdir == 1 && isempty(files(i).name, 'ants')
            filename = fullfile(session_dir, 'MPRAGE', files(i).name);
            rmdir(filename, 's');
            fprintf('    Deleted %s.\n', files(i).name);
        end
    end
    
end

%% restore functional directory
if strcmp(type, 'func') || strcmp(type, 'both')
    
    % Find bold run directories
    d = listdir(fullfile(session_dir,'*ASL*'),'dirs');
    if isempty(d) %MV
        d = listdir(fullfile(session_dir,'*asl*'),'dirs');
    end
    nruns = length(d);
    
    if nruns == 0
        msg = sprintf('No ASL directories found in %s.\n',session_dir);
        cmd = sprintf('echo "%s" >> %s', msg, LOGTXT);
        system(cmd);
        fprintf(msg);
        return;
    end
    
    for r = 1:nruns
        % initialization
        fprintf('Restore for run %s (%0.0f/%0.0f).\n', ...
            d{r}, r, nruns);
        run_dir = fullfile(session_dir, d{r});
        
        % get files
        files = dir(fullfile(run_dir));
        N_files = length(files);
        
        % delete directories
        for i = 3:N_files
            % delete if it is a directory
            if files(i).isdir == 1 && ...
                   isempty(strfind(files(i).name, 'SR')) && ...
                   isempty(strfind(files(i).name, 'ants'))
                filename = fullfile(run_dir, files(i).name);
                rmdir(filename, 's');
                fprintf('    Deleted %s.\n', files(i).name);
            end
        end
        
    end
    
end
