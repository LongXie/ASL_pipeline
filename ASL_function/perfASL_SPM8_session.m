function perfASL_SPM8_session(session_dir, infile_gz, out_folder, LOGTXT,labeltime,delaytime)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 4
    error('Not enough input argument!');
elseif nargin < 5
    % default value for Dave's PMC MCI data
    labeltime = 1.51;
    delaytime = 1.5;
elseif nargin < 6
    delaytime = 1.5;
end

% check input extension
[~, infile_filename, ext] = fileparts(infile_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(infile_filename);
    if ~strcmp(ext, '.nii')
        error('Input file`s extension must be .nii.gz');
    end
else
    error('Input file`s extension must be .nii.gz');
end


%% parameters


%% Find ASL run directories
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

%% Run 
savecurpath = pwd;
for r = 1:nruns
    
    %% initialization
    fprintf('Quantify CBF for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % generate output folder
    out_dir = fullfile(run_dir, out_folder);
    if isdir(out_dir)
        rmdir(out_dir,'s');
    end
    mkdir(out_dir);
    
    %% check and unzip input files
    cd(run_dir);
    
    % input file
    %infile_gz  = fullfile('smooth', 'scrRAWPCASL.nii.gz');
    
    % check if file exist
    if exist(infile_gz, 'file')
        % unzip file for SPM usage
        %[~, infile_filename] = fileparts(infile_gz);
        infile = fullfile(out_dir, infile_filename);
        system(sprintf('gunzip -c %s > %s', infile_gz, infile));
        pause(0.5)
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s does not exist in %s.\n', infile_gz, run_dir);
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end
    
    %% perform smooth
    % select files
    P = spm_select('ExtFPList', out_dir, infile_filename, Inf);
    
    asl_perf_subtract(P, 0, 0, ...
        1,      [1 1 1 0 0 1 0 1 1], 0.5,     1,      0.85, 1,...
        labeltime, delaytime, 55, 18, [],[],[]);     %Last [] is for maskimg;
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    system('gzip *.nii');
    
end
cd(savecurpath);
