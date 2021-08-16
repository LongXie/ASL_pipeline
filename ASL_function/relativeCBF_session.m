function relativeCBF_session...
    (session_dir, infile_gz, gbCBF, out_folder, out_prefix, type, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 7
    error('Not enough input argument!');
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

% check gbCBF file
[~, ~, ext] = fileparts(gbCBF);
if ~strcmp(ext, '.mat')
    error('Input data file`s extension must be .mat');
end

%% parameters
if ~strcmp(type, 'GM') && ~strcmp(type, 'WM') && ~strcmp(type, 'GMWM')
    error('Type %s is not supported.', type)
end

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
    fprintf('Compute relative CBF map for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});

    % output dir
    out_dir = fullfile(run_dir, out_folder);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    
    %% check and unzip CBF time series
    cd(run_dir);
    
    % input file
    %in_dir = 'clean_SCORE';
    %infile_gz  = fullfile(in_dir, 'ss_cleaned_meanCBF.nii.gz');
    %[~, infile_filename] = fileparts(infile_gz);
    %gbCBF = fullfile(in_dir, 'globalCBF.mat');
    %out_dir = in_dir;
    out_file = fullfile(out_dir, [out_prefix, infile_filename]);
    if exist(out_file, 'file')
        delete(out_file);
    end
    
    % check if file exist
    if exist(infile_gz, 'file') && exist(gbCBF, 'file')
        % unzip file for SPM usage
        infile = fullfile(out_dir, infile_filename);
        system(sprintf('gunzip -c %s > %s', infile_gz, infile));
        pause(0.5)
        load(gbCBF);
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s or %s does not exist in %s.\n', infile_gz, gbCBF, run_dir);
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end
    
    %% compute relative cbf
    % load cbf
    file_func = spm_select('FPList', out_dir, ['^', infile_filename, '$']);
    file_func = deblank(file_func(1,:));
    DataVol = spm_vol(file_func);
    Data = spm_read_vols(DataVol);
    
    % load global cbf
    if strcmp(type, 'GM')
        globalCBF = globalCBF.meanCBF_GM;
    elseif strcmp(type, 'WM')
        globalCBF = globalCBF.meanCBF_WM;
    elseif strcmp(type, 'GMWM')
        globalCBF = globalCBF.meanCBF_GMWM;
    else
        error('Type %s is not supported.', type)
    end

    % check if globalCBF is 0 or NA
    if globalCBF == 0 || isnan(globalCBF)
        msg = sprintf('ERROR: globalCBF can not be %s.\n', num2str(globalCBF));
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end
    
    % compute relative cbf
    Data = Data./ globalCBF;
    DataVol.fname = out_file;
    spm_write_vol(DataVol, Data);
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    system('gzip -f *.nii');
    
end

% back to original dir
cd(savecurpath);
