function ANTS_globalCBF_session...
    (session_dir, infile_gz, seg_gz, out_folder, out_filename, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 6
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

% check input extension
[~, seg_filename, ext] = fileparts(seg_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(seg_filename);
    if ~strcmp(ext, '.nii')
        error('Input seg file`s extension must be .nii.gz');
    end
else
    error('Input seg file`s extension must be .nii.gz');
end

% check output data file
[~, ~, ext] = fileparts(out_filename);
if ~strcmp(ext, '.mat')
    error('Output data file`s extension must be .mat');
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
    fprintf('Compute global CBF for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % generate output folder
    out_dir = fullfile(run_dir, out_folder);
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    
    %% check and unzip CBF time series
    cd(run_dir);
    
    % check if file exist
    if exist(infile_gz, 'file') && exist(seg_gz, 'file')
        % unzip file for SPM usage
        infile = fullfile(out_dir, infile_filename);
        system(sprintf('gunzip -c %s > %s', infile_gz, infile));
        segfile = fullfile(out_dir, seg_filename);
        system(sprintf('gunzip -c %s > %s', seg_gz, segfile));
        pause(0.5)
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s or %s does not exist in %s.\n', infile_gz, seg_gz, run_dir);
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end
    
    %% compute global CBF
    % load time series
    cbfloc = spm_select('ExtFPList', out_dir, ['^', infile_filename, '$'], Inf);
    CBFMap = spm_read_vols(spm_vol(cbfloc));
    CBFMapgzero = (CBFMap > 0);
    
    % load GM
    segloc = spm_select('FPlist', out_dir, ['^', seg_filename, '$']);
    seg = spm_read_vols(spm_vol(segloc));
    GMMask = (seg == 2) .* CBFMapgzero;
    WMMask = (seg == 3) .* CBFMapgzero;
    GlobalMask = (seg == 2 | seg == 3) .* CBFMapgzero;
    
    %% compute the average CBFs
    globalCBF = [];
    globalCBF.meanCBF_GM = sum(CBFMap(:) .* GMMask(:)) / sum(GMMask(:));
    globalCBF.meanCBF_WM = sum(CBFMap(:) .* WMMask(:)) / sum(WMMask(:));
    globalCBF.meanCBF_GMWM = sum(CBFMap(:) .* GlobalMask(:)) / sum(GlobalMask(:));

    % save statistics
    stat_fn = fullfile(out_dir, out_filename);
    save(stat_fn, 'globalCBF')

    %% Output status
    fprintf('....GM: %4.2f; WM %4.2f; GMWM: %4.2f.\n', ...
        globalCBF.meanCBF_GM, globalCBF.meanCBF_WM, globalCBF.meanCBF_GMWM);

    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename, seg_filename);
    
end

% back to original dir
cd(savecurpath);
