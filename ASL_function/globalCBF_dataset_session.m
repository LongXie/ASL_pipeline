function globalCBF_dataset_session...
    (session_dir, anat_dir, infile_gz, out_folder, out_filename, LOGTXT)

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

% check output data file
[~, ~, ext] = fileparts(out_filename);
if ~strcmp(ext, '.mat')
    error('Output data file`s extension must be .mat');
end

%% parameters
GMProbThresh = 0.8;
WMProbThresh = 0.9;

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

%% check and unzip tissue probability maps
%anat_dir = fullfile(session_dir, 'MPRAGE', 'segmentation');

% tissue probability mapes
gm_prob_gz = fullfile(anat_dir, 'c1rMPRAGE.nii.gz');
wm_prob_gz = fullfile(anat_dir, 'c2rMPRAGE.nii.gz');
brainmask_gz = fullfile(anat_dir, 'GMWMMask.nii.gz');

% check if exist
if exist(gm_prob_gz, 'file') && exist(wm_prob_gz, 'file') ...
    && exist(brainmask_gz, 'file')
    % unzip file for SPM usage
    % gm
    [~, gm_filename] = fileparts(gm_prob_gz);
    gm_prob = fullfile(anat_dir, gm_filename);
    system(sprintf('gunzip -c %s > %s', gm_prob_gz, gm_prob));
    % wm
    [~, wm_filename] = fileparts(wm_prob_gz);
    wm_prob = fullfile(anat_dir, wm_filename);
    system(sprintf('gunzip -c %s > %s', wm_prob_gz, wm_prob));
    % brainmask
    [~, brainmask_filename] = fileparts(brainmask_gz);
    brainmask = fullfile(anat_dir, brainmask_filename);
    system(sprintf('gunzip -c %s > %s', brainmask_gz, brainmask));
else
    % if not exist, report error
    msg = sprintf('ERROR: Tissue probability maps or brain mask do not exist in %s.\n', run_dir);
    system(sprintf('echo %s >> %s', msg, LOGTXT));
    error(msg);
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
    
    %% check and unzip CBF time series
    cd(run_dir);
    
    % input file
    %infile_gz  = fullfile('clean_SCORE', 'ss_cleaned_meanCBF.nii.gz');
    
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
    
    %% perform cleaning
    % load time series
    cbfloc = spm_select('ExtFPList', out_dir, ['^', infile_filename, '$'], Inf);
    CBFMap = spm_read_vols(spm_vol(cbfloc));
    CBFMapgzero = (CBFMap > 0);
    
    % load GM
    gmloc = spm_select('FPlist', anat_dir, ['^', gm_filename, '$']);
    GMProbMask = spm_read_vols(spm_vol(gmloc));
    GMMask = (GMProbMask > GMProbThresh).*CBFMapgzero;
    
    % load WM
    wmloc = spm_select('FPlist', anat_dir, ['^', wm_filename, '$']);
    WMProbMask = spm_read_vols(spm_vol(wmloc));
    WMMask = (WMProbMask > WMProbThresh) .* CBFMapgzero;
    
    % load GM WM mask
    bmloc = spm_select('FPlist', anat_dir , ['^', brainmask_filename, '$']);
    GlobalProbMask = spm_read_vols(spm_vol(bmloc));
    GlobalMask = (GlobalProbMask > 0) .* CBFMapgzero;
    
    %% erode wm mask
    SE = ones(2);
    WMMask_eroded = zeros(size(WMMask));
    for qq = 1:size(WMMask,3)
        WMMask_eroded(:,:,qq) = imerode(WMMask(:,:,qq), SE);
    end
    
    %% compute the average CBFs
    globalCBF = [];
    meanCBF_GM = sum(CBFMap(:) .* GMMask(:)) / sum(GMMask(:));
    meanCBF_WM = sum(CBFMap(:) .* WMMask_eroded(:)) / sum(WMMask_eroded(:));
    meanCBF_GMWM = sum(CBFMap(:) .* GlobalMask(:)) / sum(GlobalMask(:));
    globalCBF.names = {'meanCBF_GM', 'meanCBF_WM', 'meanCBF_GMWM'};
    globalCBF.measures = [meanCBF_GM, meanCBF_WM, meanCBF_GMWM];
    
    % save statistics
    stat_fn = fullfile(out_dir, out_filename);
    save(stat_fn, 'globalCBF')

    %% Output status
    fprintf('....GM: %4.2f; WM %4.2f; GMWM: %4.2f.\n', ...
        meanCBF_GM, meanCBF_WM, meanCBF_GMWM);

    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    
end

% remove the unzip anat files
cd(anat_dir)
delete(gm_filename, wm_filename, brainmask_filename);

% back to original dir
cd(savecurpath);
