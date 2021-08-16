function PVC_session...
    (session_dir, gm_prob_gz, wm_prob_gz, infile_gz, out_folder, out_prefix, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 5
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


%% parameters
GMProbThresh = 0.5;

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
%gm_prob_gz = spm_select('FPList', anat_dir, ['^c1\w*.*nii.gz']);
%wm_prob_gz = spm_select('FPList', anat_dir, ['^c2\w*.*nii.gz']);
%brainmask_gz = spm_select('FPList', anat_dir, ['GMWMMask.nii.gz']);

% check if exist
if exist(gm_prob_gz, 'file') && exist(wm_prob_gz, 'file')
    % unzip file for SPM usage
    % gm
    [gm_anat_dir, gm_filename] = fileparts(gm_prob_gz);
    gm_prob = fullfile(gm_anat_dir, gm_filename);
    system(sprintf('gunzip -c %s > %s', gm_prob_gz, gm_prob));
    % wm
    [wm_anat_dir, wm_filename] = fileparts(wm_prob_gz);
    wm_prob = fullfile(wm_anat_dir, wm_filename);
    system(sprintf('gunzip -c %s > %s', wm_prob_gz, wm_prob));
else
    % if not exist, report error
    msg = sprintf('ERROR: Tissue probability maps do not exist in %s.\n', run_dir);
    system(sprintf('echo %s >> %s', msg, LOGTXT));
    error(msg);
end

%% Run 
savecurpath = pwd;
for r = 1:nruns
    
    %% initialization
    fprintf('Performing PVC for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % generate output folder
    out_dir = fullfile(run_dir, out_folder);
    if exist(out_dir, 'dir')
        rmdir(out_dir,'s');
    end
    mkdir(out_dir);
    
    
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
    
    %% load images
    % load time series
    cbfloc = spm_select('ExtFPList', out_dir, ['^', infile_filename, '$'], Inf);
    cbfvol = spm_vol(cbfloc);
    CBFMap_ORIG = spm_read_vols(cbfvol);
    CBFMap = CBFMap_ORIG;%  zeros(size(CBFMap_ORIG));
    
    % load GM
    gmloc = spm_select('FPlist', gm_anat_dir, ['^', gm_filename, '$']);
    GMProbMask = spm_read_vols(spm_vol(gmloc));
    idx = GMProbMask >= GMProbThresh;
    
    % load WM
    wmloc = spm_select('FPlist', wm_anat_dir, ['^', wm_filename, '$']);
    WMProbMask = spm_read_vols(spm_vol(wmloc));
    
    %% PVE correction
    CBFMap(idx)=CBFMap_ORIG(idx)./(GMProbMask(idx)+0.4*WMProbMask(idx));

    %% save file
    out_file = fullfile(out_dir, [out_prefix, infile_filename]);
    cbfvol.fname = out_file;
    spm_write_vol(cbfvol, CBFMap);
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    system('gzip *.nii');
    
end

% remove the unzip anat files
cd(gm_anat_dir)
delete(gm_filename);
cd(wm_anat_dir)
delete(wm_filename);

% back to original dir
cd(savecurpath);
