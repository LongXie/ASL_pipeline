function OutlierCleaning_SCRUBdenoising_test1_session...
    (session_dir, anat_dir, infile_gz, brainmask_gz, out_folder, out_filename_gz, out_stat, LOGTXT)

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

% check output extension
[~, out_filename, ext] = fileparts(out_filename_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(out_filename);
    if ~strcmp(ext, '.nii')
        error('Output file`s extension must be .nii.gz');
    end
else
    error('Output file`s extension must be .nii.gz');
end

% check output data file
[~, ~, ext] = fileparts(out_stat);
if ~strcmp(ext, '.mat')
    error('Output data file`s extension must be .mat');
end

% check mask extension
[~, brainmask_filename, ext] = fileparts(brainmask_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(brainmask_filename);
    if ~strcmp(ext, '.nii')
        error('Mask file`s extension must be .nii.gz');
    end
else
    error('Mask file`s extension must be .nii.gz');
end

%% parameters
%thresh = 0.95;

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
gm_prob_gz = spm_select('FPList', anat_dir, ['^c1\w*.*nii.gz']);
wm_prob_gz = spm_select('FPList', anat_dir, ['^c2\w*.*nii.gz']);
csf_prob_gz = spm_select('FPList', anat_dir, ['^c3\w*.*nii.gz']);
%brainmask_gz = fullfile(anat_dir, 'BrainMask.nii.gz');

% check if exist
if exist(gm_prob_gz, 'file') && exist(wm_prob_gz, 'file') ...
    && exist(csf_prob_gz, 'file')
    % unzip file for SPM usage
    % gm
    [~, gm_filename] = fileparts(gm_prob_gz);
    gm_prob = fullfile(anat_dir, gm_filename);
    system(sprintf('gunzip -c %s > %s', gm_prob_gz, gm_prob));
    % wm
    [~, wm_filename] = fileparts(wm_prob_gz);
    wm_prob = fullfile(anat_dir, wm_filename);
    system(sprintf('gunzip -c %s > %s', wm_prob_gz, wm_prob));
    % csf
    [~, csf_filename] = fileparts(csf_prob_gz);
    csf_prob = fullfile(anat_dir, csf_filename);
    system(sprintf('gunzip -c %s > %s', csf_prob_gz, csf_prob));
    % brain mask
    brainmask = fullfile(anat_dir, brainmask_filename);
    system(sprintf('gunzip -c %s > %s', brainmask_gz, brainmask));
else
    % if not exist, report error
    msg = sprintf('ERROR: Tissue probability maps or brain mask do not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end

%% Run 
savecurpath = pwd;
for r = 1:nruns
    
    %% initialization
    fprintf('Clean CBF using SCRUB denoising for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % generate output folder
    out_dir = fullfile(run_dir, out_folder);
    if isdir(out_dir)
        rmdir(out_dir,'s');
    end
    mkdir(out_dir);
    
    %% check and unzip CBF time series
    cd(run_dir);
    
    % input file
    %infile_gz  = fullfile('perf', 'cbf_0_scrRAWPCASL_4D.nii.gz');
    
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
    
    % load tissue maps
    gmloc = spm_select('FPlist', anat_dir, ['^', gm_filename, '$']);
    wmloc = spm_select('FPlist', anat_dir, ['^', wm_filename, '$']);
    csfloc = spm_select('FPlist', anat_dir, ['^', csf_filename, '$']);
    
    % SCORE
    %[recon, oidx] = RobustBayesianASLdenoising(cbfloc, gmloc, wmloc, csfloc, brainmask, thresh);
    [~,recon,~,oidx] = SCRUBdenoising_test1(cbfloc, gmloc, wmloc, csfloc, brainmask, 'huber', 1,1,0);
    
    %recon(isnan(recon(:))) = 0;
    
    % get statistics
    stat = [];
    stat.RejectRateTotal = mean(oidx>=1);
    stat.RejectRateS1 = mean(oidx==1);
    stat.RejectRateS2 = mean(oidx==2);
    stat.TotalPairs = length(oidx);
    stat.RemainPairs = sum(oidx == 0);
    stat.odix = oidx;
    stat.names = {'RejectRateTotal', ...
                  'TotalPairs', ...
                  'RemainPairs'};
    stat.measures = [stat.RejectRateTotal, ...
                     stat.TotalPairs, ...
                     stat.RemainPairs];

    % save the reconstructed image
    vo = spm_vol(cbfloc);
    vo = vo(1);
    %out_filename = 'cleaned_meanCBF.nii';
    vo.fname = fullfile(out_dir, out_filename);
    spm_write_vol(vo, recon);

    % save statistics
    stat_fn = fullfile(out_dir, out_stat);
    save(stat_fn, 'stat')

    %% Output status
    fprintf('....Rejection rate is %1.2f. %1.0f out of %1.0f remain.).\n', ...
        stat.RejectRateTotal, stat.RemainPairs, stat.TotalPairs);

    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    system('gzip -f *.nii');
    
end

% remove the unzip anat files
cd(anat_dir)
delete(gm_filename, wm_filename, csf_filename, brainmask);

% back to original dir
cd(savecurpath);


