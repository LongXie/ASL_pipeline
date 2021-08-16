function skullstrip_SPM8_session...
    (session_dir, brainmask_gz, infile_gz, out_prefix, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 5
    error('Not enough input argument!');
end

% check input extension
[in_dir, infile_filename, ext] = fileparts(infile_gz);
out_dir = in_dir;
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(infile_filename);
    if ~strcmp(ext, '.nii')
        error('Input file`s extension must be .nii.gz');
    end
else
    error('Input file`s extension must be .nii.gz');
end

% check mask extension
[anat_dir, brainmask_filename, ext] = fileparts(brainmask_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(brainmask_filename);
    if ~strcmp(ext, '.nii')
        error('Mask file`s extension must be .nii.gz');
    end
else
    error('Mask file`s extension must be .nii.gz');
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

%% check and unzip tissue probability maps
%anat_dir = fullfile(session_dir, 'MPRAGE', 'segmentation');

% tissue probability mapes
%brainmask_gz = fullfile(anat_dir, 'BrainMask.nii.gz');

% check if exist
if exist(brainmask_gz, 'file')
    % unzip file for SPM usage
    % brainmask
    %[~, brainmask_filename] = fileparts(brainmask_gz);
    brainmask = fullfile(anat_dir, brainmask_filename);
    system(sprintf('gunzip -c %s > %s', brainmask_gz, brainmask));
else
    % if not exist, report error
    msg = sprintf('ERROR: brain mask do not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end

%% Run 
savecurpath = pwd;
for r = 1:nruns
    
    %% initialization
    fprintf('Clean CBF for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
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
    
    %% Perform skull striping
    P = spm_select('FPlist', in_dir, ['^', infile_filename, '$']);
    BM = spm_select('FPlist', anat_dir , ['^', brainmask_filename, '$']);
    matlabbatch=cell(1,1);
    matlabbatch{1,1}.spm.util.imcalc.input = {BM, P};
    matlabbatch{1,1}.spm.util.imcalc.output = [out_prefix, infile_filename];
    matlabbatch{1,1}.spm.util.imcalc.outdir = {out_dir};
    matlabbatch{1,1}.spm.util.imcalc.expression = strcat('i1.*i2');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);

    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    system('gzip *.nii');
    
end

% remove the unzip anat files
cd(anat_dir)
delete(brainmask_filename);

% back to original dir
cd(savecurpath);
