function newsegment_SPM8_session(session_dir, MPRAGE_gz, out_folder, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 4
    error('Not enough input argument!');
end

% check input extension
[anat_dir, infile_filename, ext] = fileparts(MPRAGE_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(infile_filename);
    if ~strcmp(ext, '.nii')
        error('Input file`s extension must be .nii.gz');
    end
else
    error('Input file`s extension must be .nii.gz');
end

%% check if anatomical scan exist or not and unzip

% check if this subject has already been processed 
out_dir = fullfile(anat_dir, out_folder);
c1_gz = fullfile(out_dir, 'c1MPRAGE.nii.gz');
c2_gz = fullfile(out_dir, 'c2MPRAGE.nii.gz');
c3_gz = fullfile(out_dir, 'c3MPRAGE.nii.gz');

if ( ~isempty(dir(c1_gz)) && ~isempty(dir(c2_gz)) && ~isempty(dir(c3_gz)))
    msg = 'Segmented files already found for this subject! Skipping segmentation SPM8\n';
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    fprintf(msg);
    return
end

% generate out_dir
if exist(out_dir, 'dir')
   rmdir(out_dir, 's');
end
mkdir(out_dir);

% unzip mprage
%MPRAGE_gz = fullfile(anat_dir, 'MPRAGE.nii.gz');
if ~exist(MPRAGE_gz, 'file')
    % if not exist, report error
    msg = sprintf('ERROR: MPRAGE does not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
else
    % unzip anatomical scan
    [~, anat_filename] = fileparts(MPRAGE_gz);
    [~, anat_main] = fileparts(anat_filename);
    MPRAGE = fullfile(out_dir, anat_filename);
    system(sprintf('gunzip -c %s > %s', MPRAGE_gz, MPRAGE));
end 

%% perform new segment
matlabbatch = cell(1,1);
spm('defaults', 'FMRI');
matlabbatch{1}.spm.tools.preproc8.channel.vols = {MPRAGE};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {fullfile(spm('Dir'), 'toolbox','Seg', 'TPM.nii,1')};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 1]; % for DARTEL
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {fullfile(spm('Dir'), 'toolbox','Seg', 'TPM.nii,2')};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 1]; % for DARTEL
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {fullfile(spm('Dir'), 'toolbox','Seg', 'TPM.nii,3')};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 1]; % for DARTEL
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {fullfile(spm('Dir'), 'toolbox','Seg', 'TPM.nii,4')};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {fullfile(spm('Dir'), 'toolbox','Seg', 'TPM.nii,5')};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {fullfile(spm('Dir'), 'toolbox','Seg', 'TPM.nii,6')};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0];

% run job
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);


%% zip generated images and remove unziped images
delete(MPRAGE);
cd(out_dir)
system('gzip *.nii');

