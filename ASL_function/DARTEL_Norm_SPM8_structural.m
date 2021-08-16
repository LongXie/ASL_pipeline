function DARTEL_Norm_SPM8_structural...
    (template, flowfield, infile_gz, ref_gz, out_dir, LOGTXT, sigma)

%% check input and output
if nargin < 6
    error('Not enough input argument!');
end

if nargin < 7
    sigma = 8;
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
[~, ref_filename, ext] = fileparts(ref_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(ref_filename);
    if ~strcmp(ext, '.nii')
        error('Reference file`s extension must be .nii.gz');
    end
else
    error('Reference file`s extension must be .nii.gz');
end

%% check if template and flowfield exist
if ~exist(flowfield, 'file') || ~exist(template, 'file')
    % if not exist, report error
    msg = 'ERROR: Flowfield or template do not exist in %s.\n';
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end

%% initialization
fprintf('Normalize structural .\n');
savecurpath = pwd;

% generate output folder
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% check and unzip CBF time series

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

%% perform normalization
clear matlabbatch
matlabbatch{1}.spm.tools.dartel.mni_norm.template = {template};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.flowfield = {flowfield};
matlabbatch{1}.spm.tools.dartel.mni_norm.data.subj.images = {infile};
matlabbatch{1}.spm.tools.dartel.mni_norm.vox = [2 2 2];
matlabbatch{1}.spm.tools.dartel.mni_norm.bb = [-90 -126 -72
    90   90  108];
matlabbatch{1}.spm.tools.dartel.mni_norm.preserve = 0;
matlabbatch{1}.spm.tools.dartel.mni_norm.fwhm = [sigma sigma sigma];

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

%% remove the original unzip file and zip the outputfile
cd(out_dir);
delete(infile_filename);
norm_map_nii = fullfile(out_dir, ['sw', infile_filename]);
norm_map = fullfile(out_dir, ['sw', infile_filename, '.gz']);
if exist(norm_map, 'file')
    delete(norm_map);
end
system(sprintf('gzip %s', norm_map_nii));

%% reslice the normalized map to the MNI
cmd = sprintf('c3d %s %s -reslice-identity -o %s', ref_gz, norm_map, norm_map);
display(cmd)
system(cmd);

% back to original dir
cd(savecurpath);
