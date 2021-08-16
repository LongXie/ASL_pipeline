function newsegment_lowresolution_ss_SPM8_session...
    (session_dir, MPRAGE_gz, brainmask_gz, func_gz, out_folder, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 6
    error('Not enough input argument!');
end

% check input MPRAGE extension
[anat_dir, anat_filename, ext] = fileparts(MPRAGE_gz);
if strcmp(ext, '.gz')
    [~, anat_main, ext] = fileparts(anat_filename);
    if ~strcmp(ext, '.nii')
        error('Anatomical scan`s extension must be .nii.gz');
    end
else
    error('Anatomical scan`s extension must be .nii.gz');
end

% check input func extension
[~, func_filename, ext] = fileparts(func_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(func_filename);
    if ~strcmp(ext, '.nii')
        error('Functional scan`s extension must be .nii.gz');
    end
else
    error('Functional scan`s extension must be .nii.gz');
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
global defaults;
spm_get_defaults;
defs = defaults.realign;

% Flags to pass to routine to create resliced images
% (spm_reslice)
resFlags = struct(...
    'interp', 1,...                       % trilinear interpolation
    'wrap', defs.write.wrap,...           % wrapping info (ignore...)
    'mask', defs.write.mask,...           % masking (see spm_reslice)
    'which',1,...                         % write reslice time series for later use
    'mean',0); 

thresh_brain = 0.5;
thresh_gmwm = 0.1;
savecurpath = pwd;
SmoothFWHM = 6;

%% Find bold run directories
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

%% check if anatomical scan exist or not and unzip
% generate outdir
%anat_dir = fullfile(session_dir, 'MPRAGE');
%MPRAGE_gz = fullfile(anat_dir, 'MPRAGE.nii.gz');

out_dir = fullfile(anat_dir, out_folder);
if exist(out_dir, 'dir')
   rmdir(out_dir, 's');
end
mkdir(out_dir);

% unzip mprage
if ~exist(MPRAGE_gz, 'file')
    % if not exist, report error
    msg = sprintf('ERROR: MPRAGE does not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
else
    % unzip anatomical scan
    %[~, anat_filename] = fileparts(MPRAGE_gz);
    %[~, anat_main] = fileparts(anat_filename);
    MPRAGE = fullfile(out_dir, anat_filename);
    system(sprintf('gunzip -c %s > %s', MPRAGE_gz, MPRAGE));
end 

% load anatomical image

PG = spm_select('FPList', out_dir, anat_filename);
PG = PG(1,:);

% smooth anatomical image before subsampleing
%Q = fullfile(out_dir,['s' spm_str_manip(PG,'dt')]);
%spm_smooth(PG,Q,SmoothFWHM);

%% check and unzip brain mask maps

% check if exist
if exist(brainmask_gz, 'file')
    % unzip file for SPM usage
    % brainmask
    %[~, brainmask_filename] = fileparts(brainmask_gz);
    brainmask = fullfile(out_dir, brainmask_filename);
    system(sprintf('gunzip -c %s > %s', brainmask_gz, brainmask));
else
    % if not exist, report error
    msg = sprintf('ERROR: brain mask do not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end

%% load mean functional image
% unzip and load mean functional image
%func_pathname = fullfile();
%func_gz = fullfile(func_pathname, 'cmeanRAWPCASL.nii.gz');

% find rest ASL
idx = find(cellfun(@isempty, strfind(d, 'Rest'))==0);
if isempty(idx)
    idx = 1;
end
func_gz = fullfile(session_dir, d{idx}, func_gz);

if exist(func_gz, 'file') 
    % unzip file for SPM usage
    [func_pathname, func_filename] = fileparts(func_gz);
    func = fullfile(func_pathname, func_filename);
    system(sprintf('gunzip -c %s > %s', func_gz, func));
    pause(0.5)
else
    % if not exist, report error
    msg = sprintf('ERROR: %s or %s does not exist in %s.\n', func_gz, run_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end

PF = spm_select('FPList', func_pathname, ['^', func_filename, '$']);
PF = PF(1, :);

%% reslice anatomical image
clear PO
PO = strvcat(PF, PG);
spm_reslice(PO, resFlags);

%% Perform skull striping
file_ranat = spm_select('FPlist', out_dir, ['^r', anat_filename]);
matlabbatch=cell(1,1);
matlabbatch{1,1}.spm.util.imcalc.input = {file_ranat, brainmask};
matlabbatch{1,1}.spm.util.imcalc.output = ['ss_r', anat_filename];
matlabbatch{1,1}.spm.util.imcalc.outdir = {out_dir};
matlabbatch{1,1}.spm.util.imcalc.expression = strcat('i1.*i2');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

%% new segment seg
file_to_segment = spm_select('FPlist', out_dir, ['^ss_r', anat_filename]);

matlabbatch = cell(1,1);
spm('defaults', 'FMRI');
matlabbatch{1}.spm.tools.preproc8.channel.vols = {file_to_segment};
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

fn_GM_prob = spm_select('FPlist', out_dir, ['^c1ss_r', anat_filename]);
fn_WM_prob = spm_select('FPlist', out_dir, ['^c2ss_r', anat_filename]);
fn_CSF_prob = spm_select('FPlist', out_dir, ['^c3ss_r', anat_filename]);

%% generate brain masks
cd(out_dir)
% 
% % strip skull using a strick threshold
% cmd0=['c3d ', deblank(file_ranat),' -replace nan 0 -o ', deblank(file_ranat)];
% cmd1=['bet ', deblank(file_ranat),' ', ...
%            ['ss_r', anat_main, '_0d8.nii.gz'], ...
%            ' -f 0.8 -m'];
% cmd2=['bet ', deblank(file_ranat),' ', ...
%            ['ss_r', anat_main, '_0d4.nii.gz'], ...
%            ' -f 0.4 -m'];
% cmd3=['gunzip ', ['ss_r', anat_main, '_0d8_mask.nii.gz']];
% cmd4=['gunzip ', ['ss_r', anat_main, '_0d4_mask.nii.gz']];
% cmd5=['rm -rf ', ['ss_r' anat_main, '_0d8.nii.gz']];
% cmd6=['rm -rf ', ['ss_r' anat_main, '_0d4.nii.gz']];
% system(cmd0);
% system(cmd1);
% system(cmd2);
% system(cmd3);
% system(cmd4);
% system(cmd5);
% system(cmd6);
% 
% % generate small brain mask
% matlabbatch = cell(1,3);
% fn_ranat08_mask = fullfile(out_dir, ['ss_r', anat_main, '_0d8_mask.nii']);
% matlabbatch{1,1}.spm.util.imcalc.input = ...
%     {fn_GM_prob,fn_WM_prob,fn_CSF_prob,fn_ranat08_mask};
% matlabbatch{1,1}.spm.util.imcalc.output = strcat('SmallBrainMask.nii');
% matlabbatch{1,1}.spm.util.imcalc.outdir = {out_dir};
% matlabbatch{1,1}.spm.util.imcalc.expression = strcat(['sum(X)>', num2str(thresh_brain)]);
% matlabbatch{1,1}.spm.util.imcalc.options.dmtx = double(1);
% 
% % generate normal brain mask
% fn_ranat04_mask = fullfile(out_dir, ['ss_r', anat_main, '_0d4_mask.nii']);
% matlabbatch{1,2}.spm.util.imcalc.input = ...
%     {fn_GM_prob,fn_WM_prob,fn_CSF_prob,fn_ranat04_mask};
% matlabbatch{1,2}.spm.util.imcalc.output = strcat('BrainMask.nii');
% matlabbatch{1,2}.spm.util.imcalc.outdir = {out_dir};
% matlabbatch{1,2}.spm.util.imcalc.expression = strcat(['sum(X)>', num2str(thresh_brain)]);
% matlabbatch{1,2}.spm.util.imcalc.options.dmtx = double(1);

% generate WM GM mask
matlabbatch = cell(1,1);
matlabbatch{1,1}.spm.util.imcalc.input = {fn_GM_prob, fn_WM_prob};
matlabbatch{1,1}.spm.util.imcalc.output = strcat('GMWMMask.nii');
matlabbatch{1,1}.spm.util.imcalc.outdir = {out_dir};
matlabbatch{1,1}.spm.util.imcalc.expression = strcat('sum(X)>',num2str(thresh_gmwm));
matlabbatch{1,1}.spm.util.imcalc.options.dmtx = double(1);

% run job
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

%% zip generated images and remove unziped images
delete(MPRAGE, brainmask, func);
cd(out_dir)
system('gzip *.nii');

cd(savecurpath);
