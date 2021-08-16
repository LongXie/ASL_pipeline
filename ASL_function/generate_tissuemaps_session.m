function generate_tissuemaps_session...
    (session_dir, anat_dir, func_gz, out_dir, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 4
    error('Not enough input argument!');
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

thresh_gmwm = 0.1;
SmoothFWHM = 2;
savecurpath = pwd;
defaults=spm_get_defaults;
flags = defaults.coreg;
resFlags = struct(...
    'interp', 1,...                       % 1: trilinear interpolation   0: nearest neighbour
    'wrap', flags.write.wrap,...           % wrapping info (ignore...)
    'mask', flags.write.mask,...           % masking (see spm_reslice)
    'which',1,...                         % write reslice time series for later use, don't write the first 1
    'mean',0);

%% check out directory
if exist(out_dir, 'dir')
   rmdir(out_dir, 's');
end
mkdir(out_dir);


%% check and unzip tissue probability maps
%anat_dir = fullfile(session_dir, 'MPRAGE', 'segmentation');

% tissue probability mapes
gm_prob_gz = fullfile(anat_dir, 'c1MPRAGE.nii.gz');
wm_prob_gz = fullfile(anat_dir, 'c2MPRAGE.nii.gz');
csf_prob_gz = fullfile(anat_dir, 'c3MPRAGE.nii.gz');

% check if exist
if exist(gm_prob_gz, 'file') && exist(wm_prob_gz, 'file') ...
    && exist(csf_prob_gz, 'file')
    % unzip file for SPM usage
    % gm
%     [~, gm_filename] = fileparts(gm_prob_gz);
%     gm_prob = fullfile(out_dir, gm_filename);
%     system(sprintf('gunzip -c %s > %s', gm_prob_gz, gm_prob));
%     % wm
%     [~, wm_filename] = fileparts(wm_prob_gz);
%     wm_prob = fullfile(out_dir, wm_filename);
%     system(sprintf('gunzip -c %s > %s', wm_prob_gz, wm_prob));
%     % csf
%     [~, csf_filename] = fileparts(csf_prob_gz);
%     csf_prob = fullfile(out_dir, csf_filename);
%     system(sprintf('gunzip -c %s > %s', csf_prob_gz, csf_prob));
    sprintf('check ok. \n');
else
    % if not exist, report error
    msg = sprintf('ERROR: Tissue probability maps or brain mask do not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end

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
%     [func_pathname, func_filename] = fileparts(func_gz);
%     func = fullfile(func_pathname, func_filename);
%     system(sprintf('gunzip -c %s > %s', func_gz, func));
%     pause(0.5)
    sprintf('check ok. \n');
else
    % if not exist, report error
    msg = sprintf('ERROR: %s or %s does not exist in %s.\n', func_gz, run_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end


%% use c3d to generate masks
OUTT1MASK = fullfile(anat_dir, 'BrainT1Mask.nii.gz');
OUTMASK = fullfile(out_dir, 'SmallBrainMask.nii.gz');
cmd = sprintf('c3d %s %s -add %s -add -smooth 4mm -thresh 0.5 inf 1 0 -o %s %s %s -smooth 8mm -reslice-identity -thresh 0.5 inf 1 0 -o %s', ...
    gm_prob_gz, wm_prob_gz, csf_prob_gz, OUTT1MASK, func_gz, OUTT1MASK, OUTMASK);
fprintf(cmd)
fprintf('\n')
system(cmd);

% generate a dialated mask
OUTDILATEDMASK = fullfile(out_dir, 'BrainMask.nii.gz');
cmd = sprintf('c3d %s -dilate 1 1x1x1vox -o %s', OUTMASK, OUTDILATEDMASK);
fprintf(cmd)
fprintf('\n')
system(cmd);

% generate resliced GM probability maps
GM_PROB_LOWRES = fullfile(out_dir, 'c1MPRAGE.nii.gz');
cmd = sprintf('c3d %s %s -smooth 1mm -reslice-identity -o %s', ...
    func_gz, gm_prob_gz, GM_PROB_LOWRES);
fprintf(cmd)
fprintf('\n')
system(cmd);

% generate resliced WM probability maps
WM_PROB_LOWRES = fullfile(out_dir, 'c2MPRAGE.nii.gz');
cmd = sprintf('c3d %s %s -smooth 1mm -reslice-identity -o %s', ...
    func_gz, wm_prob_gz, WM_PROB_LOWRES);
fprintf(cmd)
fprintf('\n')
system(cmd);

% generate resliced CSF probability maps
CSF_PROB_LOWRES = fullfile(out_dir, 'c3MPRAGE.nii.gz');
cmd = sprintf('c3d %s %s -smooth 1mm -reslice-identity -o %s', ...
    func_gz, csf_prob_gz, CSF_PROB_LOWRES);
fprintf(cmd)
fprintf('\n')
system(cmd);

% generate GMWM mask
OUTGMWMMASK = fullfile(out_dir, 'GMWMMask.nii.gz');
cmd = sprintf('c3d %s %s -add -thresh %s inf 1 0 -o %s', ...
    GM_PROB_LOWRES, WM_PROB_LOWRES, thresh_gmwm, OUTGMWMMASK);
fprintf(cmd)
fprintf('\n')
system(cmd);

% generate a very small mask only using GM and WM (sudipto)
OUTTINYMASK = fullfile(out_dir, 'TinyBrainMask.nii.gz');
cmd = sprintf('c3d %s -smooth 2mm %s -smooth 2mm -add -thresh 0.2 inf 1 0 -holefill 1 0 -popas SEG %s -push SEG -reslice-identity -thresh 0.1 info 1 0 -o %s', ...
    gm_prob_gz, wm_prob_gz, func_gz, OUTTINYMASK);
fprintf(cmd)
fprintf('\n')
system(cmd);

% %% use SPM to generate masks
% % gray matter prob
% Qim = fullfile(out_dir, ['s', spm_str_manip(gm_prob,'dt')]);
% spm_smooth(gm_prob, Qim, SmoothFWHM);
% P=char(func,Qim);
% spm_reslice(P,resFlags);
% 
% % white matter prob
% Qim = fullfile(out_dir, ['s', spm_str_manip(wm_prob,'dt')]);
% spm_smooth(wm_prob, Qim, SmoothFWHM);
% P=char(func,Qim);
% spm_reslice(P,resFlags);
% 
% % csf matter prob
% Qim = fullfile(out_dir, ['s', spm_str_manip(csf_prob,'dt')]);
% spm_smooth(csf_prob, Qim, SmoothFWHM);
% P=char(func,Qim);
% spm_reslice(P,resFlags);
% 
% %% zip generated images and remove unziped images
% delete(gm_prob, wm_prob, csf_prob, func);
% cd(out_dir)
% system('gzip *.nii');
% 
% cd(savecurpath);

