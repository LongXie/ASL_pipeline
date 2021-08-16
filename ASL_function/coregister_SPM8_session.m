function coregister_SPM8_session...
    (session_dir, MPRAGE_gz, meanfunc_gz, func_gz, out_folder, out_prefix, LOGTXT)


ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)


%% check input and output
if nargin < 7
    error('Not enough input argument!');
end

% check input MPRAGE extension
[anat_dir, anat_filename, ext] = fileparts(MPRAGE_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(anat_filename);
    if ~strcmp(ext, '.nii')
        error('Anatomical scan`s extension must be .nii.gz');
    end
else
    error('Anatomical scan`s extension must be .nii.gz');
end

% check input mean func extension
[~, meanfunc_filename, ext] = fileparts(meanfunc_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(meanfunc_filename);
    if ~strcmp(ext, '.nii')
        error('Mean functional scan`s extension must be .nii.gz');
    end
else
    error('Mean functional scan`s extension must be .nii.gz');
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

%% parameters
global defaults;
spm_get_defaults;
flags = defaults.coreg;

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

%% check if anatomical scan exist or not
%anat_dir = fullfile(session_dir, 'MPRAGE');
%MPRAGE_gz = fullfile(anat_dir, 'MPRAGE.nii.gz');
if ~exist(MPRAGE_gz, 'file')
    % if not exist, report error
    msg = sprintf('ERROR: MPRAGE does not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
else
    % unzip anatomical scan
    %[~, anat_filename] = fileparts(MPRAGE_gz);
    MPRAGE = fullfile(anat_dir, anat_filename);
    system(sprintf('gunzip -c %s > %s', MPRAGE_gz, MPRAGE));
end 

%% Run coregistration
savecurpath = pwd;

for r = 1:nruns   
    
    %% initialization
    fprintf('Coregistration for run %s (%0.0f/%0.0f).\n',...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % generate output folder
    out_dir = fullfile(run_dir, out_folder);
    if isdir(out_dir)
        rmdir(out_dir,'s');
    end
    mkdir(out_dir);
    
    %% check and unzip input files
    cd(run_dir);
    %meanfunc_gz  = fullfile('mc', 'meanRAWPCASL.nii.gz');
    %func_gz = fullfile('mc', 'rRAWPCASL.nii.gz');
    
    % check if file exist
    if exist(meanfunc_gz, 'file') && exist(func_gz, 'file')
        % unzip file for SPM usage
        %[~, meanfunc_filename] = fileparts(meanfunc_gz);
        meanfunc_filename_prefix = [out_prefix, meanfunc_filename];
        meanfunc = fullfile(out_dir, meanfunc_filename_prefix);
        system(sprintf('gunzip -c %s > %s', meanfunc_gz, meanfunc));
        
        %[~, func_filename] = fileparts(func_gz);
        func_filename_prefix = [out_prefix, func_filename];
        func = fullfile(out_dir, func_filename_prefix);
        system(sprintf('gunzip -c %s > %s', func_gz, func));
        pause(0.5)
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s or %s does not exist in %s.\n', ...
            meanfunc_gz, func_gz, run_dir);
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end    
    
    %% perform coregistration
    cd(out_dir)
    
    % read anatomical image
    PG = spm_select('FPList', anat_dir, anat_filename);
    PG = PG(1,:);
    VG = spm_vol(PG);

    % read mean functional image
    PF = spm_select('FPList', out_dir, meanfunc_filename_prefix);
    PF = PF(1,:);
    VF = spm_vol(PF);
    
    % Coregister Source to Target
    %this method from spm_coreg_ui.m
    %get coregistration parameters
    x  = spm_coreg(VG, VF, flags.estimate);
    %get the transformation to be applied with parameters 'x'
    M  = inv(spm_matrix(x));
    
    %transform the mean image;
    clear PO;
    PO = spm_select('ExtFPList', out_dir, func_filename_prefix, Inf);
    if isempty(PO) | PO=='/'
        PO = PF;
    else
        PO = strvcat(PF,PO);
    end
    
    %in MM we put the transformations for the 'other' images
    MM = zeros(4,4,size(PO,1));
    for j=1:size(PO,1)
        %get the transformation matrix for every image
        MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
    end
    for j=1:size(PO,1)
        %write the transformation applied to every image
        spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
    end
    
    %% remove the original unzip file and zip the outputfile
    system('gzip *.nii');
    
end

% remove unziped MPRAGE
delete(MPRAGE);

cd(savecurpath);
