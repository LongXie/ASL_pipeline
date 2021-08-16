function extractFuncROI_session...
    (session_dir, FuncROI, aCBF_gz, rCBF_gz, out_folder, out_filename, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 7
    error('Not enough input argument!');
end

% check input extension
[aCBF_dir, aCBF_filename, ext] = fileparts(aCBF_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(aCBF_filename);
    if ~strcmp(ext, '.nii')
        error('Absolute CBF`s extension must be .nii.gz');
    end
else
    error('Absolute CBF`s extension must be .nii.gz');
end

% check input extension
[rCBF_dir, rCBF_filename, ext] = fileparts(rCBF_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(rCBF_filename);
    if ~strcmp(ext, '.nii')
        error('Relative CBF`s extension must be .nii.gz');
    end
else
    error('Relative CBF`s extension must be .nii.gz');
end

% check Func ROI file
[~, ~, ext] = fileparts(FuncROI);
if ~strcmp(ext, '.nii')
    error('Func ROI`s extension must be .nii');
else
    if ~exist(FuncROI, 'file')
        error('Func ROI does not exist.');
    end
end

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
    fprintf('Extract Func ROI CBF for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});

    % output dir
    out_dir = fullfile(run_dir, out_folder);
    if exist(out_dir, 'dir')
        rmdir(out_dir, 's');
    end
    mkdir(out_dir);
    
    
    %% check and unzip aCBF and rCBF file
    cd(run_dir);
    
    % check if file exist
    if exist(aCBF_gz, 'file') && exist(rCBF_gz, 'file') 
        % unzip file for SPM usage
        aCBF = fullfile(out_dir, aCBF_filename);
        system(sprintf('gunzip -c %s > %s', aCBF_gz, aCBF));
        rCBF = fullfile(out_dir, rCBF_filename);
        system(sprintf('gunzip -c %s > %s', rCBF_gz, rCBF));
        pause(0.5)
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s or %s does not exist in %s.\n', aCBF_gz, rCBF_gz, run_dir);
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end
    
    %% extract AAL ROI CBF
    % load AAL
    Funcvol = spm_read_vols(spm_vol(FuncROI));
    Funcvol = Funcvol(:);
    
    % load aCBF
    aCBFvol = spm_read_vols(spm_vol(aCBF));
    aCBFvol = aCBFvol(:);
    
    % load rCBF
    rCBFvol = spm_read_vols(spm_vol(rCBF));
    rCBFvol = rCBFvol(:);

    % compute absolute ROI CBF
    FuncROICBF.aCBF = mean(aCBFvol(Funcvol == 1));
    
    % compute relative ROI CBF
    FuncROICBF.rCBF = mean(rCBFvol(Funcvol == 1));
    
    % save ROI CBF
    out_file = fullfile(out_dir, out_filename);
    save(out_file, 'FuncROICBF');
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(aCBF, rCBF);
    
end

% back to original dir
cd(savecurpath);
