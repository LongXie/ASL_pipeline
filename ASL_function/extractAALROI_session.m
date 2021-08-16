function extractAALROI_session...
    (session_dir, AALROI, aCBF_gz, rCBF_gz, out_folder, out_filename, LOGTXT)

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

% check AAL ROI file
[~, ~, ext] = fileparts(AALROI);
if ~strcmp(ext, '.nii')
    error('AAL ROI`s extension must be .nii');
else
    if ~exist(AALROI, 'file')
        error('AAL ROI does not exist.');
    end
end

%% parameters
LPCingulate = 4021;
RPCingulate = 4022;
LPrecuneus = 6301;
RPrecuneus = 6302;
LFusiform = 5401;
RFusiform = 5402;
LHippo = 4101;
RHippo = 4102;
LParaHippo = 4111;
RParaHippo = 4112;

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
    fprintf('Extract AAL ROI CBF for run %s (%0.0f/%0.0f).\n', ...
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
    AALvol = spm_read_vols(spm_vol(AALROI));
    AALvol = AALvol(:);
    
    % load aCBF
    aCBFvol = spm_read_vols(spm_vol(aCBF));
    aCBFvol = aCBFvol(:);
    
    % load rCBF
    rCBFvol = spm_read_vols(spm_vol(rCBF));
    rCBFvol = rCBFvol(:);

    % compute absolute ROI CBF
    AALROICBF.aCBF.LPCingulate = mean(aCBFvol(AALvol == LPCingulate));
    AALROICBF.aCBF.RPCingulate = mean(aCBFvol(AALvol == RPCingulate));
    AALROICBF.aCBF.MPCingulate = mean([AALROICBF.aCBF.LPCingulate, AALROICBF.aCBF.RPCingulate]);
    AALROICBF.aCBF.LPrecuneus = mean(aCBFvol(AALvol == LPrecuneus));
    AALROICBF.aCBF.RPrecuneus = mean(aCBFvol(AALvol == RPrecuneus));
    AALROICBF.aCBF.MPrecuneus = mean([AALROICBF.aCBF.LPrecuneus, AALROICBF.aCBF.RPrecuneus]);
    AALROICBF.aCBF.LFusiform = mean(aCBFvol(AALvol == LFusiform));
    AALROICBF.aCBF.RFusiform = mean(aCBFvol(AALvol == RFusiform));
    AALROICBF.aCBF.MFusiform = mean([AALROICBF.aCBF.LFusiform, AALROICBF.aCBF.RFusiform]);
    AALROICBF.aCBF.LHippo = mean(aCBFvol(AALvol == LHippo));
    AALROICBF.aCBF.RHippo = mean(aCBFvol(AALvol == RHippo));
    AALROICBF.aCBF.MHippo = mean([AALROICBF.aCBF.LHippo, AALROICBF.aCBF.RHippo]);
    AALROICBF.aCBF.LParaHippo = mean(aCBFvol(AALvol == LParaHippo));
    AALROICBF.aCBF.RParaHippo = mean(aCBFvol(AALvol == RParaHippo));
    AALROICBF.aCBF.MParaHippo = mean([AALROICBF.aCBF.LParaHippo, AALROICBF.aCBF.RParaHippo]);

    % compute relative ROI CBF
    AALROICBF.rCBF.LPCingulate = mean(rCBFvol(AALvol == LPCingulate));
    AALROICBF.rCBF.RPCingulate = mean(rCBFvol(AALvol == RPCingulate));
    AALROICBF.rCBF.MPCingulate = mean([AALROICBF.rCBF.LPCingulate, AALROICBF.rCBF.RPCingulate]);
    AALROICBF.rCBF.LPrecuneus = mean(rCBFvol(AALvol == LPrecuneus));
    AALROICBF.rCBF.RPrecuneus = mean(rCBFvol(AALvol == RPrecuneus));
    AALROICBF.rCBF.MPrecuneus = mean([AALROICBF.rCBF.LPrecuneus, AALROICBF.rCBF.RPrecuneus]);
    AALROICBF.rCBF.LFusiform = mean(rCBFvol(AALvol == LFusiform));
    AALROICBF.rCBF.RFusiform = mean(rCBFvol(AALvol == RFusiform));
    AALROICBF.rCBF.MFusiform = mean([AALROICBF.rCBF.LFusiform, AALROICBF.rCBF.RFusiform]);
    AALROICBF.rCBF.LHippo = mean(rCBFvol(AALvol == LHippo));
    AALROICBF.rCBF.RHippo = mean(rCBFvol(AALvol == RHippo));
    AALROICBF.rCBF.MHippo = mean([AALROICBF.rCBF.LHippo, AALROICBF.rCBF.RHippo]);
    AALROICBF.rCBF.LParaHippo = mean(rCBFvol(AALvol == LParaHippo));
    AALROICBF.rCBF.RParaHippo = mean(rCBFvol(AALvol == RParaHippo));
    AALROICBF.rCBF.MParaHippo = mean([AALROICBF.rCBF.LParaHippo, AALROICBF.rCBF.RParaHippo]);
    
    % save ROI CBF
    out_file = fullfile(out_dir, out_filename);
    save(out_file, 'AALROICBF');
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(aCBF, rCBF);
    
end

% back to original dir
cd(savecurpath);
