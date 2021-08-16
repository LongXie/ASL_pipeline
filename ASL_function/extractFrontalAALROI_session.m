function extractFrontalAALROI_session...
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

Frontal_Sup_L = 2101;
Frontal_Sup_R = 2102;
Frontal_Sup_Orb_L = 2111;
Frontal_Sup_Orb_R = 2112;
Frontal_Mid_L = 2201;
Frontal_Mid_R = 2202;
Frontal_Mid_Orb_L = 2211;
Frontal_Mid_Orb_R = 2212;
Frontal_Inf_Oper_L = 2301;
Frontal_Inf_Oper_R = 2302;
Frontal_Inf_Tri_L = 2311;
Frontal_Inf_Tri_R = 2312;
Frontal_Inf_Orb_L = 2321;
Frontal_Inf_Orb_R = 2322;
Rolandic_Oper_L = 2331;
Rolandic_Oper_R = 2332;
Frontal_Sup_Medial_L = 2601;
Frontal_Sup_Medial_R = 2602;
Frontal_Med_Orb_L = 2611;
Frontal_Med_Orb_R = 2612;
Rectus_L = 2701;
Rectus_R = 2702;
Cingulum_Ant_L = 4001;
Cingulum_Ant_R = 4002;

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
    AALROICBF.aCBF.Frontal_Sup_L = mean(aCBFvol(AALvol == Frontal_Sup_L));
    AALROICBF.aCBF.Frontal_Sup_R = mean(aCBFvol(AALvol == Frontal_Sup_R));
    AALROICBF.aCBF.Frontal_Sup_M = mean([AALROICBF.aCBF.Frontal_Sup_L, AALROICBF.aCBF.Frontal_Sup_R]);
    
    AALROICBF.aCBF.Frontal_Sup_Orb_L = mean(aCBFvol(AALvol == Frontal_Sup_Orb_L));
    AALROICBF.aCBF.Frontal_Sup_Orb_R = mean(aCBFvol(AALvol == Frontal_Sup_Orb_R));
    AALROICBF.aCBF.Frontal_Sup_Orb_M = mean([AALROICBF.aCBF.Frontal_Sup_Orb_L, AALROICBF.aCBF.Frontal_Sup_Orb_R]);
    
    AALROICBF.aCBF.Frontal_Mid_L = mean(aCBFvol(AALvol == Frontal_Mid_L));
    AALROICBF.aCBF.Frontal_Mid_R = mean(aCBFvol(AALvol == Frontal_Mid_R));
    AALROICBF.aCBF.Frontal_Mid_M = mean([AALROICBF.aCBF.Frontal_Mid_L, AALROICBF.aCBF.Frontal_Mid_R]);
    
    AALROICBF.aCBF.Frontal_Mid_Orb_L = mean(aCBFvol(AALvol == Frontal_Mid_Orb_L));
    AALROICBF.aCBF.Frontal_Mid_Orb_R = mean(aCBFvol(AALvol == Frontal_Mid_Orb_R));
    AALROICBF.aCBF.Frontal_Mid_Orb_M = mean([AALROICBF.aCBF.Frontal_Mid_Orb_L, AALROICBF.aCBF.Frontal_Mid_Orb_R]);
    
    AALROICBF.aCBF.Frontal_Inf_Oper_L = mean(aCBFvol(AALvol == Frontal_Inf_Oper_L));
    AALROICBF.aCBF.Frontal_Inf_Oper_R = mean(aCBFvol(AALvol == Frontal_Inf_Oper_R));
    AALROICBF.aCBF.Frontal_Inf_Oper_M = mean([AALROICBF.aCBF.Frontal_Inf_Oper_L, AALROICBF.aCBF.Frontal_Inf_Oper_R]);

    AALROICBF.aCBF.Frontal_Inf_Tri_L = mean(aCBFvol(AALvol == Frontal_Inf_Tri_L));
    AALROICBF.aCBF.Frontal_Inf_Tri_R = mean(aCBFvol(AALvol == Frontal_Inf_Tri_R));
    AALROICBF.aCBF.Frontal_Inf_Tri_M = mean([AALROICBF.aCBF.Frontal_Inf_Tri_L, AALROICBF.aCBF.Frontal_Inf_Tri_R]);
    
    AALROICBF.aCBF.Frontal_Inf_Orb_L = mean(aCBFvol(AALvol == Frontal_Inf_Orb_L));
    AALROICBF.aCBF.Frontal_Inf_Orb_R = mean(aCBFvol(AALvol == Frontal_Inf_Orb_R));
    AALROICBF.aCBF.Frontal_Inf_Orb_M = mean([AALROICBF.aCBF.Frontal_Inf_Orb_L, AALROICBF.aCBF.Frontal_Inf_Orb_R]);
    
    AALROICBF.aCBF.Rolandic_Oper_L = mean(aCBFvol(AALvol == Rolandic_Oper_L));
    AALROICBF.aCBF.Rolandic_Oper_R = mean(aCBFvol(AALvol == Rolandic_Oper_R));
    AALROICBF.aCBF.Rolandic_Oper_M = mean([AALROICBF.aCBF.Rolandic_Oper_L, AALROICBF.aCBF.Rolandic_Oper_R]);
    
    AALROICBF.aCBF.Frontal_Sup_Medial_L = mean(aCBFvol(AALvol == Frontal_Sup_Medial_L));
    AALROICBF.aCBF.Frontal_Sup_Medial_R = mean(aCBFvol(AALvol == Frontal_Sup_Medial_R));
    AALROICBF.aCBF.Frontal_Sup_Medial_M = mean([AALROICBF.aCBF.Frontal_Sup_Medial_L, AALROICBF.aCBF.Frontal_Sup_Medial_R]);
    
    AALROICBF.aCBF.Frontal_Med_Orb_L = mean(aCBFvol(AALvol == Frontal_Med_Orb_L));
    AALROICBF.aCBF.Frontal_Med_Orb_R = mean(aCBFvol(AALvol == Frontal_Med_Orb_R));
    AALROICBF.aCBF.Frontal_Med_Orb_M = mean([AALROICBF.aCBF.Frontal_Med_Orb_L, AALROICBF.aCBF.Frontal_Med_Orb_R]);

    AALROICBF.aCBF.Rectus_L = mean(aCBFvol(AALvol == Rectus_L));
    AALROICBF.aCBF.Rectus_R = mean(aCBFvol(AALvol == Rectus_R));
    AALROICBF.aCBF.Rectus_M = mean([AALROICBF.aCBF.Rectus_L, AALROICBF.aCBF.Rectus_R]);
    
    AALROICBF.aCBF.Cingulum_Ant_L = mean(aCBFvol(AALvol == Cingulum_Ant_L));
    AALROICBF.aCBF.Cingulum_Ant_R = mean(aCBFvol(AALvol == Cingulum_Ant_R));
    AALROICBF.aCBF.Cingulum_Ant_M = mean([AALROICBF.aCBF.Cingulum_Ant_L, AALROICBF.aCBF.Cingulum_Ant_R]);

    
    
    % compute relative ROI CBF
    AALROICBF.rCBF.Frontal_Sup_L = mean(rCBFvol(AALvol == Frontal_Sup_L));
    AALROICBF.rCBF.Frontal_Sup_R = mean(rCBFvol(AALvol == Frontal_Sup_R));
    AALROICBF.rCBF.Frontal_Sup_M = mean([AALROICBF.rCBF.Frontal_Sup_L, AALROICBF.rCBF.Frontal_Sup_R]);
    
    AALROICBF.rCBF.Frontal_Sup_Orb_L = mean(rCBFvol(AALvol == Frontal_Sup_Orb_L));
    AALROICBF.rCBF.Frontal_Sup_Orb_R = mean(rCBFvol(AALvol == Frontal_Sup_Orb_R));
    AALROICBF.rCBF.Frontal_Sup_Orb_M = mean([AALROICBF.rCBF.Frontal_Sup_Orb_L, AALROICBF.rCBF.Frontal_Sup_Orb_R]);
    
    AALROICBF.rCBF.Frontal_Mid_L = mean(rCBFvol(AALvol == Frontal_Mid_L));
    AALROICBF.rCBF.Frontal_Mid_R = mean(rCBFvol(AALvol == Frontal_Mid_R));
    AALROICBF.rCBF.Frontal_Mid_M = mean([AALROICBF.rCBF.Frontal_Mid_L, AALROICBF.rCBF.Frontal_Mid_R]);
    
    AALROICBF.rCBF.Frontal_Mid_Orb_L = mean(rCBFvol(AALvol == Frontal_Mid_Orb_L));
    AALROICBF.rCBF.Frontal_Mid_Orb_R = mean(rCBFvol(AALvol == Frontal_Mid_Orb_R));
    AALROICBF.rCBF.Frontal_Mid_Orb_M = mean([AALROICBF.rCBF.Frontal_Mid_Orb_L, AALROICBF.rCBF.Frontal_Mid_Orb_R]);
    
    AALROICBF.rCBF.Frontal_Inf_Oper_L = mean(rCBFvol(AALvol == Frontal_Inf_Oper_L));
    AALROICBF.rCBF.Frontal_Inf_Oper_R = mean(rCBFvol(AALvol == Frontal_Inf_Oper_R));
    AALROICBF.rCBF.Frontal_Inf_Oper_M = mean([AALROICBF.rCBF.Frontal_Inf_Oper_L, AALROICBF.rCBF.Frontal_Inf_Oper_R]);

    AALROICBF.rCBF.Frontal_Inf_Tri_L = mean(rCBFvol(AALvol == Frontal_Inf_Tri_L));
    AALROICBF.rCBF.Frontal_Inf_Tri_R = mean(rCBFvol(AALvol == Frontal_Inf_Tri_R));
    AALROICBF.rCBF.Frontal_Inf_Tri_M = mean([AALROICBF.rCBF.Frontal_Inf_Tri_L, AALROICBF.rCBF.Frontal_Inf_Tri_R]);
    
    AALROICBF.rCBF.Frontal_Inf_Orb_L = mean(rCBFvol(AALvol == Frontal_Inf_Orb_L));
    AALROICBF.rCBF.Frontal_Inf_Orb_R = mean(rCBFvol(AALvol == Frontal_Inf_Orb_R));
    AALROICBF.rCBF.Frontal_Inf_Orb_M = mean([AALROICBF.rCBF.Frontal_Inf_Orb_L, AALROICBF.rCBF.Frontal_Inf_Orb_R]);
    
    AALROICBF.rCBF.Rolandic_Oper_L = mean(rCBFvol(AALvol == Rolandic_Oper_L));
    AALROICBF.rCBF.Rolandic_Oper_R = mean(rCBFvol(AALvol == Rolandic_Oper_R));
    AALROICBF.rCBF.Rolandic_Oper_M = mean([AALROICBF.rCBF.Rolandic_Oper_L, AALROICBF.rCBF.Rolandic_Oper_R]);
    
    AALROICBF.rCBF.Frontal_Sup_Medial_L = mean(rCBFvol(AALvol == Frontal_Sup_Medial_L));
    AALROICBF.rCBF.Frontal_Sup_Medial_R = mean(rCBFvol(AALvol == Frontal_Sup_Medial_R));
    AALROICBF.rCBF.Frontal_Sup_Medial_M = mean([AALROICBF.rCBF.Frontal_Sup_Medial_L, AALROICBF.rCBF.Frontal_Sup_Medial_R]);
    
    AALROICBF.rCBF.Frontal_Med_Orb_L = mean(rCBFvol(AALvol == Frontal_Med_Orb_L));
    AALROICBF.rCBF.Frontal_Med_Orb_R = mean(rCBFvol(AALvol == Frontal_Med_Orb_R));
    AALROICBF.rCBF.Frontal_Med_Orb_M = mean([AALROICBF.rCBF.Frontal_Med_Orb_L, AALROICBF.rCBF.Frontal_Med_Orb_R]);

    AALROICBF.rCBF.Rectus_L = mean(rCBFvol(AALvol == Rectus_L));
    AALROICBF.rCBF.Rectus_R = mean(rCBFvol(AALvol == Rectus_R));
    AALROICBF.rCBF.Rectus_M = mean([AALROICBF.rCBF.Rectus_L, AALROICBF.rCBF.Rectus_R]);
    
    AALROICBF.rCBF.Cingulum_Ant_L = mean(rCBFvol(AALvol == Cingulum_Ant_L));
    AALROICBF.rCBF.Cingulum_Ant_R = mean(rCBFvol(AALvol == Cingulum_Ant_R));
    AALROICBF.rCBF.Cingulum_Ant_M = mean([AALROICBF.rCBF.Cingulum_Ant_L, AALROICBF.rCBF.Cingulum_Ant_R]);

    % save ROI CBF
    out_file = fullfile(out_dir, out_filename);
    save(out_file, 'AALROICBF');
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(aCBF, rCBF);
    
end

% back to original dir
cd(savecurpath);
