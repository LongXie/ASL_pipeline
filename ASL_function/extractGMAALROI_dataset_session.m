function extractGMAALROI_dataset_session...
    (session_dir, AALROI, aCBF_gz, rCBF_gz, gm_gz, out_folder, out_filename, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 8
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

% check input extension
[gm_dir, gm_filename, ext] = fileparts(gm_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(gm_filename);
    if ~strcmp(ext, '.nii')
        error('Gray matter probability map`s extension must be .nii.gz');
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
ROI_names = {'LPostCingulate', 'RPostCingulate', 'MPostCingulate', ...
             'LPrecuneus', 'RPrecuneus', 'MPrecuneus' ...
             'LFusiform', 'RFusiform', 'MFusiform', ...
             'LHippo', 'RHippo', 'MHippo', ...
             'LParaHippo', 'RParaHippo', 'MParaHippo'};
ROI_values = {[4021], [4022], [4021, 4022], ...
              [6301], [6302], [6301, 6302], ...
              [5401], [5402], [5401, 5402], ...
              [4101], [4102], [4101, 4102], ...
              [4111], [4112], [4111, 4112]};
ROI_combine_type = ...
            [1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1];
          
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
    if exist(aCBF_gz, 'file') && exist(rCBF_gz, 'file') && exist(gm_gz, 'file')
        % unzip file for SPM usage
        aCBF = fullfile(out_dir, aCBF_filename);
        system(sprintf('gunzip -c %s > %s', aCBF_gz, aCBF));
        rCBF = fullfile(out_dir, rCBF_filename);
        system(sprintf('gunzip -c %s > %s', rCBF_gz, rCBF));
        gm = fullfile(out_dir, gm_filename);
        system(sprintf('gunzip -c %s > %s', gm_gz, gm));
        pause(0.5)
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s, %s or %s does not exist in %s.\n', aCBF_gz, rCBF_gz, gm_gz, run_dir);
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

    % load gm prob
    gmvol = spm_read_vols(spm_vol(gm));
    gmvol = gmvol(:);
    gmvol_mask = gmvol >= 0.4;
    
    % compute absolute ROI CBF
%     AALROICBF.names{1} = 'aLPostCingulate';
%     AALROICBF.measures(1) = mean(aCBFvol(AALvol == LPCingulate));
%     AALROICBF.names{2} = 'aRPostCingulate';
%     AALROICBF.measures(2) = mean(aCBFvol(AALvol == RPCingulate));
%     AALROICBF.names{3} = 'aMPostCingulate';
%     AALROICBF.measures(3) = mean(AALROICBF.measures(1:2));
%     AALROICBF.names{4} = 'aLPrecuneus';
%     AALROICBF.measures(4) = mean(aCBFvol(AALvol == LPrecuneus));
%     AALROICBF.names{5} = 'aRPrecuneus';
%     AALROICBF.measures(5) = mean(aCBFvol(AALvol == RPrecuneus));
%     AALROICBF.names{6} = 'aMPrecuneus';
%     AALROICBF.measures(6) = mean(AALROICBF.measures(4:5));
%     AALROICBF.names{7} = 'aLFusiform';
%     AALROICBF.measures(7) = mean(aCBFvol(AALvol == LFusiform));
%     AALROICBF.names{8} = 'aRFusiform';
%     AALROICBF.measures(8) = mean(aCBFvol(AALvol == RFusiform));
%     AALROICBF.names{9} = 'aMFusiform';
%     AALROICBF.measures(9) = mean(AALROICBF.measures(7:8));
%     AALROICBF.names{10} = 'aLHippo';
%     AALROICBF.measures(10) = mean(aCBFvol(AALvol == LHippo));
%     AALROICBF.names{11} = 'aRHippo';
%     AALROICBF.measures(11) = mean(aCBFvol(AALvol == RHippo));
%     AALROICBF.names{12} = 'aMHippo';
%     AALROICBF.measures(12) = mean(AALROICBF.measures(10:11));
%     AALROICBF.names{13} = 'aLParaHippo';
%     AALROICBF.measures(13) = mean(aCBFvol(AALvol == LParaHippo));
%     AALROICBF.names{14} = 'aRParaHippo';
%     AALROICBF.measures(14) = mean(aCBFvol(AALvol == RParaHippo));
%     AALROICBF.names{15} = 'aMParaHippo';
%     AALROICBF.measures(15) = mean(AALROICBF.measures(13:14));
%     
%     % compute relative ROI CBF
%     AALROICBF.names{16} = 'rLPostCingulate';
%     AALROICBF.measures(16) = mean(rCBFvol(AALvol == LPCingulate));
%     AALROICBF.names{17} = 'rRPostCingulate';
%     AALROICBF.measures(17) = mean(rCBFvol(AALvol == RPCingulate));
%     AALROICBF.names{18} = 'rMPostCingulate';
%     AALROICBF.measures(18) = mean(AALROICBF.measures(16:17));
%     AALROICBF.names{19} = 'rLPrecuneus';
%     AALROICBF.measures(19) = mean(rCBFvol(AALvol == LPrecuneus));
%     AALROICBF.names{20} = 'rRPrecuneus';
%     AALROICBF.measures(20) = mean(rCBFvol(AALvol == RPrecuneus));
%     AALROICBF.names{21} = 'rMPrecuneus';
%     AALROICBF.measures(21) = mean(AALROICBF.measures(19:20));
%     AALROICBF.names{22} = 'rLFusiform';
%     AALROICBF.measures(22) = mean(rCBFvol(AALvol == LFusiform));
%     AALROICBF.names{23} = 'rRFusiform';
%     AALROICBF.measures(23) = mean(rCBFvol(AALvol == RFusiform));
%     AALROICBF.names{24} = 'rMFusiform';
%     AALROICBF.measures(24) = mean(AALROICBF.measures(22:23));
%     AALROICBF.names{25} = 'rLHippo';
%     AALROICBF.measures(25) = mean(rCBFvol(AALvol == LHippo));
%     AALROICBF.names{26} = 'rRHippo';
%     AALROICBF.measures(26) = mean(rCBFvol(AALvol == RHippo));
%     AALROICBF.names{27} = 'rMHippo';
%     AALROICBF.measures(27) = mean(AALROICBF.measures(25:26));
%     AALROICBF.names{28} = 'rLParaHippo';
%     AALROICBF.measures(28) = mean(rCBFvol(AALvol == LParaHippo));
%     AALROICBF.names{29} = 'rRParaHippo';
%     AALROICBF.measures(29) = mean(rCBFvol(AALvol == RParaHippo));
%     AALROICBF.names{30} = 'rMParaHippo';
%     AALROICBF.measures(30) = mean(AALROICBF.measures(28:29));
    
    % compute absolute
    n = 1;
    AALROICBF = [];
    CBFvol = aCBFvol;
    for ii = 1:length(ROI_names)
        AALROICBF.names{n} = ['a' ROI_names{ii}];
        ROI_value = ROI_values{ii};
        if ROI_combine_type(ii) == 1
            temp = zeros(length(ROI_value), 1);
            for jj = 1:length(ROI_value)
                temp(jj) = mean(CBFvol(AALvol == ROI_value(jj) & gmvol_mask));
            end
            AALROICBF.measures(n) = mean(temp);
        elseif ROI_combine_type(ii) == 2
            temp = zeros(size(AALvol));
            for jj = 1:length(ROI_value)
                temp = temp | (AALvol == ROI_value(jj));
            end
            AALROICBF.measures(n) = mean(CBFvol(temp & gmvol_mask));
        else
            AALROICBF.measures(n) = mean([AALROICBF.measures(n-1), ...
                                          AALROICBF.measures(n-2)]);
        end
        n = n + 1;
    end
    
    % compute relative
    CBFvol = rCBFvol;
    for ii = 1:length(ROI_names)
        AALROICBF.names{n} = ['r' ROI_names{ii}];
        ROI_value = ROI_values{ii};
        if ROI_combine_type(ii) == 1
            temp = zeros(length(ROI_value), 1);
            for jj = 1:length(ROI_value)
                temp(jj) = mean(CBFvol(AALvol == ROI_value(jj) & gmvol_mask));
            end
            AALROICBF.measures(n) = mean(temp);
        elseif ROI_combine_type(ii) == 2
            temp = zeros(size(AALvol));
            for jj = 1:length(ROI_value)
                temp = temp | (AALvol == ROI_value(jj));
            end
            AALROICBF.measures(n) = mean(CBFvol(temp & gmvol_mask));
        else
            AALROICBF.measures(n) = mean([AALROICBF.measures(n-1), ...
                                          AALROICBF.measures(n-2)]);
        end
        n = n + 1;
    end
    

    % save ROI CBF
    out_file = fullfile(out_dir, out_filename);
    save(out_file, 'AALROICBF');
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(aCBF, rCBF,gm);
    
end

% back to original dir
cd(savecurpath);
