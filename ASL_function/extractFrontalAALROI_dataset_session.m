function extractFrontalAALROI_dataset_session...
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



ROI_names = {'LFrontal_Sup', 'RFrontal_Sup', 'MFrontal_Sup', ...
             'LFrontal_Sup_Orb', 'RFrontal_Sup_Orb', 'MFrontal_Sup_Orb' ...
             'LFrontal_Mid', 'RFrontal_Mid', 'MFrontal_Mid', ...
             'LFrontal_Mid_Orb', 'RFrontal_Mid_Orb', 'MFrontal_Mid_Orb', ...
             'LFrontal_Inf_Oper', 'RFrontal_Inf_Oper', 'MFrontal_Inf_Oper', ...
             'LFrontal_Inf_Tri', 'RFrontal_Inf_Tri', 'MFrontal_Inf_Tri', ...
             'LFrontal_Inf_Orb', 'RFrontal_Inf_Orb', 'MFrontal_Inf_Orb', ...
             'LRolandic_Oper', 'RRolandic_Oper', 'MRolandic_Oper', ...
             'LFrontal_Sup_Medial', 'RFrontal_Sup_Medial', 'MFrontal_Sup_Medial', ...
             'LFrontal_Med_Orb', 'RFrontal_Med_Orb', 'MFrontal_Med_Orb', ...
             'LRectus', 'RRectus', 'MRectus', ...
             'LCingulum_Ant', 'RCingulum_Ant', 'MCingulum_Ant', ...
             'LFrontal_Orbital', 'RFrontal_Orbital', ...
             'MFrontal_Orbital'};
ROI_values = {[2101], [2102], [2101, 2102], ...
              [2111], [2112], [2111, 2112], ...
              [2201], [2202], [2201, 2202], ...
              [2211], [2212], [2211, 2212], ...
              [2301], [2302], [2301, 2302], ...
              [2311], [2312], [2311, 2312], ...
              [2321], [2322], [2321, 2322], ...
              [2331], [2332], [2331, 2332], ...
              [2601], [2602], [2601, 2602], ...
              [2611], [2612], [2611, 2612], ...
              [2701], [2702], [2701, 2702], ...
              [4001], [4002], [4001, 4002], ...
              [2111, 2211, 2321, 2611], [2112, 2212, 2322, 2612], ...
              [-1]};
ROI_combine_type = ...
            [1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             1 1 1 ...
             2 2 ...
             3 ...
             ];
          
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
                temp(jj) = mean(CBFvol(AALvol == ROI_value(jj)));
            end
            AALROICBF.measures(n) = mean(temp);
        elseif ROI_combine_type(ii) == 2
            temp = zeros(size(AALvol));
            for jj = 1:length(ROI_value)
                temp = temp | (AALvol == ROI_value(jj));
            end
            AALROICBF.measures(n) = mean(CBFvol(temp));
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
                temp(jj) = mean(CBFvol(AALvol == ROI_value(jj)));
            end
            AALROICBF.measures(n) = mean(temp);
        elseif ROI_combine_type(ii) == 2
            temp = zeros(size(AALvol));
            for jj = 1:length(ROI_value)
                temp = temp | (AALvol == ROI_value(jj));
            end
            AALROICBF.measures(n) = mean(CBFvol(temp));
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
    delete(aCBF, rCBF);
    
end

% back to original dir
cd(savecurpath);
