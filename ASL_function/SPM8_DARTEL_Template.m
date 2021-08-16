function SPM8_DARTEL_Template(study_dir, subj_txt, timepoint_txt, anat_folder, LOGTXT)
    
ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

if nargin < 5
    error('Not enough input argument!');
end

%% check whether subj_txt and timepoint_txt exist
if exist(subj_txt, 'file')
    IDs = importsubj(subj_txt);
    nsubj = length(IDs);
    if nsubj <= 0
        error('%s is empty.', subj_txt);
    end
else
    error('%s does not exist.', subj_txt);
end

if exist(timepoint_txt, 'file')
    TimePoint = importsubj(timepoint_txt);
    ntp = length(TimePoint);
    if ntp <= 0
        error('%s is empty.', timepoint_txt);
    end
else
    error('%s does not exist.', timepoint_txt);
end

%% check if all the data exist
for ii = 1:nsubj
    for jj = 1:ntp
       
        % initialization
        id = IDs{ii};
        tp = TimePoint{jj};

        session_dir = fullfile(study_dir, id, tp);  
        if ~exist(session_dir, 'dir')
            continue;
        end
        
        anat_dir = fullfile(session_dir, 'MPRAGE', anat_folder);
        rc1_file_gz = fullfile(anat_dir, 'rc1MPRAGE.nii.gz');
        rc2_file_gz = fullfile(anat_dir, 'rc2MPRAGE.nii.gz');
        
        % check if rc1 exist and unzip
        if ~exist(rc1_file_gz, 'file')
            % if not exist, report error
            msg = sprintf('ERROR: rc1 of subject %s timepoint %s does not exist.\n', id, tp);
            system(sprintf('echo "%s" >> %s', msg, LOGTXT));
            error(msg);
        end
        
        % check if rc2 exist and unzip
        if ~exist(rc2_file_gz, 'file')
            % if not exist, report error
            msg = sprintf('ERROR: rc2 of subject %s timepoint %s does not exist.\n', id, tp);
            system(sprintf('echo "%s" >> %s', msg, LOGTXT));
            error(msg);
        end
        
    end
end

%% unzip and prepare for DARTEL
%rc1_files = cell(nsubj * ntp, 1);
%rc2_files = cell(nsubj * ntp, 1);
idx = 1;

for ii = 1:nsubj
    for jj = 1:ntp
        
        % initialization
        id = IDs{ii};
        tp = TimePoint{jj};
        
        session_dir = fullfile(study_dir, id, tp);
        if ~exist(session_dir, 'dir')
            continue;
        end
        
        anat_dir = fullfile(session_dir, 'MPRAGE', anat_folder);
        rc1_file_gz = fullfile(anat_dir, 'rc1MPRAGE.nii.gz');
        rc2_file_gz = fullfile(anat_dir, 'rc2MPRAGE.nii.gz');
        
        % unzip file for SPM usage
        [~, rc1_filename] = fileparts(rc1_file_gz);
        rc1_file = fullfile(anat_dir, rc1_filename);
        system(sprintf('gunzip -c %s > %s', rc1_file_gz, rc1_file));
        rc1_files{idx, 1} = fullfile(anat_dir, rc1_filename);
        
        % unzip file for SPM usage
        [~, rc2_filename] = fileparts(rc2_file_gz);
        rc2_file = fullfile(anat_dir, rc2_filename);
        system(sprintf('gunzip -c %s > %s', rc2_file_gz, rc2_file));
        rc2_files{idx, 1} = fullfile(anat_dir, rc2_filename);
        
        % prepapre for next iteration
        idx = idx + 1;
        
    end
end
pause(0.5)

%% run DARTEL
clear matlabbatch
%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

matlabbatch{1}.spm.tools.dartel.warp.images = {
    rc1_files
    rc2_files
    }';
matlabbatch{1}.spm.tools.dartel.warp.settings.template = ['Template'];
matlabbatch{1}.spm.tools.dartel.warp.settings.rform = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).rparam = [4 2 1e-006];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(1).slam = 16;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).rparam = [2 1 1e-006];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).K = 0;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(2).slam = 8;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).rparam = [1 0.5 1e-006];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).K = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(3).slam = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).rparam = [0.5 0.25 1e-006];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).K = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(4).slam = 2;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).rparam = [0.25 0.125 1e-006];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).K = 4;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(5).slam = 1;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).its = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).rparam = [0.125 0.0625 1e-006];
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).K = 6;
matlabbatch{1}.spm.tools.dartel.warp.settings.param(6).slam = 0.5;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.lmreg = 0.01;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.cyc = 3;
matlabbatch{1}.spm.tools.dartel.warp.settings.optim.its = 3;

spm_jobman('initcfg');
spm_jobman('run',matlabbatch);

%% delete all the unziped files
for ii = 1:length(rc1_files)
    delete(rc1_files{ii,1}, rc2_files{ii,1});
end

end
