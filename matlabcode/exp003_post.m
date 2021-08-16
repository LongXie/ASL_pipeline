%% 0. Initialization
clc
clear all
close all

ROOT = '/data/picsl/longxie/WolkMCI/PMC_MCI/ASLPET';
CODEDIR = fullfile(ROOT, 'code', 'matlabcode');
ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
DATADIR = fullfile(ROOT, 'recon');

addpath(CODEDIR)
addpath(ASLFUNCDIR)

% scan information
SUBJ_TXT = fullfile(ROOT, 'info', 'subj.txt');
%DEMOGRAPHIC_TXT = fullfile(ROOT, 'analysis_input', 'METF_demog.csv');
fid = fopen(SUBJ_TXT);
IDs = textscan(fid, '%s');
IDs = IDs{1};
fclose(fid);
nsubj = length(IDs);

TIMEPOINT_TXT = fullfile(ROOT, 'analysis_input', 'timepoint.txt');
fid = fopen(TIMEPOINT_TXT);
TimePoint = textscan(fid, '%s');
TimePoint = TimePoint{1};
fclose(fid);
ntp = length(TimePoint);

% experiment ID
expid = '003';
EXPDIR = fullfile(ROOT, 'exp', ['exp', expid, '_PVT']);
if ~exist(EXPDIR, 'dir')
    mkdir(EXPDIR);
end
cd(EXPDIR)

% run information
RUN_TXT = fullfile(ROOT, 'analysis_input', 'run.txt');
fid = fopen(RUN_TXT);
RUNs = textscan(fid, '%s');
RUNs = RUNs{1};
fclose(fid);
nruns = length(RUNs);

%% 1. summarize global CBF information
% save the result in result directory

% global cbf
fprintf('\n\nSummarizing global CBF.\n');

idx = 1;
IDList = [];
TPList = [];
IDTPList = [];
globalCBF_temp = [];

for ii = 1:nsubj
    for jj = 1:ntp
        
        % init
        id = IDs{ii};
        tp = TimePoint{jj};
        
        if exist(fullfile(DATADIR, id, tp), 'dir')
        
            % ID information
            IDList{idx,1} = id;
            TPList{idx,1} = num2str(jj);
            IDTPList{idx,1} = [id, '_', num2str(jj)];
            
            for kk = 1:nruns
                
                % initialization
                run = RUNs{kk};
                run_dir = fullfile(DATADIR, id, tp, run);
                
                if exist(run_dir, 'dir')
                
                    clear globalCBF
                    gbcbf = fullfile(run_dir, 'clean_SCORE', 'globalCBF.mat');
                    load(gbcbf);
                    nmeasures = length(globalCBF.measures);
                    start = nmeasures*(kk-1)+1;
                    globalCBF_temp(idx, start:(start+nmeasures-1)) = globalCBF.measures;
                    
                    
                else
                    
                    start = nmeasures*(kk-1)+1;
                    globalCBF_temp(idx, start:(start+nmeasures-1)) = nan;
                    fprintf('%s %s does not have %s.\n',...
                        id, tp, run);
                end
                
            end
            
            % next iteration
            idx = idx + 1;
            
        end
    end    
end

% make array names
ARRAYNAME = [];
for kk = 1:nruns
    for ll = 1:nmeasures
        ARRAYNAME = [ARRAYNAME, ',''', globalCBF.names{ll}, '_', RUNs{kk}, ''''];               
    end
end

% construct dataset
eval(['globalCBF_all = dataset({IDList, ''ID''},', ...
      '{TPList,''Time''},{IDTPList,''ID_TIME''},', ...
      '{globalCBF_temp', ARRAYNAME, '});']);
    

% save summarized global CBF
save(fullfile(EXPDIR, 'globalCBF_all.mat'), 'globalCBF_all')

%% 2. summarize outlier cleaning information 
% global cbf
fprintf('\n\nSummarizing outlier cleaning.\n');

idx = 1;
IDList = [];
TPList = [];
IDTPList = [];
ocstat_temp = [];

for ii = 1:nsubj
    for jj = 1:ntp
        
        % init
        id = IDs{ii};
        tp = TimePoint{jj};
        
        if exist(fullfile(DATADIR, id, tp), 'dir')
        
            % ID information
            IDList{idx,1} = id;
            TPList{idx,1} = num2str(jj);
            IDTPList{idx,1} = [id, '_', num2str(jj)];
            
            for kk = 1:nruns
                
                % initialization
                run = RUNs{kk};
                run_dir = fullfile(DATADIR, id, tp, run);
                
                if exist(run_dir, 'dir')
                
                    clear stat
                    ocstat = fullfile(run_dir, 'clean_SCORE', 'stat.mat');
                    load(ocstat);
                    nmeasures = length(stat.measures);
                    start = nmeasures*(kk-1)+1;
                    ocstat_temp(idx, start:(start+nmeasures-1)) = stat.measures;
                    
                    
                else
                    
                    start = nmeasures*(kk-1)+1;
                    ocstat_temp(idx, start:(start+nmeasures-1)) = nan;
                    fprintf('%s %s does not have %s.\n',...
                        id, tp, run);
                end
                
            end
            
            % next iteration
            idx = idx + 1;
            
        end
    end    
end

% make array names
ARRAYNAME = [];
for kk = 1:nruns
    for ll = 1:nmeasures
        ARRAYNAME = [ARRAYNAME, ',''', stat.names{ll}, '_', RUNs{kk}, ''''];               
    end
end

% construct dataset
eval(['ocstat_all = dataset({IDList, ''ID''},', ...
      '{TPList,''Time''},{IDTPList,''ID_TIME''},', ...
      '{ocstat_temp', ARRAYNAME, '});']);

% save summarized global CBF
save(fullfile(EXPDIR, 'ocstat_all.mat'), 'ocstat_all')

%% generate rejection matrix for QC
tp = 1;
idx = 1;
ocstat_QC = ocstat_all(1, :);
for i = 1:size(ocstat_all,1)
    if strcmp(ocstat_all.Time(i), num2str(tp))
    
        ocstat_QC(idx, :) = ocstat_all(i, :);
        idx = idx + 1;
        
    end
end

%% 3. extract AAL ROI CBF
% global cbf
fprintf('\n\nSummarizing AAL ROI CBF.\n');

idx = 1;
IDList = [];
TPList = [];
IDTPList = [];
AALROICBF_temp = [];

for ii = 1:nsubj
    for jj = 1:ntp
        
        % init
        id = IDs{ii};
        tp = TimePoint{jj};
        
        if exist(fullfile(DATADIR, id, tp), 'dir')
        
            % ID information
            IDList{idx,1} = id;
            TPList{idx,1} = num2str(jj);
            IDTPList{idx,1} = [id, '_', num2str(jj)];
            
            for kk = 1:nruns
                
                % initialization
                run = RUNs{kk};
                run_dir = fullfile(DATADIR, id, tp, run);
                
                if exist(run_dir, 'dir')
                
                    clear AALROICBF
                    aalcbf = fullfile(run_dir, 'AALROI', 'AALROI.mat');
                    load(aalcbf);
                    nmeasures = length(AALROICBF.measures);
                    start = nmeasures*(kk-1)+1;
                    AALROICBF_temp(idx, start:(start+nmeasures-1)) = AALROICBF.measures;
                    
                else
                    
                    start = nmeasures*(kk-1)+1;
                    AALROICBF_temp(idx, start:(start+nmeasures-1)) = nan;
                    fprintf('%s %s does not have %s.\n',...
                        id, tp, run);
                end
                
            end
            
            % next iteration
            idx = idx + 1;
            
        end
    end    
end

% make array names
ARRAYNAME = [];
for kk = 1:nruns
    for ll = 1:nmeasures
        ARRAYNAME = [ARRAYNAME, ',''', AALROICBF.names{ll}, '_', RUNs{kk}, ''''];               
    end
end

% construct dataset
eval(['AALROICBF_all = dataset({IDList, ''ID''},', ...
      '{TPList,''Time''},{IDTPList,''ID_TIME''},', ...
      '{AALROICBF_temp', ARRAYNAME, '});']);
    
% save summarized global CBF
save(fullfile(EXPDIR, 'AALROICBF_all.mat'), 'AALROICBF_all')

%% 4. print CBF maps into images for QC purpose
fprintf('\n\nPrint images for QC.\n');

% generate figure
scrsz=get(0,'Screensize');
h = figure('Position',[1 1 scrsz(4) scrsz(4)],'Numbertitle','off');

for ii = 1:nsubj
    for jj = 1:ntp
        
        id = IDs{ii};
        tp = TimePoint{jj};
        
        if exist(fullfile(DATADIR, id, tp), 'dir')
            
            for kk = 1:nruns
                
                % initialization
                run = RUNs{kk};
                run_dir = fullfile(DATADIR, id, tp, run);
                maskfile_subj = fullfile(DATADIR, id, tp, 'MPRAGE', 'segmentation', 'SmallBrainMask.nii.gz');
                maskfile_temp = fullfile(DATADIR, id, tp, 'MPRAGE', 'normalized', 'swSmallBrainMask.nii.gz');
                
                if exist(run_dir, 'dir')
                    
                    % absolute CBF in subject space
                    infile = fullfile(run_dir, 'clean_SCORE', 'cleaned_meanCBF.nii.gz');
                    outfile = fullfile(EXPDIR, 'images', 'absolute_subj', ...
                        [id, '_', num2str(jj), tp, '_', num2str(kk), run, '.jpg']);
                    makeimages(h, infile, outfile, [0, 150], 1, 0.5, maskfile_subj);
                    
                    % relative CBF in subject space
                    infile = fullfile(run_dir, 'clean_SCORE', 'rel_cleaned_meanCBF.nii.gz');
                    outfile = fullfile(EXPDIR, 'images', 'relative_subj', ...
                        [id, '_', num2str(jj), tp, '_', num2str(kk), run, '.jpg']);
                    %makeimages(h, infile, outfile, [], 1, 0.5, maskfile_subj);
                    
                    infile = fullfile(run_dir, 'normalized', 'swcleaned_meanCBF.nii.gz');
                    outfile = fullfile(EXPDIR, 'images', 'abs_template', ...
                        [id, '_', num2str(jj), tp, '_', num2str(kk), run, '.jpg']);
                    makeimages(h, infile, outfile, [0, 150], 1, 0.5, maskfile_temp);
                    
                    infile = fullfile(run_dir, 'normalized', 'rel_swcleaned_meanCBF.nii.gz');
                    outfile = fullfile(EXPDIR, 'images', 'rel_template', ...
                        [id, '_', num2str(jj), tp, '_', num2str(kk), run, '.jpg']);
                    %makeimages(h, infile, outfile, [], 1, 0.5, maskfile_temp);
                    
                else
                    fprintf('%s %s does not have %s.\n',...
                        id, tp, run);
                end
                
            end
        end
    end
end

close(h)

%% 5. copy cbf images (absolute, relative) to result directory
fprintf('\n\nPrint images for QC.\n');

for ii = 1:nsubj
    for jj = 1:ntp
        for kk = 1:nruns
            
            % initialization
            id = IDs{ii};
            tp = TimePoint{jj};
            run = RUNs{kk};
            run_dir = fullfile(DATADIR, id, tp, run);
            
            if exist(run_dir, 'dir')
            
                % out_dir
                out_dir = fullfile(EXPDIR, 'CBFmaps', id);
                if ~exist(out_dir, 'dir')
                    mkdir(out_dir)
                end

                % copy absolute cbf in subject space
                infile = fullfile(run_dir, 'clean_SCORE', 'cleaned_meanCBF.nii.gz');
                outfile = fullfile(out_dir, ...
                    [id, '_', tp, '_', run, '_abs_subj.nii.gz']);
                copyfile(infile, outfile, 'f');
            
                % copy relative cbf in subject space
                infile = fullfile(run_dir, 'clean_SCORE', 'rel_cleaned_meanCBF.nii.gz');
                outfile = fullfile(out_dir, ...
                    [id, '_', tp, '_', run, '_rel_subj.nii.gz']);
                copyfile(infile, outfile, 'f');
                
                % copy absolute cbf in template space
                infile = fullfile(run_dir, 'normalized', 'swcleaned_meanCBF.nii.gz');
                outfile = fullfile(out_dir, ...
                    [id, '_', tp, '_', run, '_abs_template.nii.gz']);
                copyfile(infile, outfile, 'f');
                
                % copy relative cbf in template space
                infile = fullfile(run_dir, 'normalized', 'rel_swcleaned_meanCBF.nii.gz');
                outfile = fullfile(out_dir, ...
                    [id, '_', tp, '_', run, '_rel_template.nii.gz']);
                copyfile(infile, outfile, 'f');
                
            else
                 fprintf('%s %s does not have %s.\n',...
                    id, tp, run);
            end
            
        end
    end
end

%% 6. QC selection
QCmatrix_score = importQCmatrix_4groups(fullfile(EXPDIR, 'QCmatrix_samesubject.txt'));
QCmatrix = QCmatrix_score;
QCmatrix(QCmatrix_score <= 4) = 1;
QCmatrix(QCmatrix_score > 4) = 0;
QCmatrix(QCmatrix_score == -1) = -1;
mkdir(EXPDIR, 'QC')

% save QCmatrix to files
for ii = 1:ntp
    for jj = 1:nruns
        
        % initialization
        tp = TimePoint{ii};
        run = RUNs{jj};
        
        filename = fullfile(EXPDIR, 'QC', [tp, '_', run, '.txt']);
        dlmwrite(filename, QCmatrix(:, (ii-1)*nruns+jj));

    end
end

%% save QC result
idx = 1;
IDList = [];
TPList = [];
IDTPList = [];
qc_temp = [];

for ii = 1:nsubj
    for jj = 1:ntp
        
        % init
        id = IDs{ii};
        tp = TimePoint{jj};
        
        if exist(fullfile(DATADIR, id, tp), 'dir')
        
            % ID information
            IDList{idx,1} = id;
            TPList{idx,1} = num2str(jj);
            IDTPList{idx,1} = [id, '_', num2str(jj)];
            
            for kk = 1:nruns
                qc_temp(idx, kk) = QCmatrix(ii, (jj-1)*nruns+kk);
            end
            
            % next iteration
            idx = idx + 1;
            
        end
    end    
end

% make array names
ARRAYNAME = [];
for kk = 1:nruns
    ARRAYNAME = [ARRAYNAME, ',''', 'include', '_', RUNs{kk}, ''''];               
end

% construct dataset
eval(['qc_all = dataset({IDList, ''ID''},', ...
      '{TPList,''Time''},{IDTPList,''ID_TIME''},', ...
      '{qc_temp', ARRAYNAME, '});']);

save(fullfile(EXPDIR, 'QCmatrix.mat'), 'QCmatrix', 'qc_all')

%% 9. save all the information to csv file for statistical analysis
% global cbf
fprintf('\n\nSaving everything to csv.\n');
load(fullfile(EXPDIR, 'globalCBF_all.mat'));
load(fullfile(EXPDIR, 'ocstat_all.mat'));
load(fullfile(EXPDIR, 'QCmatrix.mat'));
load(fullfile(EXPDIR, 'AALROICBF_all.mat'));
%DEMOG = dataset('File', DEMOGRAPHIC_TXT, 'Delimiter', ',');

% INIT
PRECOLUME = 3;
N_qc = size(qc_all, 2) - PRECOLUME;
N_globalCBF = size(globalCBF_all, 2) - PRECOLUME;
N_AALROICBF = size(AALROICBF_all, 2) - PRECOLUME;
N_ocstat = size(ocstat_all, 2) - PRECOLUME;

% join the information
for kk = 1:nruns
    
    start = PRECOLUME+(N_qc/nruns)*(kk-1);
    qc_temp = qc_all(:, [1:PRECOLUME, start+1:(start+N_qc/nruns)]);
    start = PRECOLUME+(N_globalCBF/nruns)*(kk-1);
    globalCBF_temp = globalCBF_all(:, [1:PRECOLUME, start+1:(start+N_globalCBF/nruns)]);
    start = PRECOLUME+(N_AALROICBF/nruns)*(kk-1);
    AALROICBF_temp = AALROICBF_all(:, [1:PRECOLUME, start+1:(start+N_AALROICBF/nruns)]);
    start = PRECOLUME+(N_ocstat/nruns)*(kk-1);
    ocstat_temp = ocstat_all(:, [1:PRECOLUME, start+1:(start+N_ocstat/nruns)]);
    
    fullinfo_temp = join(qc_temp, globalCBF_temp);
    fullinfo_temp = join(fullinfo_temp, AALROICBF_temp);
    fullinfo_temp = join(fullinfo_temp, ocstat_temp);
    
    if kk == 1
        fullinfo = fullinfo_temp;
    else
        fullinfo = join(fullinfo, fullinfo_temp);
    end

end
%fullinfo = join(fullinfo, DEMOG(:,3:end));

% write to file
filename = fullfile(EXPDIR, 'allinfo.csv');
export(fullinfo, 'File', filename, 'Delimiter', ',');
save(fullfile(EXPDIR, 'allinfo.mat'), 'fullinfo');

fprintf('done!\n')


