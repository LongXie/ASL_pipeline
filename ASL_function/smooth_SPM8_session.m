function smooth_SPM8_session(session_dir, infile_gz, out_folder, out_prefix, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 5
    error('Not enough input argument!');
end

% check input extension
[~, infile_filename, ext] = fileparts(infile_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(infile_filename);
    if ~strcmp(ext, '.nii')
        error('Input file`s extension must be .nii.gz');
    end
else
    error('Input file`s extension must be .nii.gz');
end

%% parameters
SmoothFWHM = 4;

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

%% Run 
savecurpath = pwd;

for r = 1:nruns
    
    %% initialization
    fprintf('Smooth for run %s (%0.0f/%0.0f).\n', ...
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
    
    % input file
    %infile_gz  = fullfile('coreg', 'crRAWPCASL.nii.gz');
    
    % check if file exist
    if exist(infile_gz, 'file')
        % unzip file for SPM usage
        %[~, infile_filename] = fileparts(infile_gz);
        infile = fullfile(out_dir, infile_filename);
        system(sprintf('gunzip -c %s > %s', infile_gz, infile));
        pause(0.5)
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s does not exist in %s.\n', infile_gz, run_dir);
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end
    
    %% perform smooth
    % select files
    P = spm_select('ExtFPList', out_dir, infile_filename, Inf);
    
    %take image per image smooth
    for im = 1:size(P,1)
        %this is actual image
        Pim = P(im,:);
        %this the new image name
        Qim = fullfile(out_dir, [out_prefix, spm_str_manip(Pim,'dt')]);
        %now call spm_smooth with kernel defined at PAR
        spm_smooth(Pim, Qim, SmoothFWHM);
    end 
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    system('gzip *.nii');
    
end

% go back to original folder
cd(savecurpath);
