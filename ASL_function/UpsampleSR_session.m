function UpsampleSR_session(session_dir, infile_gz, out_folder, out_prefix, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
SRDIR = '/home/longxie/pkg/superresolution';
NLMUpsampleDir = fullfile(SRDIR, 'NLMUpsample');
addpath(ASLFUNCDIR)
addpath(NLMUpsampleDir)

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
lf_x = 2;
lf_y = 2;
lf_z = 2;

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

%% Run MCFLIRT
savecurpath = pwd;

for r = 1:nruns   
    
    %% initialization
    fprintf('Upsample for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % generate output folder
    out_dir = fullfile(run_dir, out_folder);
    if exist(out_dir, 'dir')
        rmdir(out_dir,'s');
    end
    mkdir(out_dir);
    
    %% check and unzip input files
    cd(run_dir);
    
    % check if file exist
    if exist(infile_gz, 'file')
        % unzip file for SPM usage
        infile = fullfile(out_dir, infile_filename);
        system(sprintf('gunzip -c %s > %s', infile_gz, infile));
        pause(0.5)
    else
        % if not exist, report error
        msg = sprintf('ERROR: %s does not exist in %s.\n', infile_gz, run_dir);
        system(sprintf('echo "%s" >> %s', msg, LOGTXT));
        error(msg);
    end
    
    %% split the 4D file to 3D
    V = spm_select('FPlist', out_dir, infile_filename);
    Vo = spm_file_split(V, out_dir);
    
    %% perform superresolution for each measurement
    for i = 1:length(Vo)
        
        fprintf('Running frame %s/%s.\n', ...
            num2str(i), num2str(length(Vo)));
        
        % filenames
        file3D = Vo(i).fname;
        [dir, filename3D, ext] = fileparts(file3D);
        out_file3D = fullfile(dir, [out_prefix, filename3D, ext]);
        
        % use NLMUpsample
        NLMUpsample2(file3D, out_file3D, lf_x, lf_y, lf_z);
        
        % get rid of negative value
        cmd = sprintf('c3d %s -clip 0 inf -o %s', out_file3D, out_file3D);
        system(cmd);
        
    end
    
    %% merge upsampled 3d files
    files = spm_select('FPlist', out_dir, ['^', out_prefix, '\w*.*nii$']);
    out_file4D = fullfile(out_dir, [out_prefix, infile_filename]);
    spm_file_merge(cellstr(files), out_file4D, 4);
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir)
    system(sprintf('gzip %s', out_file4D));
    delete(fullfile(out_dir, '*.nii'))
    
end

% go back to original directory
cd(savecurpath);
