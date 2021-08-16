function mc_SPM8_session(session_dir, infile_gz, out_folder, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 4
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
global defaults;
spm_get_defaults;

% Get realignment defaults
defs = defaults.realign;

% Flags to pass to routine to calculate realignment parameters
% (spm_realign)

%as (possibly) seen at spm_realign_ui,
% -fwhm = 5 for fMRI
% -rtm = 0 for fMRI
% for this particular data set, we did not perform a second realignment to the mean
% the coregistration between the reference control and label volume is also omitted
reaFlags = struct(...
    'quality', defs.estimate.quality,...  % estimation quality
    'fwhm', 5,...                         % smooth before calculation
    'rtm', 1,...                          % whether to realign to mean
    'PW',''...                            %
    );

% Flags to pass to routine to create resliced images
% (spm_reslice)
resFlags = struct(...
    'interp', 1,...                       % trilinear interpolation
    'wrap', defs.write.wrap,...           % wrapping info (ignore...)
    'mask', defs.write.mask,...           % masking (see spm_reslice)
    'which',2,...                         % write reslice time series for later use
    'mean',1);                            % do write mean image

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
    fprintf('Motion correction for run %s (%0.0f/%0.0f).\n', ...
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
    
    %% perform motion correction
    
    P=spm_select('ExtFPlist', out_dir, infile_filename, Inf);
    spm_realign_asl(P, reaFlags);
    
    % Run reslice
    spm_reslice(P, resFlags);
    
    %% check motion parameter
    % mc param file
%     [~, filename] = fileparts(infile_filename);
%     parfile = fullfile(out_dir, ['rp_', filename, '.txt']);
%     fid = fopen(parfile,'r');
%     mcparam = (fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g\r\n',[12 inf]))';
%     fclose(fid);
% 
%     hmc=figure(); 
%     clf
%     subplot(211),
%     plot(mcparam(:,[1:3,7:9]),'--','LineWidth',2.5),
%     subplot(212),plot(mcparam(:,[4:6,10:12]),'--','LineWidth',2.5)
% 
%     % Save image to disk
%     realignfig = fullfile(out_dir,[filename '.png']);
%     defs.printstr = ['print -dpng -painters -append -noui ' realignfig];
%     eval(defs.printstr);
%     close(hmc);
    
    %% remove the original unzip file and zip the outputfile
    cd(out_dir)
    delete(infile_filename);
    system('gzip *.nii');
    
end

% go back to original directory
cd(savecurpath);
