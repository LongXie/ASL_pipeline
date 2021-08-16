function check_mc_SPM8_session(session_dir, mc_param)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 2
    error('Not enough input argument!');
end

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


%% show parameter
for r = 1:nruns
    
    %% initialization
    fprintf('Check motion correction for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % mc param file
    parfile = fullfile(run_dir, mc_param);
    fid = fopen(parfile,'r');
    mcparam = (fscanf(fid, '%g %g %g %g %g %g %g %g %g %g %g %g\r\n',[12 inf]))';
    fclose(fid);
    
    %%
    hmc=figure('visible', 'off'); 
    %hmc=figure(); 
    clf
    subplot(211),
    plot(mcparam(:,[1:3,7:9]),'--','LineWidth',2.5),
    subplot(212),plot(mcparam(:,[4:6,10:12]),'--','LineWidth',2.5)
    pause(0.5)
        
    %% Save image to disk
    [outdir, filename, ext] = fileparts(parfile);
    realignfig = fullfile(outdir,[filename '.png']);
    defs.printstr = ['print -dpng -painters -append -noui ' realignfig];
    eval(defs.printstr);
    close(hmc);
end

end