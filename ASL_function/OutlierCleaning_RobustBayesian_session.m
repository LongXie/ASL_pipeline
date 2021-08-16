function OutlierCleaning_RobustBayesian_session...
    (session_dir, anat_dir, infile_gz, brainmask_gz, out_folder, out_filename_gz, out_stat, LOGTXT)

ASLFUNCDIR = fullfile('/data/picsl/longxie/WolkMCI', 'code', 'ASL_function');
addpath(ASLFUNCDIR)

%% check input and output
if nargin < 7
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

% check output extension
[~, out_filename, ext] = fileparts(out_filename_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(out_filename);
    if ~strcmp(ext, '.nii')
        error('Output file`s extension must be .nii.gz');
    end
else
    error('Output file`s extension must be .nii.gz');
end

% check output data file
[~, ~, ext] = fileparts(out_stat);
if ~strcmp(ext, '.mat')
    error('Output data file`s extension must be .mat');
end

% check mask extension
[~, brainmask_filename, ext] = fileparts(brainmask_gz);
if strcmp(ext, '.gz')
    [~, ~, ext] = fileparts(brainmask_filename);
    if ~strcmp(ext, '.nii')
        error('Mask file`s extension must be .nii.gz');
    end
else
    error('Mask file`s extension must be .nii.gz');
end

%% parameters
thresh = 0.95;

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

%% check and unzip tissue probability maps
%anat_dir = fullfile(session_dir, 'MPRAGE', 'segmentation');

% tissue probability mapes
gm_prob_gz = spm_select('FPList', anat_dir, ['^c1\w*.*nii.gz']);
wm_prob_gz = spm_select('FPList', anat_dir, ['^c2\w*.*nii.gz']);
csf_prob_gz = spm_select('FPList', anat_dir, ['^c3\w*.*nii.gz']);
%brainmask_gz = fullfile(anat_dir, 'BrainMask.nii.gz');

% check if exist
if exist(gm_prob_gz, 'file') && exist(wm_prob_gz, 'file') ...
    && exist(csf_prob_gz, 'file')
    % unzip file for SPM usage
    % gm
    [~, gm_filename] = fileparts(gm_prob_gz);
    gm_prob = fullfile(anat_dir, gm_filename);
    system(sprintf('gunzip -c %s > %s', gm_prob_gz, gm_prob));
    % wm
    [~, wm_filename] = fileparts(wm_prob_gz);
    wm_prob = fullfile(anat_dir, wm_filename);
    system(sprintf('gunzip -c %s > %s', wm_prob_gz, wm_prob));
    % csf
    [~, csf_filename] = fileparts(csf_prob_gz);
    csf_prob = fullfile(anat_dir, csf_filename);
    system(sprintf('gunzip -c %s > %s', csf_prob_gz, csf_prob));
    % brain mask
    brainmask = fullfile(anat_dir, brainmask_filename);
    system(sprintf('gunzip -c %s > %s', brainmask_gz, brainmask));
else
    % if not exist, report error
    msg = sprintf('ERROR: Tissue probability maps or brain mask do not exist in %s.\n', session_dir);
    system(sprintf('echo "%s" >> %s', msg, LOGTXT));
    error(msg);
end

%% Run 
savecurpath = pwd;
for r = 1:nruns
    
    %% initialization
    fprintf('Clean CBF using RoobustBayesian SCORE for run %s (%0.0f/%0.0f).\n', ...
        d{r}, r, nruns);
    run_dir = fullfile(session_dir, d{r});
    
    % generate output folder
    out_dir = fullfile(run_dir, out_folder);
    if isdir(out_dir)
        rmdir(out_dir,'s');
    end
    mkdir(out_dir);
    
    %% check and unzip CBF time series
    cd(run_dir);
    
    % input file
    %infile_gz  = fullfile('perf', 'cbf_0_scrRAWPCASL_4D.nii.gz');
    
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
    
    %% perform cleaning
    % load time series
    cbfloc = spm_select('ExtFPList', out_dir, ['^', infile_filename, '$'], Inf);
    
    % load tissue maps
    gmloc = spm_select('FPlist', anat_dir, ['^', gm_filename, '$']);
    wmloc = spm_select('FPlist', anat_dir, ['^', wm_filename, '$']);
    csfloc = spm_select('FPlist', anat_dir, ['^', csf_filename, '$']);
    
    % SCORE
    [recon, oidx] = RobustBayesianASLdenoising(cbfloc, gmloc, wmloc, csfloc, brainmask, thresh);
    %recon(isnan(recon(:))) = 0;
    
    % get statistics
    stat = [];
    stat.RejectRateTotal = mean(oidx>=1);
    stat.RejectRateS1 = mean(oidx==1);
    stat.RejectRateS2 = mean(oidx==2);
    stat.TotalPairs = length(oidx);
    stat.RemainPairs = sum(oidx == 0);
    stat.odix = oidx;
    stat.names = {'RejectRateTotal', ...
                  'TotalPairs', ...
                  'RemainPairs'};
    stat.measures = [stat.RejectRateTotal, ...
                     stat.TotalPairs, ...
                     stat.RemainPairs];

    % save the reconstructed image
    vo = spm_vol(cbfloc);
    vo = vo(1);
    %out_filename = 'cleaned_meanCBF.nii';
    vo.fname = fullfile(out_dir, out_filename);
    spm_write_vol(vo, recon);

    % save statistics
    stat_fn = fullfile(out_dir, out_stat);
    save(stat_fn, 'stat')

    %% Output status
    fprintf('....Rejection rate is %1.2f. %1.0f out of %1.0f remain.).\n', ...
        stat.RejectRateTotal, stat.RemainPairs, stat.TotalPairs);

    %% remove the original unzip file and zip the outputfile
    cd(out_dir);
    delete(infile_filename);
    system('gzip -f *.nii');
    
end

% remove the unzip anat files
cd(anat_dir)
delete(gm_filename, wm_filename, csf_filename, brainmask);

% back to original dir
cd(savecurpath);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = myrobustfit(y,priorw,mu,bprior,lmd,localprior)


% if ~exist('wfun','var')||isempty(wfun)
%     wfun = @bisquare;
tune = 4.685;
% end


if ~exist('lmd','var')||isempty(lmd)
    lmd=0;
    bprior=1;
end

% if (nargin<3 || isempty(wfun)), wfun = 'bisquare'; end


[n,p]=size(y);

if ~exist('priorw','var')||isempty(priorw)
    priorw=ones(size(y));
else
    priorw=priorw./repmat(max(priorw),n,1);
end



X=ones(n,p);

if ~all(priorw==1)
    sw = sqrt(priorw);
    X = bsxfun(@times,X,sw);
    y = y.*sw;
else
    sw = 1;
end

b=(sum(X.*y)+mu.*bprior+lmd.*localprior)./(sum(X.*X)+mu+lmd);%X\y;

% b=(     (sum(X.*y)+lmd.*localprior)./(sum(X.*X)+lmd)  +    mu.*bprior  )./(1+mu);%X\y;


b0 = zeros(size(b));
h = min(.9999, (X./  repmat(sqrt(sum(X.^2)),n,1)     ).^2);
adjfactor = 1 ./ sqrt(1-h./priorw);
% dfe = n-xrank;




% If we get a perfect or near perfect fit, the whole idea of finding
% outliers by comparing them to the residual standard deviation becomes
% difficult.  We'll deal with that by never allowing our estimate of the
% standard deviation of the error term to get below a value that is a small
% fraction of the standard deviation of the raw response values.
tiny_s = 1e-6 * std(y);
tiny_s(tiny_s==0)=1;

D = sqrt(eps(class(X)));
iter = 0;
iterlim = 50;
while((iter==0) || any(abs(b-b0) > D*max(abs(b),abs(b0))))
    iter = iter+1;
    if (iter>iterlim)
        %       warning(message('stats:statrobustfit:IterationLimit'));
        break;
    end
    
    % Compute residuals from previous fit, then compute scale estimate
    r = y - X.*repmat(b,n,1);
    radj = r .* adjfactor ./ sw;
    
    rs=sort(abs(radj));
    
    s = median(rs)/0.6745;
    
    % Compute new weights from these residuals, then re-fit
    w = bisquare(radj./ repmat((max(s,tiny_s)*tune),n,1)   );
    b0 = b;
    %    [b,wxrank] = wfit(y,X,w);
    z=sqrt(w);
    x=X.*z;
    yz=y.*z;
    b=(sum(x.*yz)+mu.*bprior+lmd.*localprior)./(sum(x.*x)+mu+lmd);
end

function w = bisquare(r)
w = (abs(r)<1) .* (1 - r.^2).^2;

function s = madsigma(r)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
rs = sort(abs(r));
s = median(rs) / 0.6745;

function [recon,dat,indx]=SCORE_den(dat,gm,wm,csf,thresh)

% SCORE: performs Structural Correlation based Outlier REjection
% Input:
%   dat: CBF time series (4D)
%   gmtpm: Grey matter tissue probability map
%   wmtpm: White matter tissue probability map
%   csftpm: CSF tissue probability map
%   thresh: Threshold to create gm, wm and csf mask
% 
% Output:
%   recon: Estimated mean CBF map
%   noimg: Number of volumes retained
%   indx: Index file;  0: volumes retained for the mean CBF computation
%                      1: volumes rejected based on GM threshold
%                      2: volumes rejected based on structural correlation
% Report any bug to Sudipto Dolui: sudiptod@mail.med.upenn.edu
% Copyright Sudipto Dolui
                    
      
msk=(gm+wm+csf)>0;
gm=gm>thresh;
wm=wm>thresh;
csf=csf>thresh;


nogm=sum(gm(:)>0)-1;
nowm=sum(wm(:)>0)-1;
nocsf=sum(csf(:)>0)-1;


TD=size(dat,4);
MnGM=zeros(TD,1);



for tdim=1:TD
    tmp=dat(:,:,:,tdim);
    MnGM(tdim)=mean(tmp(gm));
end
MedMnGM=median(MnGM);           % Robust estimation of mean
SDMnGM=mad(MnGM,1)/0.675;       % Robust estimation of standard deviation
indx=double(abs(MnGM-MedMnGM)>2.5*SDMnGM);    % Volumes outside 2.5 SD of Mean are discarded
%indx=zeros(1,size(dat,4));

R=mean(dat(:,:,:,indx==0),4);
V=nogm*var(R(gm))+nowm*var(R(wm))+nocsf*var(R(csf));
V_prev=V+1;

while V<V_prev
    V_prev=V;
    R_prev=R;
    indx_prev=indx;
    CC=-2*ones(TD,1);
    for tdim=1:TD
        if(indx(tdim)~=0)
            continue;
        end
        tmp=dat(:,:,:,tdim);
        CC(tdim)=corr(R(msk),tmp(msk));
    end
    [~,inx]=max(CC);
    indx(inx)=2;
    R=mean(dat(:,:,:,indx==0),4);
    V=nogm*var(R(gm))+nowm*var(R(wm))+nocsf*var(R(csf));
end

recon=R_prev;
% indx=indx_prev;
% noimg=sum(indx==0);
dat=dat(:,:,:,indx_prev==0);

function [R, indx] = RobustBayesianASLdenoising(cbfloc,gmloc,wmloc,csfloc,maskloc,thresh)

% RobustBayesianASLdenoising: performs ASL denoising using a Robust
% Bayesian approach
% Input:
%   P: File name of 4D CBF time series
%   roiloc: file names of GM, WM and CSF tissue probability maps
%   maskloc: file name of brain mask
% 
% Output:
%   R: Estimated mean CBF map
% Report any bug to Sudipto Dolui: sudiptod@mail.med.upenn.edu

%v=spm_vol(P);
%dat=spm_read_vols(v);

dat=(spm_read_vols(spm_vol(cbfloc)));
gm=(spm_read_vols(spm_vol(gmloc)));
wm=(spm_read_vols(spm_vol(wmloc)));
csf=(spm_read_vols(spm_vol(csfloc)));
[gm,wm,csf]=deal(gm./(gm+wm+csf),wm./(gm+wm+csf),csf./(gm+wm+csf));
[x,y,z,~]=size(dat);



%roidat=spm_read_vols(spm_vol(roiloc));
%gm=roidat(:,:,:,1);
%wm=roidat(:,:,:,2);
%csf=roidat(:,:,:,3);

mask=spm_read_vols(spm_vol(maskloc))>0;
%mask=(gm+wm+csf)>0;
%mask=mask&(mean(dat,4)~=0);
gm=gm.*mask;
wm=wm.*mask;
csf=csf.*mask;


[Score,dat,indx]=SCORE_den(dat,gm,wm,csf,thresh);
clear T

Y=zeros(sum(mask(:)),size(dat,4));

for k=1:size(dat,4)
    tmp=dat(:,:,:,k);
    Y(:,k)=tmp(mask);
end

R=Score;

for iter=1
    GMCBF=median(R(gm>0.9));
    WMCBF=median(R(wm>0.9));
    % CSFCBF=median(R(csf>0.9));
    GlobalPriorfull=GMCBF*gm+WMCBF*wm;%+CSFCBF*csf;
    GlobalPrior=GlobalPriorfull(mask);
    VV=(mad(Y,1,2)/0.6745).^2;
    globalvar=(mad(Y(:),1)/0.6745)^2;
%     V=zeros(size(Score));
%     V(mask)=VV;


    mu=(VV'/globalvar);
    mu=(mu<1).*mu.^4+(mu>=1).*(4*mu-3);%(mu<1).*(mu.^4)+(mu>=1).*(2*mu.^2-1);
    localprior=0;
    lmd=0;

    b1=myrobustfit(Y',[],mu,GlobalPrior',lmd,localprior);


    R=zeros(x,y,z);
    R(mask)=b1;
end
