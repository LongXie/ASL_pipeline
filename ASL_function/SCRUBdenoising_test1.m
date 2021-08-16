function [prefix,R,Score,indx]= SCRUBdenoising_test1(filename,gmloc,wmloc,csfloc,maskimg,wfun,flagprior,flagmodrobust,flagstd,prefix)

% RobustBayesianASLdenoising: performs ASL denoising using a Robust
% Bayesian approach
% Input:
%   filename: File name of 4D CBF time series
%   tpmimgs: file names of GM, WM and CSF tissue probability maps in order
%   maskimg: file name of brain mask
%   wfun: robust function used in robust regression
%   flagprior: 1 implies use prior, 0 no prior
%   flagmodrobust: Use modified robust method
%   flagstd: 1 implies std in robust regression computed using regular std,
%   0 implies robust estimate
%   prefix: prefix to be added to the name of the saved file
%
% Output:
%   prefix: prefix added to the name of the file
%   R: Estimated mean CBF map
%   Score: Score output
% Report any bug to Sudipto Dolui: sudiptod@mail.med.upenn.edu


if ~ismember(flagstd,[0,1])
    flagstd=1;       % flagstd = 1 implies std, 0 implies robust method
end
if ~ismember(flagprior,[0,1])
    flagprior=1;       % flagprior = 1 implies global prior, 0 implies no prior
end
if ~ismember(flagmodrobust,[0,1])
    flagprior=1;       % flagprior = 1 implies global prior, 0 implies no prior
end


v=spm_vol(filename);
dat=spm_read_vols(v);
[x,y,z,~]=size(dat);

%roidat=spm_read_vols(spm_vol(tpmimgs));
%gm=roidat(:,:,:,1);
%wm=roidat(:,:,:,2);
%csf=roidat(:,:,:,3);
gm=spm_read_vols(spm_vol(gmloc));
wm=spm_read_vols(spm_vol(wmloc));
csf=spm_read_vols(spm_vol(csfloc));

mask=spm_read_vols(spm_vol(maskimg))>0;
mask=mask & (max(abs(dat),[],4)>0) & any(~isnan(dat),4);
gm=gm.*mask;
wm=wm.*mask;
csf=csf.*mask;


[Score,dat,indx]=SCORE_den(dat,gm,wm,csf,0.95);
clear T

Y=zeros(sum(mask(:)),size(dat,4));

for k=1:size(dat,4)
    tmp=dat(:,:,:,k);
    Y(:,k)=tmp(mask);
end

R=Score;

for iter=1
    if flagprior==0
        mu=0;
        GlobalPrior=0;
    else
        GMCBF=median(R(gm>0.9));
        WMCBF=median(R(wm>0.9));
        % CSFCBF=median(R(csf>0.9));
        GlobalPriorfull=GMCBF*gm+WMCBF*wm;%+CSFCBF*csf;
        GlobalPrior=GlobalPriorfull(mask);
        VV=(mad(Y,1,2)/0.6745).^2;
        % %     globalvar=(mad(Y(:),1)/0.6745)^2;
        % %     globalvar=(median(VV)+3*mad(VV,1)/0.6745);
        %      globalvar=(median(VV));
        % %     V=zeros(size(Score));
        % %     V(mask)=VV;
        
        globalvar=mean(var(Y,[],2))/size(Y,2)+var(mean(Y,2));
        
        mu1=(VV'/globalvar);
        mu=max(0,mu1-mean(mu1));
        % %     mu=((mu>=1)&(mu<10)).*(mu-1)+(mu>=10).*(0.05*mu.^2+4);% +  ((R(mask))'<0).*abs(R(mask))';
        %     mu=((mu>=3)&(mu<10)).*(mu-3)+(mu>=10).*(0.05*mu.^2+2);% +  ((R(mask))'<0).*abs(R(mask))';
        % %     mu=mu*size(Y,2)/9;
    end
    localprior=0;
    lmd=0;
    
    
    b1=myrobustfit(Y',[],mu,GlobalPrior',lmd,localprior,wfun,flagstd,flagmodrobust);
    
    
    R=zeros(x,y,z);
    R(mask)=b1;
end

if ~exist('prefix','var')||isempty(prefix)
    prefix = [];
    if flagprior==1
        prefix = strcat(prefix,'prior_');
    else
        prefix = strcat(prefix,'noprior_');
    end
    
    if flagmodrobust==1
        prefix = strcat(prefix,'modrob_');
    else
        prefix = strcat(prefix,'rob_');
    end
    
    prefix = strcat(prefix,wfun,'_');
    
    if flagstd==1
        prefix = strcat(prefix,'std');
    else
        prefix = strcat(prefix,'rstd');
    end
end
% if prefix(end)~='_'
%     prefix=strcat(prefix,'_');
% end

%vo=v(1);
%vo.fname=fullfile(spm_str_manip(vo.fname,'H'),[prefix '_' spm_str_manip(vo.fname,'t')]);
%spm_write_vol(vo,R);



function b = myrobustfit(y,priorw,mu,bprior,lmd,localprior,wfun,flagstd,flagmodrobust)


if ~exist('wfun','var')||isempty(wfun)
    wfun = 'huber';
end


if ischar(wfun)
    switch(wfun)
        case 'andrews'
            wfun = @andrews;
            t = 1.339;
        case 'bisquare'
            wfun = @bisquare;
            t = 4.685;
        case 'cauchy'
            wfun = @cauchy;
            t= 2.385;
        case 'fair'
            wfun = @fair;
            t = 1.400;
        case 'huber'
            wfun = @huber;
            t = 1.345;
        case 'logistic'
            wfun = @logistic;
            t = 1.205;
        case 'ols'
            wfun = @ols;
            t = 1;
        case 'talwar'
            wfun = @talwar;
            t = 2.795;
        case 'welsch'
            wfun = @welsch;
            t = 2.985;
    end
end
tune = t;



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

b=(sum(X.*y)+mu.*bprior+lmd.*localprior)./(sum(X.*X)+mu+lmd);

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
    
    if flagstd==1
        s = sqrt(mean(radj.^2));
    elseif flagstd==0
        rs=sort(abs(radj));
        s = median(rs)/0.6745;
    else
        error('Unknown flagstd')
    end
    
    % Compute new weights from these residuals, then re-fit
    w = wfun(radj.*(1-flagmodrobust*exp(-repmat(mu,n,1)))./ repmat((max(s,tiny_s)*tune),n,1)   );
    b0 = b;
    %    [b,wxrank] = wfit(y,X,w);
    z=sqrt(w);
    x=X.*z;
    yz=y.*z;
    b=(sum(x.*yz)+mu.*bprior+lmd.*localprior)./(sum(x.*x)+mu+lmd);
end



% function w = bisquare(r)
% w = (abs(r)<1) .* (1 - r.^2).^2;

% function w = huber(r)
% w = 1 ./ max(1, abs(r));

% function w = talwar(r)
% w = 1*(abs(r)<1);

function w = andrews(r)
r = max(sqrt(eps(class(r))), abs(r));
w = (abs(r)<pi) .* sin(r) ./ r;

function w = bisquare(r)
w = (abs(r)<1) .* (1 - r.^2).^2;

function w = cauchy(r)
w = 1 ./ (1 + r.^2);

function w = fair(r)
w = 1 ./ (1 + abs(r));

function w = huber(r)
w = 1 ./ max(1, abs(r));

function w = logistic(r)
r = max(sqrt(eps(class(r))), abs(r));
w = tanh(r) ./ r;

function w = ols(r)
w = ones(size(r));

function w = talwar(r)
w = 1 * (abs(r)<1);

function w = welsch(r)
w = exp(-(r.^2));


function s = madsigma(r)
%MADSIGMA    Compute sigma estimate using MAD of residuals from 0
rs = sort(abs(r));
s = median(rs) / 0.6745;


function [recon,dat, indx]=SCORE_den(dat,gmtpm,wmtpm,csftpm,thresh)

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



gm=gmtpm>thresh;
wm=wmtpm>thresh;
csf=csftpm>thresh;
msk=(gm+wm+csf)>0;

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



