function R= RobustBayesianASLdenoising(P,roiloc,maskloc)

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

v=spm_vol(P);
dat=spm_read_vols(v);
[x,y,z,~]=size(dat);

roidat=spm_read_vols(spm_vol(roiloc));
gm=roidat(:,:,:,1);
wm=roidat(:,:,:,2);
csf=roidat(:,:,:,3);

mask=spm_read_vols(spm_vol(maskloc))>0;
mask=mask&(mean(dat,4)~=0);
gm=gm.*mask;
wm=wm.*mask;
csf=csf.*mask;


[Score,dat]=SCORE_den(dat,gm,wm,csf,0.95);
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


function [recon,dat]=SCORE_den(dat,gmtpm,wmtpm,csftpm,thresh)

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









%
%
% w=3;
% [localprior,localvar]=deal(zeros(size(Score)));
% for xdim=1:x
%     for ydim=1:y
%         for zdim=1:z
%             if mask(xdim,ydim,zdim)==0
%                 continue;
%             end
%             wt=0;
%             tmplocalprior=0;
%             tmpvar=0;
%             for windx=-w:w
%                 if (xdim+windx)<1||(xdim+windx)>x
%                     continue;
%                 end
%                 for windy=-w:w
%                     if (ydim+windy)<1||(ydim+windy)>y
%                         continue;
%                     end
% 
%                     if mask(xdim+windx,ydim+windy,zdim)==0
%                         continue;
%                     end
%                     wt1=exp(-0.5*((gm(xdim,ydim,zdim)-gm(xdim+windx,ydim+windy,zdim))^2+...
%                         (wm(xdim,ydim,zdim)-wm(xdim+windx,ydim+windy,zdim))^2+...
%                         (csf(xdim,ydim,zdim)-csf(xdim+windx,ydim+windy,zdim))^2));
% 
%                     wt2=exp(-(windx^2+windy^2)/18);
%                     wt3=exp(-V(xdim+windx,ydim+windy,zdim)/globalvar);
% 
%                     tmpwt=wt1*wt2*wt3;
%                     wt=wt+tmpwt;
% 
% 
%                     tmplocalprior=tmplocalprior+tmpwt*Score(xdim+windx,ydim+windy,zdim);
%                     tmpvar=tmpvar+tmpwt*V(xdim+windx,ydim+windy,zdim);
% 
% 
%                 end
%             end
%             localprior(xdim,ydim,zdim)=tmplocalprior/wt;
%             localvar(xdim,ydim,zdim)=tmpvar/wt;
% 
%         end
%     end
% end



%         b1=myrobustfit(Y',[],[],exp(-(Y-repmat(Prior,1,size(dat,4))).^2/2/globalstd^2)',(S'/globalstd).^3,Prior');


% mu=(harmmean([localvar(mask)';V(mask)'])/globalvar).^2;
% lmd=(V(mask)'./(localvar(mask)'+(localvar(mask)'/globalvar).^2)).^2;
%
%
%
%
%
%
% L=zeros(x,y,z);L(mask)=lmd;
% M=zeros(size(Score)); M(mask)=mu;


