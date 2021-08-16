function [recon,noimg,indx]=...
    SCORE_core(cbfloc,gmloc,wmloc,csfloc,thresh, range)

if nargin < 6
    range = 0;
end

dat=(spm_read_vols(spm_vol(cbfloc)));
gm=(spm_read_vols(spm_vol(gmloc)));
wm=(spm_read_vols(spm_vol(wmloc)));
csf=(spm_read_vols(spm_vol(csfloc)));

msk=(gm+wm+csf)>0;
[gm,wm,csf]=deal(gm./(gm+wm+csf),wm./(gm+wm+csf),csf./(gm+wm+csf));
[gm(~msk),wm(~msk),csf(~msk)]=deal(0);

level1=csf>thresh;                        
level2=wm>thresh;
level3=gm>thresh;

TD=size(dat,4);
indx=zeros(TD,1);
if range == 1
    for tdim=1:TD
        tmp=dat(:,:,:,tdim);
        if (mean(tmp(level3))<0)||(mean(tmp(level3))>150)
            mean(tmp(level3));
            indx(tdim)=1;
        end
        
    end
end

A=mean(dat(:,:,:,(indx==0)),4);
%Crr=sqrt((std(A(level1))/mean(A(level1)))^2+(std(A(level2))/mean(A(level2)))^2+(std(A(level3))/mean(A(level3)))^2);
%Crr=var(A(level1))+var(A(level2))+var(A(level3));
Crr=sum(level1(:))*var(A(level1))+sum(level2(:))*var(A(level2))+sum(level3(:))*var(A(level3));

flag=0;
while flag==0
    prevcorr=Crr;
    prevA=A;
    previndx=indx;
    
    CC=-2*ones(TD,1);
    
    for tdim=1:TD
        if indx(tdim)==0            
            tmp=dat(:,:,:,tdim);
            CC(tdim)=corr(A(msk),tmp(msk));
        end
    end
    [~,in]=max(CC);
    
    indx(in)=2;
    
    A=mean(dat(:,:,:,(indx==0)),4);
    %             Crr=sqrt((std(A(level1))/mean(A(level1)))^2+(std(A(level2))/mean(A(level2)))^2+(std(A(level3))/mean(A(level3)))^2);
    Crr=sum(level1(:))*var(A(level1))+sum(level2(:))*var(A(level2))+sum(level3(:))*var(A(level3));
    if (Crr>prevcorr)||(sum(indx==0)==0)
        flag=1;
        A=prevA;
        indx=previndx;

    end
end

recon=A;
noimg=sum(indx==0);

