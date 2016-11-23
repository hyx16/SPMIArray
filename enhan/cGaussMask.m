% calculate masks for time-frequency bin using complex Gauss model
% corresponding to 'ME" in  "The THU-SPMI CHiME-4 system : Lightweight design with
% advanced multi-channel processing, feature enhancement, and language modeling."

% implemented based on Higuchi T, Ito N, Yoshioka T, et al. 
% "Robust MVDR beamforming using time-frequency masks for online/offline ASR in noise[C]"
% IEEE International Conference on Acoustics, Speech and Signal Processing. 2016.

% Copyright 2016 Tsinghua University SPMI Lab, Hongyu Xiang

% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

function softmask = cGaussMask(ftbin,Nsource,XX,EMITERNUM)

[Nchan,Nbin,Nframe] = size(ftbin);
%EM algorithm
softmask = zeros(Nbin, Nframe, Nsource);
%init parameters
phiCor = ones(Nbin, Nframe, Nsource);
Rcor = zeros(Nchan, Nchan, Nbin, Nsource);
Rcor(:,:,:,1) = mean(XX, 4);% + repmat(1 * eye(Nchan),[1,1,Nfft]);
Rcor(:,:,:,2) = mean(XX(:,:,:,[1:10,Nframe-10:Nframe]),4);
%Rcor(:,:,:,2) = bsxfun(@times, mean(mean(mean(XX, 4),1),2), repmat(eye(Nchan),[1,1,Nfft]));
RcorDet = zeros(Nbin, Nsource);
RcorInv = zeros(Nchan,Nchan,Nbin,Nsource);
for EMiter = 1:EMITERNUM
    for flp = 1:Nbin
        for slp = 1:Nsource
            RcorDet(flp,slp) = abs(det(Rcor(:,:,flp,slp)));
            condRcor = rcond(Rcor(:,:,flp,slp));
            if (condRcor < eps)
                return;
            end
            RcorInv(:,:,flp,slp) = inv(Rcor(:,:,flp,slp));
        end
    end
    
    %the LOG coefficient before exp term of gauss
    coeff = - bsxfun(@plus, Nchan * log(phiCor), permute(log(RcorDet),[1,3,2]));
    
    %ignore mu for simplity, yh * Rcor * y
    retr = abs((squeeze(sum(bsxfun(@times, sum(bsxfun(@times, conj(permute(ftbin, [1,4,2,3])), permute(RcorInv,[1,2,3,5,4])),1),  permute(ftbin,[4,1,2,3])), 2))));
    
    ex = -retr ./ phiCor;
    
    exc = coeff + ex;
    %for calculate cost
    %ct = exc;
    exc = bsxfun(@minus, exc, max(exc,[],3)); 
    softmask = bsxfun(@rdivide, exp(exc), sum(exp(exc), 3));

    %calculate cost to ensure converge
    %tt = ct .* softmask;
    %cost = sum(tt(:))   
    
    %update Rcor
%     Rcor =  bsxfun(@rdivide, squeeze(sum(bsxfun(@times, permute(softmask./((phiCor)),[4,5,1,2,3]),XX),4)),...
%                              permute(squeeze(sum(softmask,2)),[3,4,1,2]));
                         
    Rcor =  bsxfun(@rdivide, permute(sum(bsxfun(@times, permute(softmask./((phiCor)),[4,5,1,2,3]),XX),4), [1,2,3,5,4]),...
                             permute(squeeze(sum(softmask,2)),[3,4,1,2]));                         
    phiCor = retr / Nchan;
    %max(phiCor(:))
    %disp(['EMiter',num2str(EMiter)]);
end