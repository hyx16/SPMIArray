% Copyright 2016 Tsinghua University SPMI Lab, Hongyu Xiang
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

function Y = PMWF(ftbin, Gcor, Ncor, refMic, gamma)
% X_n(t,w) = H(w)S_n(t,w) + N_n(t,w)
% X_n(t,w) = G_n(t,x) + N_n(t,w)
% target is G_ref
% PMWF beamforming
%Y = zeros(size(X,1), size(X,2));
[Nchan, Nbin, ~] = size(ftbin);
u0 = zeros(Nchan,1);
u0(refMic) = 1;

Tcor = Gcor+gamma*Ncor;
TcorInv = zeros(size(Tcor));    
    
for flp=1:Nbin
    TcorInv(:,:,flp) = inv(Tcor(:,:,flp));
end
    
%Gcor : nchan * nchan * nbin 
%Ncor : nchan * nchan * nbin
%u0 : nhcan
%ftbin : nchan * nbin * nframe
Y = squeeze(sum(bsxfun(@times, conj(sum(bsxfun(@times, TcorInv, permute(sum(bsxfun(@times, Gcor, u0.'),2),[2,1,3])),2)), permute(ftbin,[1,4,2,3])),1));
