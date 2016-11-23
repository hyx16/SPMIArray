% Copyright 2016 Tsinghua University SPMI Lab, Hongyu Xiang
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

function noise = NOISE(ftbin, Gcor, Ncor, refMic)
% estimate noises, not used yet
[Nchan, Nbin, ~] = size(ftbin);
u0 = zeros(Nchan,1);
u0(refMic) = 1;

Tcor = Gcor + Ncor;
TcorInv = zeros(size(Tcor));    
    
for flp=1:Nbin
    TcorInv(:,:,flp) = inv(Tcor(:,:,flp));
end
noise = squeeze(sum(bsxfun(@times, conj(sum(bsxfun(@times, TcorInv, permute(sum(bsxfun(@times, Ncor, u0.'),2),[2,1,3])),2)), permute(ftbin,[1,4,2,3])),1));