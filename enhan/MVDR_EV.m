% Copyright 2016 Tsinghua University SPMI Lab, Hongyu Xiang
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

function Y = MVDR_EV(ftbin, Gcor, Ncor)
    [Nchan,Nbin,~] = size(ftbin);
    Df = zeros(Nchan, Nbin);
    NcorInv = zeros(Nchan,Nchan,Nbin);
    for flp = 1:Nbin
        NcorInv(:,:,flp) = inv(Ncor(:,:,flp));
        [vv,dd] = eig(Gcor(:,:,flp));
        [~,ddind] = max(diag(dd));
        Df(:,flp) = vv(:,ddind);
    end
    %tdt : Nchan * 1 * Nbin
    tdt = sum(bsxfun(@times,NcorInv,permute(Df,[3,1,2])),2);
    Y = squeeze(sum(bsxfun(@times,conj(bsxfun(@rdivide, tdt, sum(bsxfun(@times,tdt,conj(permute(Df,[1,3,2]))),1))),permute(ftbin,[1,4,2,3])),1));
end