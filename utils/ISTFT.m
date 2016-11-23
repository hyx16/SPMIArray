% Copyright 2016 Tsinghua University SPMI Lab, Hongyu Xiang
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

function y = ISTFT(FTbin, Lwindow, overlap)
%FTbin nbin * nframe * nchan
FTbin = cat(1,FTbin, conj(FTbin(end-1:-1:2,:,:)));
[Nfft, Nframe, Nchan] = size(FTbin);
Lgap = Lwindow - Lwindow * overlap;
y = zeros(Nframe * Lgap + Lwindow * overlap, Nchan);
for clp = 1:Nchan
    for tlp = 1:Nframe
        y((tlp - 1)*Lgap+1:(tlp-1)*Lgap+Lwindow,clp) = y((tlp - 1)*Lgap+1:(tlp-1)*Lgap+Lwindow,clp)+ sqrt(Nfft) * ifft(FTbin(:,tlp,clp), Lwindow);
    end
end