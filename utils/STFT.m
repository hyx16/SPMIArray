% Copyright 2016 Tsinghua University SPMI Lab, Hongyu Xiang
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

function [ftbin,Nframe,Nbin,Lspeech,speechFrame] =  STFT(x, Lwindow, overlap, Nfft)
%ftbin : nchan * nbin * nframe 

%separate frame
[lgth,Nchan] = size(x);
Loverlap = Lwindow * overlap;
Lgap = Lwindow - Loverlap;
Nframe = ceil((lgth - Loverlap)/ Lgap);
Lspeech = Nframe * Lgap + Loverlap;
x = [x; zeros(Lspeech - lgth, Nchan)];
speechFrame = zeros(Lwindow, Nframe, Nchan);
%win = 0.5*(1-cos(2*pi*(0:Lwindow-1)'/(Lwindow-1)));
win = window(@hann, Lwindow);

for clp = 1:Nchan
    for flp = 1:Nframe
        speechFrame(:, flp, clp) = x((flp - 1)*Lgap+1:(flp-1)*Lgap+Lwindow,clp) .* win;
    end
end

%speechFrame : Nfft * Nframe * Nchan
%ftbin = zeros(Nchan, Nbin, Nframe);
Nbin = Nfft / 2 + 1;
ftbin = permute(fft(speechFrame, Nfft, 1), [3,1,2]) / sqrt(Nfft);
%ftbin = permute(fft(speechFrame, Nfft, 1), [3,1,2]);
ftbin = ftbin(:,1:Nbin,:);