% speech enhancement for CHiME4 6-ch data
% corresponding to "The THU-SPMI CHiME-4 system : Lightweight design with
% advanced multi-channel processing, feature enhancement, and language modeling"

% based on CHiME3 official code
% Jon Barker, Ricard Marxer, Emmanuel Vincent, and Shinji Watanabe, The
% third 'CHiME' Speech Separation and Recognition Challenge: Dataset,
% task and baselines, submitted to IEEE 2015 Automatic Speech Recognition
% and Understanding Workshop (ASRU), 2015.

% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

clear all;
CHiMEInit;

Nsource = 2;
Nchan = 6;
Lwindow = 256;
Nfft = Lwindow;
overlap = 0.75;
powThresh = -10;
cmin = 6400; % minimum context duration (400 ms)
cmax = 12800; % maximum context duration (800 ms)
EMITERNUM = 20;
GAMMA = 20;

sets = {'et05' 'dt05' 'tr05'};
modes = {'real' 'simu'};
micFailFid = fopen([workRoot 'micfail.txt'],'a');


setIndBegin = 1;
setIndEnd = 2;
modeIndBegin = [1 1 1];
modeIndEnd = [1 1 1];
uttIndBegin = [1    1;
               1    1;
               1    1];
    
%for setInd = setIndBegin:length(sets)
for setInd = setIndBegin:setIndEnd    
    set = sets{setInd};
    %for modeInd = modeIndBegin(setInd):length(modes)
    for modeInd = modeIndBegin(setInd):modeIndEnd(setInd)        
        mode = modes{modeInd};
        
        % Read annotations
        mat = json2mat([Path.annotations set '_' mode '.json']);
        realMat = json2mat([Path.annotations set '_real.json']);
        
        for uttInd = uttIndBegin(setInd,modeInd):length(mat)
            Nchan = 6;
            disp([set, ' ', mode, ' ', 'uttInd ',num2str(uttInd)]);
            sem = [set '_' lower(mat{uttInd}.environment) '_' mode];
            udir = [Path.isolated sem '/'];
            uname = [mat{uttInd}.speaker '_' mat{uttInd}.wsj_name '_' mat{uttInd}.environment];
            
            % Load WAV files
            xsize = size(wavread([udir uname '.CH1.wav']));
            x = zeros(xsize(1),Nchan);
            xc = zeros(xsize(1),Nchan);
            for clp = 1:Nchan,
                [x(:,clp),fs] = wavread([udir uname '.CH' int2str(clp) '.wav']);
                xc(:,clp) = x(:,clp) / norm(x(:,clp)) * norm(x(:,1));
            end
            
            % Check microphone failure
            CLwindow = 2048;
            Coverlap = 0.5;
            CNfft = CLwindow;
            [~,~,~,~,CspeechFrame] =  STFT(xc, CLwindow, Coverlap, CNfft);
            CspeechFrame = permute(CspeechFrame, [3,1,2]);
            NcorFFT = 2 * (CNfft-1);
            c1 = fft(repmat(permute(CspeechFrame, [1,4,2,3]),[1 Nchan 1 1]), NcorFFT, 3);
            c2 = fft(repmat(permute(CspeechFrame, [4,1,2,3]),[Nchan 1 1 1]), NcorFFT, 3);
            corr = real(ifft(c1 .* conj(c2),NcorFFT,3));
            corr = squeeze(max(corr, [], 3));
            corrjudge = sum(corr, 2) - sum(bsxfun(@times,corr,eye(Nchan)),2);
            corrjudge = squeeze(corrjudge);
            corrjudge_rate = bsxfun(@rdivide, corrjudge, median(corrjudge,1));
            fail = corrjudge_rate < 0.6;
            fail = (sum(fail,2) > 1);
            fail = fail';
            if all(fail)
                fail = (sum(corrjudge_rate,2) == max(sum(corrjudge_rate,2)));
                fail = ~fail;
                fail = fail';
            end
            if (any(fail))
                disp(['find one failure: ',set,' ',mode,' ','uttInd ',num2str(uttInd),' ',uname,' : ',num2str(fail)]);
                fprintf(micFailFid, [set,' ',mode,' ','uttInd ',num2str(uttInd),' ',uname,' : ',num2str(fail),'\n']);
            end
            
            % choose reference mic
            [~,refMic] = max(sum(corrjudge(~fail,:),2));          
            
            Nchan = sum(~fail);
            x = x(:,~fail);
                        
            [ftbin,Nframe,Nbin,Lspeech] =  STFT(x, Lwindow, overlap, Nfft);   
            % for GSC fixed beamformer
            targetY = squeeze(mean(ftbin,1));
            
            % correlation of X
            XX = bsxfun(@times, permute(ftbin,[1,4,2,3]), conj(permute(ftbin,[4,1,2,3])));
            
            Xcor = mean(XX, 4);
            % Ncor = mean(XX(:,:,:,[1:10,Nframe-10:Nframe]),4);
            % Load context (up to 5 s immediately preceding the utterance)
            % noise = read_context(nchan, c_ind, Path, mode, mat, real_mat, uttInd, cmax, cmin);            
            softmask = cGaussMask(ftbin,Nsource,XX,EMITERNUM);
            Ncor = bsxfun(@rdivide, mean(bsxfun(@times, XX, permute(softmask(:,:,2),[3,4,1,2])),4),permute(mean(softmask(:,:,2),2), [2,3,1]));                       
            Gcor = Xcor - Ncor;
            
            %for GSC
            yGscFt= GSC(ftbin, targetY);
            yGsc = ISTFT(yGscFt, Lwindow, overlap);
            yGsc = yGsc / max(abs(yGsc));
            wavwrite(yGsc, 16000, 16, [Path.enhGsc sem '/' uname '.wav']);                 
            
            
            %for bss
            yBssFt = bsxfun(@times, ftbin, permute(softmask(:,:,1),[3,1,2]));  
            yBssFt = squeeze(yBssFt(refMic,:,:));
            nBssFt = bsxfun(@times, ftbin, permute(softmask(:,:,2),[3,1,2]));  
            nBssFt = squeeze(nBssFt(refMic,:,:));
            yBss = ISTFT(yBssFt,Lwindow,overlap);
            nBss = ISTFT(nBssFt,Lwindow,overlap);
            
            SNR = sum(yBss.^2) / sum(nBss.^2);
            yBss = yBss / max(abs(yBss));
            nBss = nBss / max(abs(nBss));
            
            wavwrite(yBss, 16000, 16, [Path.enhBss sem '/' uname '.wav']);
            wavwrite(nBss, 16000, 16, [Path.enhBssNoise sem '/' uname '.wav']);
            
            %for bssPmwf
            yBssPmwfFt = PMWF(ftbin, Gcor, Ncor, refMic, min(GAMMA/SNR,50));
            nBssPmwfFt = NOISE(ftbin, Gcor, Ncor, refMic);
            yBssPmwf = ISTFT(yBssPmwfFt,Lwindow,overlap);
            nBssPmwf = ISTFT(nBssPmwfFt,Lwindow,overlap);
            yBssPmwf = yBssPmwf / max(abs(yBssPmwf));
            nBssPmwf = nBssPmwf / max(abs(nBssPmwf));
             
            wavwrite(yBssPmwf, 16000, 16, [Path.enhBssPmwf sem '/' uname '.wav']);
            wavwrite(nBssPmwf, 16000, 16, [Path.enhBssPmwfNoise sem '/' uname '.wav']);
            
            %for bssMvdrEg
            yBssMvdrEgFt = MVDR_EV(ftbin, Gcor, Ncor);
            yBssMvdrEg = ISTFT(yBssMvdrEgFt,Lwindow,overlap);
            yBssMvdrEg = yBssMvdrEg / max(abs(yBssMvdrEg));
            wavwrite(yBssMvdrEg, 16000, 16, [Path.enhBssMvdrEg sem '/' uname '.wav']);
            
            %for bssMvdrRtf
            yBssMvdrRtfFt = MVDR_RTF(ftbin, Gcor, Ncor, refMic);
            yBssMvdrRtf = ISTFT(yBssMvdrRtfFt,Lwindow,overlap);
            yBssMvdrRtf = yBssMvdrRtf / max(abs(yBssMvdrRtf));
            wavwrite(yBssMvdrRtf, 16000, 16, [Path.enhBssMvdrRtf sem '/' uname '.wav']);            
            
         end
     end
end
fclose(micFailFid);