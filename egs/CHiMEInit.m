% initialize path and directories for CHiME4 enhancent

% Copyright 2016 Tsinghua University SPMI Lab, Hongyu Xiang
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)

chimeRoot = 'G:/CHIME/CHiME3/';
enhRoot = 'G:/CHIME/CHiME3/data/BSS_ENH_test/';
workRoot = 'G:/CHIME/CHiME3/TEST/BSS_OPI/';
addpath utils
addpath enhan

Path.isolated = [chimeRoot,'data/audio/16kHz/isolated/'];
Path.enhanced = [chimeRoot,'data/audio/16kHz/enhanced/']; 
Path.embedded = [chimeRoot,'data/audio/16kHz/embedded/']; 
Path.backgrounds = [chimeRoot,'data/audio/16kHz/backgrounds/'];
Path.annotations = [chimeRoot,'data/annotations/']; % path to JSON annotations

Path.enhBss = [enhRoot,'enhBss/'];
Path.enhBssMvdrRtf = [enhRoot, 'enhBssMvdrRtf/'];
Path.enhBssMvdrEg = [enhRoot, 'enhBssMvdrEg/'];
Path.enhBssPmwf = [enhRoot, 'enhBssPmwf/'];
Path.enhGsc = [enhRoot, 'enhGsc/'];

Path.enhBssNoise = [enhRoot,'enhBssNoise/'];
Path.enhBssMvdrRtfNoise = [enhRoot, 'enhBssMvdrRtfNoise/'];
Path.enhBssMvdrEgNoise = [enhRoot, 'enhBssMvdrEgNoise/'];
Path.enhBssPmwfNoise = [enhRoot, 'enhBssPmwfNoise/'];

enhDirs = {'enhBss','enhBssMvdrRtf','enhBssMvdrEg','enhBssPmwf','enhGsc','enhBssNoise','enhBssMvdrRtfNoise','enhBssMvdrEgNoise','enhBssPmwfNoise'};
sets={'et05','dt05','tr05'};
envirs = {'bus','caf','ped','str'};
modes={'real','simu'};
for hlp = 1:length(enhDirs)
    for slp = 1:length(sets)
        for elp = 1:length(envirs)
            for mlp = 1:length(modes)
                tDir = [enhRoot enhDirs{hlp} '/' sets{slp} '_' envirs{elp} '_' modes{mlp} '/'];
                if ~exist(tDir,'dir')
                    mkdir(tDir);
                end
            end
        end
    end
end