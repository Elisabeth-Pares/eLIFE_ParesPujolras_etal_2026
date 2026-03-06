function [exp, plt, pars] = GT_setup(exp, expN)

%% EXPERIMENT & PARTICIPANT INFO 
exp.code            = {'EP1', 'EP2'};
exp.name            = {'Experiment 1', 'Experiment 2'};
exp.numbers         = [1,2];
% exp.nsub              = length(exp.sub_id);

exp.thisExpName= exp.name{exp.thisExp};
exp.thisExpCode= exp.code{exp.thisExp};

%% PATH INFO 
addpath(genpath(['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/']));
addpath(genpath('C:/Users/elisa/Desktop/utils'));
cd('C:/Users/elisa/Desktop/GapsTask/');

exp.behpath           = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/behavioural data'];
exp.convpath          = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName  '/Data/converted data'];
exp.database          = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/raw data']
exp.filepath          = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/preprocessed data']%['C:/Users/elisa/Desktop/GapsTask/Data/preprocessed data'];
exp.tfapath           = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/tfa data'];
exp.plotpath          = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/plot data'];
exp.figurepath        = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/figures'];
exp.figPath           = ['C:/Users/elisa/Desktop/GapsTask/FinalFigures/'];
exp.mepath            = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/me data'];
exp.regpath           = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/reg data'];
exp.simpath           = ['C:/Users/elisa/Desktop/GapsTask/' exp.thisExpName '/Data/sim data'];

addpath('C:/Users/elisa/Desktop/GapsTask/');
addpath('C:/Users/elisa/Desktop/Data back-up/GitHub/GitHub/utils');
%% SETUP INFO 
exp.fs                = 512; %EEG sampling rate 

% External channel info 
exp.chan.veog1 = 129; 
exp.chan.veog2 = 130; 
exp.chan.heog1 = 131; 
exp.chan.heog2 = 132; 
exp.chan.emgr1 = 133; 
exp.chan.emgr2 = 134; 
exp.chan.emgl1 = 135; 
exp.chan.emgl2 = 136; 

exp.nEEGchans = 128;
exp.nchans = 137; % Last channel is empty; 

exp.eegChans = [1:128];
exp.eogChans = [129,130,131,132];

exp.srate    = 512; %sampling rate 

%% Experimental parameters 
exp.pulseDur = 0.2; % in s

%Gap duration is the only parameter varying in the two tasks
exp.gapDur{1} = [0,0.5,1]; %in seconds
exp.gapDur{2} = [0,0.120, 0.360, 1.080]; %in seconds

exp.gapLabs{1} = {'gap0', 'gap05', 'gap1'};
exp.gapLabs{2} = {'gap0', 'gap120', 'gap360', 'gap1080'};

exp.respKeys = [67,77]; %Left (c), right (m)
exp.accuracies = [0,1]; %Wrong, Right
exp.dotDirs = [0,-180]; %0: - Right, -180: left 

exp.cohLabs = {'high', 'low'};
exp.p1Labs = {'p1coh_high', 'p1coh_low'};
exp.p2Labs = {'p2coh_high', 'p2coh_low'};
exp.p1cohLabs = {'p1high', 'p1low'};
exp.p2cohLabs = {'p2high', 'p2low'};
exp.dirLabs = {'dotRight', 'dotLeft'};
exp.movLabs = {'Left', 'Right'};
exp.accLabs = {'Incorrect', 'Correct'};
exp.npulseLabs = {'onePulse', 'twoPulse'};
exp.cppLabs = {'Early', 'Late'};
exp.alphaLabs = {'Low', 'High'}; %PreStimulus power 
exp.repLabs = {'Alt', 'Rep'}; %Alternation/Repetition response

exp.movChans = [115,54]; %Central electrodes over left & right motor cortex 

exp.movCluster = [108,107,115,116,124,123,114;63,64,54,55,50,49,53]; %Left, right  

exp.parietalFocus = [4,32,31,30,21,17,18,5,19,20];

exp.blockLength = 70; 
%% Processing parameters
%% FILTERING 
exp.filter.lowerbound = 0.1; % When 0.1, files are labelled f01; when 0.01, files are just labelled f.
exp.filter.upperbound = 40;
exp.filter.order = 8;

if exp.filter.lowerbound == 0.01
    exp.filterLab = 'f';
elseif exp.filter.lowerbound == 0.1
    exp.filterLab = 'f01';
end

% Epochs 
exp.respEpoch= [-2, 1];
exp.epochs = {'P1', 'RL'};% Trigger parameters

% Analysis parameters & labels 
exp.filterLabs = {'_f', '_f01'}; %Filter lab = _f if filter .01, _f01 if .1
exp.icaLabs = {'', 'ab_cICAtrial_eyeB_rej_ri'; '', 'ab_cICAtrial_ri'}; %Rows: SL, RL; Columns: no ica, ica labels. 
exp.nanLabs = {'', 'NaN_'}; 
exp.csdLabs = {'', 'csd_prePulse'};
exp.tfaLabs = {'', 'Power'};
exp.dataLabs = {'ERP', 'TFA'};

exp.pulseLabs = {'p1', 'p2'};

% EEG parameters
exp.trigg.p1 = 5; exp.stimEpoch = [-1.5, 3]; % 3s is the max. duration of a trial, given max durations of LI & LO dots (0.75), two 0.2 pulses and a max gap of 1s.
exp.trigg.response = 13; 

exp.pulseTrigger = [5,8]; %p1, p2
% ICA channels
exp.icaChans = [1:132];

% Plotting parameters 
load chanlocsBiosemi_128_ext8 

exp.ImportantChans = [4 19 85 54 115]; % electrodes around CPz (where CPP is), Fz (CNV) and left and right motor cortex (see 'cap_128_layout_medium.jpg')
% Also mark front-most channels:
exp.FrontChans = [71 72 80 81 93 94 103 70 73 79 82 92 95 102]; % at the very front edge SD can be high not necessarily because channels are bad but because person is blinking /moving eyes. We don't want to interpolate GOOD electrodes that later will HELP us pick up artifacts that we want to throw out.
exp.rerefchan = [23]; % pick a channel to re-reference the data to and get a second picture of variance across channels - important just because sometimes the reference used to load the EEG might itself have been bad during the recording...

exp.cppChannels = [4,32,31,30,21,17,18,5,19,20];

exp.cppChans = {'A1','A2', 'A3','A4', 'A5', 'A6','A18', 'A19', 'A20', 'A31', 'A32', 'B1', 'B2', 'B3', 'C1', 'D1', 'D15', 'D16'};
exp.cppCentral = {'A1','A2', 'A3', 'B1', 'B2','C1', 'D1', 'D15', 'D16'};
exp.cppParietal= {'A4', 'A5', 'A6','A18', 'A19', 'A20', 'A31', 'A32', 'B3'};
exp.cnvChans = {'C28', 'C18', 'C15', 'C27', 'C19', 'C14','C26','C20','C13'}; %'C25', 'C21', 'C12','C24', 'C22', 'C11'};

%% ARTEFACT REJECTION 
exp.lowthresh = -150; 
exp.upthresh = 150;

% for el = 1:length(exp.cppChans)
%     thisCh = exp.cppChans(el)
%     exp.cppIdx(el) = find(strcmp({chanlocs.labels}, thisCh))
% end
% 
% for el = 1:length(exp.cppCentral)
%     thisCh = exp.cppCentral(el)
%     exp.cppCentralIdx(el) = find(strcmp({chanlocs.labels}, thisCh))
% end
% 
% for el = 1:length(exp.cppParietal)
%     thisCh = exp.cppParietal(el)
%     exp.cppParietalIdx(el) = find(strcmp({chanlocs.labels}, thisCh))
% end
% 
% for el = 1:length(exp.cnvChans)
%     thisCh = exp.cnvChans(el)
%     exp.cnvIdx(el) = find(strcmp({chanlocs.labels}, thisCh))
% end


%% TFA 
exp.tfa.min_freq =  1;     % Minimum frequency
exp.tfa.max_freq = 40;     % Maximum frequency
exp.tfa.num_frex = exp.tfa.max_freq - exp.tfa.min_freq+1;     % Specifies how many frequencies to compute. For a 1Hz resolution, do (max_freq-min_freq).
exp.tfa.frex = linspace(exp.tfa.min_freq,exp.tfa.max_freq,exp.tfa.num_frex); % Define vector with centre frequencies to compute.

exp.tfa.alphaBand =  [8:12]; %Alpha band
exp.tfa.betaBand  =  [13:30]; %Beta band
exp.tfa.alphaBetaBands = [8:30]; %Alpha & beta bands

exp.freqOfInterest = 'betaBand';

%% SIMULATION PARAMETERS 
pars.sr = 100;      %Simulation sampling rate
pars.nSim= 1000;    %Number of simulations
pars.maxIter =100;  %Max iterations for fitting 
pars.nRepMat=100;   
pars.stepAfterbndT = 0.1*pars.sr; %%When is activity reset back to 0 after hitting bound? 0.1s*sr. Originally: 50 (for 500 sr)
pars.accType = 'erf'; %How to compute accuracy 

%% PLOTTING PARAMETERS 
plt.smooth = 10; 
plt.alpha = 0.05; %transparency for SE shades

% Dimensions are always for SL, RL epochs 
plt.xlims = [1.3, 4.5; 0.5, 2.5] * exp.srate; % Epoch duration to plot, in s, and transformed to samples

plt.zero = [0.2;1.5];
plt.xline = [0.2; 1.5] * exp.srate; % to mark 0

plt.xticks = [0.2,1.2,2.2,3.2; 0.5, 1, 1.5, 2]*exp.srate;  
plt.xlabs = [0,1,2,3; -1,-0.5,0,0.5];
plt.xlabel = {'Time (s)'};

plt.ERP.timeLim = [1:512*2.2]; %2.2secs
plt.ERP.chan = [19];
plt.ERP.ylabel = {'EEG amplitude (\muV/m^{2})'}; 
plt.ERP.ylim = [-10 30];
plt.ERP.colorlimits_csd = [-13 13]; %For csd
plt.ERP.colorlimits = [-1.5 1.5]; %For nocsd

plt.TFA.chan = [86];
plt.TFA.ylabel = {'Beta power (\muV/m^{2})'}
plt.TFA.ylim = [8 16];
plt.TFA.colorlimits_csd = [-5 10]; %For csd
plt.TFA.colorlimits = [-0.2 0.2]; %For nocsd

plt.lineStyle = {'-', '--'}; %For left/right

plt.cohColors = [3]; %Coherence indicated in Blue shades; light --> dark (low, mid, high)
plt.gapColors = [5]; %Gap indicated in red shades; light --> dark (long, mid, zero gap)

plt.order = [1,3,5]; %Take 1st, 3rd and 5th colors for each condition 

plt.colors = {'k', 'b', 'r', 'g', 'm'}; 

plt.linetypes = {':', '-'}; %Linetypes for error & correct trials 

% COLORS
% 4 rows with GREEN, BLUE, YELLOW, RED palettes (5 colours each) - consistent with Shushruth's color scheme 
plt.palettes = {[0 60 48]/255, [1 102 94]/255, [53 151 143]/255, [128 205 193]/255, [199 234 229]/255;... % green
    [1,31,75]/255, [3,57,108]/255, [0,91,150]/255, [100,151,177]/255, [179,205,224]/255;... % blue % RF = 1 - CHECK!
    [84 48 5]/255, [140 81 10]/255, [191 129 45]/255, [223 194 125]/255, [246 232 195]/255; ... %yellow % RF = 2 - CHECK!
    [142,1,82]/255, [197,27,125]/255, [222,119,174]/255, [241,182,218]/255, [253,224,239]/255}; % red

%Generated with https://mdigi.tools/color-shades/
plt.satPalettes = {[13,70,7]/255, [22,116,12]/255, [31,162,16]/255, [39,209,21]/255;... % green -DONE
    [7,39,70]/255, [16,90,162]/255, [46,141,234]/255, [139, 192, 243]/255;... % blue  - DONE!
    [253 184 51]/255,[140 101 25]/255,[253 196 63]/255,[255 218 61]/255; ... %yellow  - CHECK!
    [70,7,7]/255, [162,16,16]/255, [234,46,46]/255, [243,139,139]/255;...}; % red- DONE!
    [70,7,63]/255, [162,16,148]/255, [234,46,216]/255, [243,139,233]/255;...%purple
    [70,35,0]/255, [162,81,0]/255, [234,117,0]/255, [243,121,0]/255;... %ORANGE - DONE
    [70,70,0]/256, [162,162,0]/255, [234,234,0]/255, [243,207,0]/255};%}; % GREENISH DONE!

%Generated using "shades" in https://mycolor.space/
plt.newPalettes = {    {'#09325A','#395682', '#637DAC', '#8DA7D8','#B9D2FF'},... %Blue
    {'#09115A', '#413788', '#7161B8','#A28EEA','#D4BEFF'},...%Purple
    {'#5A0909','#88352D','#B85E53','#E9897B','#FFB6A7'},... %Red
    {'#5A3909','#815B2B','#A97F4D','#D3A571','#FECD97'},... %Orange
    {'#115A09','#397C2C','#5D9F4D','#81C46F','#A6EA92'}} %Green

plt.Rcolors = {'#ffa600', '#ff6361', '#bc5090', '#58508d'};

plt.gapShades{1} = [plt.zero(1) plt.zero(1) plt.zero(1)+exp.pulseDur plt.zero(1)+exp.pulseDur; ... %First pulse, common across gaps
    plt.zero(1)+exp.pulseDur plt.zero(1)+exp.pulseDur plt.zero(1)+exp.pulseDur*2 plt.zero(1)+exp.pulseDur*2; ... %0-gap second pulse
    plt.zero(1)+exp.pulseDur+exp.gapDur{1}(2) plt.zero(1)+exp.pulseDur+exp.gapDur{1}(2) plt.zero(1)+exp.gapDur{1}(2)+exp.pulseDur*2 plt.zero(1)+exp.gapDur{1}(2)+exp.pulseDur*2; ... %0.5 gap second pulse
    plt.zero(1)+exp.pulseDur+exp.gapDur{1}(3) plt.zero(1)+exp.pulseDur+exp.gapDur{1}(3) plt.zero(1)+exp.gapDur{1}(3)+exp.pulseDur*2 plt.zero(1)+exp.gapDur{1}(3)+exp.pulseDur*2;...  %1 gap second pulse
    ]* exp.srate;  %0.5-gap second pulse *exp.srate;

plt.gapShades{2} = [plt.zero(1) plt.zero(1) plt.zero(1)+exp.pulseDur plt.zero(1)+exp.pulseDur; ... %First pulse, common across gaps
    plt.zero(1)+exp.pulseDur plt.zero(1)+exp.pulseDur plt.zero(1)+exp.pulseDur*2 plt.zero(1)+exp.pulseDur*2; ... %0-gap second pulse
    plt.zero(1)+exp.pulseDur+exp.gapDur{2}(2) plt.zero(1)+exp.pulseDur+exp.gapDur{2}(2) plt.zero(1)+exp.gapDur{2}(2)+exp.pulseDur*2 plt.zero(1)+exp.gapDur{2}(2)+exp.pulseDur*2; ... %0.12 gap second pulse
    plt.zero(1)+exp.pulseDur+exp.gapDur{2}(3) plt.zero(1)+exp.pulseDur+exp.gapDur{2}(3) plt.zero(1)+exp.gapDur{2}(3)+exp.pulseDur*2 plt.zero(1)+exp.gapDur{2}(3)+exp.pulseDur*2;...  %0.36 gap second pulse
    plt.zero(1)+exp.pulseDur+exp.gapDur{2}(4) plt.zero(1)+exp.pulseDur+exp.gapDur{2}(4) plt.zero(1)+exp.gapDur{2}(4)+exp.pulseDur*2 plt.zero(1)+exp.gapDur{2}(4)+exp.pulseDur*2;...  %0.36 gap second pulse
    ]* exp.srate;  %0.5-gap second pulse *exp.srate;

%TOPOPLOTS
plt.topoTimes = [250 500,1000,1500;-1500 -1000 -500 0];


set(gca, 'fontsize', 20)

exp.cppChans = [4,32,31,30,21,17,18,5,19,20];

