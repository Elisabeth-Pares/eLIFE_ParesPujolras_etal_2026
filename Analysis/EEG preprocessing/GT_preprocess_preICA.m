function [] = GT_preprocess_preICA(sub,exp)
% ===============================================================================
% GapsTask EEG preprocessing pipeline up to ICA 
% INPUT: fc_EPX_PX. raw data, whole experiment, concatenated blocks. 
% OUTPUT: ICA_fc_EPX_PX
%
% Preprocessing steps include: 
% 1. Filtering (Butterworth bandpass filter)
% 2. Epoching (SL, RL)
% 3. Bad channel detection & interpolation - requires manual confirmation 
% 4. Re-referencing (CAR)
% 5. Artefact flagging - mark epochs containing eyeblinks during evidence pulses 
% 6. ICA for eye movement artefact cleaning
% ===============================================================================

% Load raw converted data
EEG = pop_loadset( [exp.convpath '\fc_' exp.thisExpCode '_P' num2str(sub) '.set'] ); %fc_

% Make sure event field is numeric
% if isa(EEG.event(2).type, 'char')
%     [EEG.event.type] = str2double({EEG.event.type})';
% end
       
% Filter - band-pass filter 
% Butterworth 8th order bandpass filter - relies on ERP lab EEGlab plugin
EEG  = pop_basicfilter( EEG,  1:exp.nEEGchans , 'Cutoff', [exp.filter.lowerbound  exp.filter.upperbound], 'Design', 'butter', 'Filter', 'bandpass', 'Boundary', [], 'order', exp.filter.order);

EEG.data(137,:) = []; %Remove empty signal channel

load('chanlocsBioSemi_128_ext8.mat') %Load channel locations
EEG.chanlocs = chanlocs(1:136); %Remove channel locations for eye & emg channels
EEG.nchan = 136;

% Exp1, P9: trials 319 to 352 had no signal - remove data.
if exp.thisExp == 1 & sub == 9 
    allTrials = find([EEG.event.type] == 5);
    removeEvents = [allTrials(319) : 2263]; 
    [EEG.event(removeEvents).type] = deal(NaN);
    EEG.sumData(319:352,:) = []; %Align to beh data
end

% Epoch data (SL, RL)
EEG1 = pop_epoch( EEG, {exp.trigg.p1}, exp.stimEpoch, 'epochinfo', 'yes'); %Aligned to P1 onset
EEG2 = pop_epoch(EEG, {exp.trigg.response}, exp.respEpoch, 'epochinfo', 'yes'); %Aligned to response

% Check bad channels manually and interpolate 
[EEG1, EEG2] = GT_manualChanCheck(sub,exp,EEG1,EEG2); %Interpolate channels with big variability (>50 SD)

% Re-reference
EEG1 = pop_reref(EEG1, [1:128] ,'exclude',[129:136] ,'keepref','on');
EEG2 = pop_reref(EEG2, [1:128] ,'exclude',[129:136] ,'keepref','on');
EEG1 = eeg_checkset( EEG1 ); EEG2 = eeg_checkset( EEG2 );
    
% Mark trials with eyeblinks during gaps to reject later
[EEG1] = GT_removeEyeblinkEpochs(sub,exp,EEG1); %Reject epochs with eyeblinks during pulses
EEG2.eyeBlink_rej = [EEG1.eyeBlink_rej]; % Reject same trials from RL epochs

% Run ICA
if exp.thisExp == 1 %Estimate on RL data; ICA decomp was cleaner in this dataset.
    EEG2 = GT_runICA(sub,exp, EEG2, 'downsample'); %Downsample data for ICA decomposition only 
    %Upload ICA weights to SL data
    EEG1.icawinv = EEG2.icawinv;
    EEG1.icasphere = EEG2.icasphere;
    EEG1.icaweights = EEG2.icaweights;
    EEG1.icachansind = EEG2.icachansind;
elseif exp.thisExp == 2 %Estimate on SL data; ICA decomp was cleaner in this dataset.
    EEG1 = GT_runICA(sub,exp, EEG1, 'downsample'); %Downsample data for ICA decomposition only 
    %Upload ICA weights to RL data
    EEG2.icawinv = EEG1.icawinv;
    EEG2.icasphere = EEG1.icasphere;
    EEG2.icaweights = EEG1.icaweights;
    EEG2.icachansind = EEG1.icachansind;
end

%Save data
EEG1 = pop_saveset(EEG1,['ICA_P1_' exp.thisExpCode '_P' num2str(sub)], exp.filepath);
EEG2 = pop_saveset(EEG2,['ICA_RL_' exp.thisExpCode '_P' num2str(sub)], exp.filepath);

end


