function [EEG] = GT_cleanEpochs(sub,exp, EEG)

EEG.chanlocs(129).labels = 'VEOG1';
EEG.chanlocs(130).labels = 'VEOG2';
EEG.chanlocs(131).labels = 'HEOG1';
EEG.chanlocs(132).labels = 'HEOG2';

% Remove trials marked for rejection during EB search 
if ~isempty(EEG.eyeBlink_rej)
    EEG.data(:,:,EEG.eyeBlink_rej) = []; %delete rejected trials; 
    EEG.epoch(EEG.eyeBlink_rej) = [];
%     EEG.event = EEG1.event; EEG.urevent = EEG1.urevent;
    EEG.trials = size(EEG.data,3);
    EEG.sumData(EEG.eyeBlink_rej,:) = [];
else 
    disp('Debugging')
end

%Change search epoch according to epoch times
if ~contains(EEG.filename, 'RL')
    EEG1 = pop_eegthresh(EEG, 1, 1:exp.nEEGchans,-150 , 150 ,-0.2 , 2, 1, 0); [EEG1,comrej, badlist]= pop_TBT(EEG1,EEG1.reject.rejthreshE,20,0.3,1);
elseif contains(EEG.filename, 'RL')
    EEG1 = pop_eegthresh(EEG, 1, 1:exp.nEEGchans,-150 , 150 ,-2 , 0.5, 1, 0); [EEG1,comrej, badlist]= pop_TBT(EEG1,EEG1.reject.rejthreshE,20,0.3,1);
end

% %For P1, for example, see correction on  Epoch 14 for ch 95,101
% %figure; plot(EEG.data(50,:,14)); hold on ; plot(EEG1.data(50,:,14))

if ~isempty(badlist.bTrial_num)
    EEG.data(:,:,badlist.bTrial_num) = []; %delete rejected trials; 
    EEG.epoch(badlist.bTrial_num) = [];
%     EEG.event = EEG1.event; EEG.urevent = EEG1.urevent;
    EEG.trials = size(EEG.data,3);
    EEG.sumData(badlist.bTrial_num,:) = [];
else 
    disp('Debugging')
end
    
%Rewrite EEG channels on original datafile; keep externals as they were.
EEG.data(1:128,:,:) = EEG1.data(1:128,:,:)
EEG.badList_TBT = badlist; 

end
