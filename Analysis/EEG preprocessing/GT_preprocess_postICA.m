function [] = GT_preprocess_postICA(sub,exp)
% ===============================================================================
% GapsTask EEG preprocessing pipeline after ICA correction 
% INPUT: cICA_EPX_PX. preprocessed data with ICA weights
% OUTPUT: preproc_EPX_PX
%
% Preprocessing steps include: 
% 1. Baseline data 
% 2. Trial-by-trial interpolation (requires manual check)
% ===============================================================================

% Load ICA corrected data
EEG1 = pop_loadset( [exp.filepath '\cICA_P1_' exp.thisExpCode '_P' num2str(sub) '.set'] );
EEG2 = pop_loadset( [exp.filepath '\cICA_RL_' exp.thisExpCode '_P' num2str(sub) '.set'] );

% Baseline data
t.prePulse_baseline_idx = ([EEG1.times] >= -200 & [EEG1.times] <= 0);
t.prePulse_baseline = squeeze(mean(EEG1.data(:,t.prePulse_baseline_idx,:),2));

for e = 1:size(EEG1.data,3)
    EEG1.data(:,:,e) = EEG1.data(:,:,e) - t.prePulse_baseline(:,e); %Baseline data
    EEG2.data(:,:,e) = EEG2.data(:,:,e) - t.prePulse_baseline(:,e); %Baseline data
end

% Trial by trial artefact correction 
[EEG1] = GT_cleanEpochs(sub,exp, EEG1)
[EEG2] = GT_cleanEpochs(sub,exp, EEG2)

%Save data
EEG1 = pop_saveset(EEG1,['preproc_P1_' exp.thisExpCode '_P' num2str(sub)], exp.filepath);
EEG2 = pop_saveset(EEG2,['preproc_RL_' exp.thisExpCode '_P' num2str(sub)], exp.filepath);


end

