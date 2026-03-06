function [] = GT_applyCSD(sub,exp, epoch)

%Relies on CSD toolbox
%See https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/tutorial.html

if epoch == 1
    EEG = pop_loadset([exp.filepath '\preproc_P1_' exp.thisExpCode '_P' num2str(sub) '.set']);
elseif epoch == 2
    EEG = pop_loadset([exp.filepath '\preproc_RL_' exp.thisExpCode '_P' num2str(sub) '.set']);
end

%Load CSD transform
load('CSD_coords_biosemi128.mat');
       
figure; plot(mean(EEG.data(1:128,:,:),3)');
EEG = eeg_checkset( EEG );

%% Apply CSD to epoched, BASELINED data, but only to EEG channels
X = [];
tic();
disp(['Running CSD transform for P' num2str(sub) '...']);
for e = 1:length(EEG.epoch)
    X(:,:,e) = CSD(EEG.data(1:128,:,e), G,H);
end
toc();
%This takes about 3 minutes per participant!
EEG.data(1:128,:,:) = X;

EEG = eeg_checkset( EEG );
filename = ['csd_' EEG.filename];
EEG = pop_saveset( EEG, filename, exp.filepath);

end