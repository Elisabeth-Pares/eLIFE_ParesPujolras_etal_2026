function [EEG] = GT_runICA(sub,exp,EEG, speedMethod)

%Options to speed up ICA computation: 
%-Downsample ICA
%-Run on data subset

switch speedMethod
    case 'downsample'
        EEG1 = EEG;
        [EEG1] = pop_resample(EEG1,100); %Downsampling just for ICA
    case 'subset'
        %Select subset of epochs to run ICA faster
        EEG1 = EEG;
        EEG1.data = EEG1.data(:,:,1:100);
        EEG1.epoch = EEG1.epoch(:,1:100);
    case 'none'
        EEG1 = EEG;
end
    
% EEG1 = EEG;
tic()
EEG1 = pop_runica(EEG1, 'icatype', 'runica', 'extended',0,'interrupt','on', 'chanind', exp.eegChans) % exp.icaChans);
toc()

EEG.icawinv = EEG1.icawinv;
EEG.icasphere = EEG1.icasphere;
EEG.icaweights = EEG1.icaweights;
EEG.icachansind = EEG1.icachansind;

end