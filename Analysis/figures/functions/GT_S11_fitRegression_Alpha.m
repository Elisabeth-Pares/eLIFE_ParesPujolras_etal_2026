function [GA, allP] = GT_fitRegression(exp, par,id, GA, allP)

%Is weight of pulses on choice modulated by alpha power?
%Choice ~ p1coh + p2coh + p1 x alpha + p2 x alpha

sub = exp.sub_id{exp.thisExp}(id);
exp.parietalFocus = [4,32,31,30,21,17,18,5,...
    19,20];
%% Load data
%     Load alpha power
if exp.thisExp == 1
    filename = [exp.tfapath, '/Power_csd_preproc_P1_' exp.thisExpCode '_P' num2str(sub) '..mat'];
else
    filename = [exp.tfapath, '/Power_DS_csd_preproc_P1_' exp.thisExpCode '_P' num2str(sub) '..mat'];
end
load([filename]);

occCluster = [9,10,17,16,21,30,29,38,39,22]; %Occipital channels for alpha analysis

occAlpha = nanmean(sumDat.alphaBand([occCluster],:,:),1); %Time-resolved occipital alpha
baseline = nanmean(nanmean(sumDat.alphaBand([occCluster],1.3*200:1.5*200,:),1),2);
occAlpha = 10.*log10(occAlpha./baseline);

tonicAlpha = nanmean(occAlpha(:,1.3*200:1.5*200,:),2); %Tonic alpha: preStimulus alpha
ph_occAlpha = occAlpha - tonicAlpha; %Phasic alpha: alpha power - preStimulus alpha

%Load EEG data & behaviour
filename = ['csd_preproc_P1_' exp.thisExpCode '_P' num2str(sub) '.set'];
EEG = pop_loadset([exp.filepath '\' filename]);

allP.sumData{id} = EEG.sumData;
allP.occAlpha{id} = occAlpha;

if strcmp(par.modelName, 'PxCPPxAlpha') || strcmp(par.modelName, 'PxCPPxphAlpha')
    %Load and downsample CPP data
    % Which channels to use?
    if par.optiChan == 1 %Get chan with max. CPP amplitude for each participant
        theseChans = par.bestChan;
    else
        theseChans = exp.parietalFocus;
    end
    
    % Load EEG data & behaviour
    [EEG] = pop_resample(EEG,200); %Downsample data for speed and to make it match TFA data

    [bfilt,afilt] = butter(3, 6*2/EEG.srate, 'low');   % lo-pass
    for t = 1:size(EEG.data,3)
        for e = 1:size(EEG.data,1)
            EEG.data(e,:,t) = filtfilt(bfilt,afilt, double(EEG.data(e,:,t)));
        end
    end
    
    CPP = nanmean(EEG.data(theseChans ,:,:),1);
    allP.sumData{id} = EEG.sumData;
    allP.CPP{id} = CPP;
    zCPP = CPP;
    allP.zCPP{id} = zCPP;
end


if strcmp(par.modelName, 'PxCPP')
    timeVector = -1.5:1/512:3;
    timeRange = find(timeVector>-0.2 & timeVector < 2.2);
elseif strcmp(par.modelName, 'PxdeltaAlpha') || strcmp(par.modelName, 'Pxtonic_PxdeltaAlpha')
    timeRange = 1; %Only one value. DeltaAlpha calculated as alpha(P2+0.5) - alpha(P2-0.2)
else
    timeVector = -1.5:1/200:3;
    timeRange = find(timeVector>-0.7 & timeVector < 2.2);
end


%Run for each gap separately
for g = 1:length(exp.gapDur{exp.thisExp})
    
    thisGap = exp.gapDur{exp.thisExp}(g);
    idx = EEG.sumData.gapdur == thisGap & EEG.sumData.npulses == 2;
    gapLab = exp.gapLabs{exp.thisExp}{g};
    
    c = 0;
    %Time resolved
    for times = 1:length(timeRange)
        c = c+1;
        t = timeRange(times);
        
        
        if strcmp(par.modelName, 'PxAlpha')
            X = [nanzscore(EEG.sumData.p1coh(idx)).*sign(EEG.sumData.pdir(idx)+1),... %p1coh
                nanzscore(EEG.sumData.p2coh(idx)).*sign(EEG.sumData.pdir(idx)+1),...  %p2coh
                nanzscore(EEG.sumData.p1coh(idx)).*sign(EEG.sumData.pdir(idx)+1).*-nanzscore(squeeze(occAlpha(:,t,idx))),... %p1coh x alpha
                nanzscore(EEG.sumData.p2coh(idx)).*sign(EEG.sumData.pdir(idx)+1).*-nanzscore(squeeze(occAlpha(:,t,idx)))];   %p2coh x alpha 
        end
        
        responses = EEG.sumData.keyCode(idx);
        responses(responses == 67) = 0;
        responses(responses == 77) = 1;

        Y = responses;
        X = nanzscore(X);

        [B,dev,stats] = glmfit(X,Y,'binomial');
         
        if strcmp(par.modelName, 'PxAlpha')
            GA.(gapLab).p1coh(id,:,c) = B(2);          
            GA.(gapLab).p2coh(id,:,c) =  B(3);       
            GA.(gapLab).p1cohxAlpha(id,:,c) = B(4);         
            GA.(gapLab).p2cohxAlpha(id,:,c) = B(5); 
        end
        
        GA.time = timeVector(timeRange);
    end
    
end

end