function [allP, cppAll, tag] = GT_averageBygroup_tidy(id,exp, par, allP, cppAll)

sub = exp.sub_id{exp.thisExp}(id);

filename = ['csd_preproc_' par.epoch '_' exp.thisExpCode '_P' num2str(sub) '.set'];
EEG = pop_loadset([exp.filepath '\' filename]);

if strcmp(par.dataType, 'ERP')
    thisData = EEG.data;
    % Re-baseline trials for extra safety - trial-by-trial interpolation may
    % have changed some single trial, single-channel baselines.
    baseTime = [round(1.3*512) : round(1.499*512)];
    if strcmp(par.epoch, 'SL')
        baselines = nanmean(thisData(:,baseTime,:),2);
        thisData = thisData - baselines;
    end
    sr = 512; %Original sampling rate
    
else

    filename = [exp.tfapath, '/Power_DS_csd_preproc_' par.epoch '_' exp.thisExpCode '_P' num2str(sub) '..mat'];
    
    load([filename]);
    if strcmp(par.band, 'alpha');
        thisData = sumDat.alphaBand;
    elseif strcmp(par.band, 'beta')
        thisData = sumDat.betaBand;
    elseif strcmp(par.band, 'alphaBeta')
        thisData = sumDat.alphaBetaBand;
    end
    sr = 200; %TFA data were downsampled
    
end
%Check that sumData and epochs align
% assert(size(thisData, 3) == height(EEG.sumData))
if strcmp(par.epoch, 'RL')
    if size(thisData, 3) > height(EEG.sumData) %One extra response
        filename = ['csd_preproc_P1_' exp.thisExpCode '_P' num2str(sub) '.set'];
        EEG2 = pop_loadset([exp.filepath '\' filename]);
        for e = 1:height(EEG.sumData)
            if EEG2.epoch(e).eventurevent{1} > EEG.epoch(e).eventurevent{1}
                thisData(:,:,e) = [];
                EEG.epoch(e) = [];
                if size(thisData, 3) > height(EEG.sumData) %One extra response
                    e = e-1;
                    continue;
                else
                    break
                end
            end
        end
        
    end
    clear EEG2; clear e;
end
assert(size(thisData, 3) == height(EEG.sumData))

allRemoved = [];
indices = unique(EEG.sumData.trial);
rl_thisData = [];
if strcmp(par.epoch, 'P1')
for i = 1:size(thisData,3)
    %Align to response
    trialOnset = 1.5*sr;
    RT = EEG.metadata(EEG.metadata.trial == indices(i) & EEG.metadata.event == 1,:).RT*sr;
    trialDur = (EEG.metadata(EEG.metadata.trial == indices(i) & EEG.metadata.event == 4,:).evOffset - EEG.metadata(EEG.metadata.trial == indices(i) & EEG.metadata.event == 21,:).evOnset) *sr;
    thisResp = trialOnset + RT + trialDur;
    if (thisResp + 0.2*sr) <= 900
        rl_thisData(:,:,i) = thisData(:,round(thisResp - 1*sr):round(thisResp +0.2*sr),i);
    else
        rl_thisData(:,:,i) = repmat(NaN,size(thisData, 1), 241,1);
    end
end
allP.P1.allCond_rl(:,:,id) = nanmean(rl_thisData,3);
end
exp.resLabs = {'lowP1w', 'highP1w'};
if par.average == 1
    for group = 1:length(par.groupings)
        thisGroup= par.groupings{group};
        sortBy = GT_getConditions(thisGroup); %Get sorting parameters
        
        if sortBy.gap == 1; gaps = exp.gapLabs{exp.thisExp}; else, gaps = 1; end
        if sortBy.p1coh == 1; p1cohs = exp.p1Labs; else, p1cohs = 1; end
        if sortBy.p2coh == 1; p2cohs = exp.p2Labs; else, p2cohs = 1; end
        if sortBy.p1p2coh == 1; p1cohs = exp.p1Labs; p2cohs = exp.p2Labs; end
        if sortBy.dir == 1; dirs = exp.dirLabs; else, dirs = 1; end
        if sortBy.mov == 1; movs = exp.respKeys; else, movs = 1; end
        if sortBy.acc == 1; accs = exp.accuracies; else accs = 1; end
        
        
        for np = 1:2
            nLab = exp.npulseLabs{np}; nTag = 'np_';
            
            for g = 1:length(gaps) %gap duration
                thisGap = exp.gapDur{exp.thisExp}(g);
                
                for p1coh = 1:length(p1cohs) % p1 coherence
                    this_p1Coh = exp.cohLabs{p1coh};
                    
                    for p2coh = 1:length(p2cohs) % p2 coherence
                        this_p2Coh = exp.cohLabs{p2coh};
                        
                        for d = 1:length(dirs) % direction
                            thisDir = exp.dotDirs(d);
                            
                            for m = 1:length(movs)
                                thisMov = exp.respKeys(m);
                                
                                for acc = 1:length(accs)
                                    thisAcc = exp.accuracies(acc);
                                    
                                    
                                    %% SELECT TRIALS
                                    
                                    %All trials - ntrial vector with ones to use when nogrouping per condition is done
                                    allTrials = repmat(1,size(thisData,3),1);
                                    
                                    %1p/2p trials
                                    nSorting = EEG.sumData.npulses == np;
                                    
                                    %Gap sorting
                                    if sortBy.gap == 1, gapSorting = (EEG.sumData.gapdur == thisGap); gapTag = 'gap_'; gapLab = exp.gapLabs{exp.thisExp}{g};
                                    else gapSorting = allTrials; gapTag= ''; gapLab = '';
                                    end
                                    
                                    %P1 coh sorting
                                    if sortBy.p1coh == 1, p1Sorting = strcmp(EEG.sumData.p1hlcoh, this_p1Coh); p1Tag = 'p1coh_'; p1Lab = this_p1Coh;
                                    else p1Sorting = allTrials; p1Tag= ''; p1Lab = '';
                                    end
                                    
                                    %P2 coh sorting
                                    if sortBy.p2coh == 1, p2Sorting = strcmp(EEG.sumData.p2hlcoh, this_p2Coh); p2Tag = 'p2coh_';p2Lab = this_p2Coh;
                                    else p2Sorting = allTrials; p2Tag= ''; p2Lab = '';
                                    end
                                    
                                    %P1P2 coh sorting
                                    if sortBy.p1p2coh == 1, p1p2Sorting = (strcmp(EEG.sumData.p1hlcoh, this_p1Coh) & strcmp(EEG.sumData.p2hlcoh, this_p2Coh)); p1p2Tag = 'p1p2coh_'; p1p2Lab = [this_p1Coh '_' this_p2Coh];
                                    else p1p2Sorting = allTrials; p1p2Tag= ''; p1p2Lab = '';
                                    end
                                    
                                    %Direction sorting
                                    if sortBy.dir == 1, dirSorting = (EEG.sumData.pdir == thisDir); dirTag = 'dir_'; dirLab = exp.dirLabs{d};
                                    else dirSorting = allTrials; dirTag = ''; dirLab = '';
                                    end
                                    
                                    %Action sorting
                                    if sortBy.mov == 1, actSorting = (EEG.sumData.keyCode == thisMov); actTag = 'mov_'; actLab = exp.movLabs{m};
                                    else actSorting = allTrials; actTag = ''; actLab = '';
                                    end
                                    
                                    %Accuracy sorting
                                    if sortBy.acc == 1, accSorting = (EEG.sumData.accuracy == thisAcc); accTag = 'acc_'; accLab = exp.accLabs{acc};
                                    else accSorting = allTrials; accTag = ''; accLab = '';
                                    end
                                    
                                    %Currently selected trials
                                    idx = nSorting & gapSorting & p1Sorting & p2Sorting & p1p2Sorting & dirSorting & actSorting & accSorting;
                                    
                                    %Current sorting tag
                                    tag = [nTag, gapTag, p1Tag, p2Tag, p1p2Tag, dirTag, actTag, accTag];
                                    
                                    %Current condition
                                    thisTag = [nLab, gapLab, p1Lab, p2Lab, p1p2Lab, dirLab, actLab, accLab];
                                    
                                    %% AVERAGE
                                    % Get mean (w/o nan) of selected trials
                                    if ~isfield(allP, par.epoch); allP.(par.epoch) = struct(); end
                                    if ~isfield(allP.(par.epoch), tag); allP.(par.epoch).(tag) = []; end
                                    if ~isfield(allP.(par.epoch).(tag), thisTag); allP.(par.epoch).(tag).(thisTag) = []; end
                                    if ~isfield(allP.(par.epoch).(tag).(thisTag), 'sl'); allP.(par.epoch).(tag).(thisTag).sl = repmat(NaN,size(thisData,1),size(thisData,2)); end
                                    
                                    if sum(idx)>1
                                        allP.(par.epoch).(tag).(thisTag).sl(:,:,id) = nanmean(thisData(:,:,idx),3); %Average across selected trials
                                    else
                                        allP.(par.epoch).(tag).(thisTag).sl(:,:,id) = repmat(NaN,size(thisData,1),size(thisData,2)); %If only one trial, set to NaN
                                    end
                                    
                                    allP.(par.epoch).(tag).(thisTag).nTrials(:,id) = sum(idx); %Save number of trials per condition
                                    
                                                                     
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
end
