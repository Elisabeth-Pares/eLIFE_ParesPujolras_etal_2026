function [GA, allP] = GT_fitBehRegression_CPP(exp, par,id, GA, allP)

%Is weight of pulses on choice modulated by CPP amplitude? 
%Choice ~ p1coh + p2coh + p1 x CPP + p2 x CPP


sub = exp.sub_id{exp.thisExp}(id);

% Which channels to use?
theseChans = exp.parietalFocus;


% Load EEG data & behaviour
filename = ['csd_preproc_P1_' exp.thisExpCode '_P' num2str(sub) '.set'];
EEG = pop_loadset([exp.filepath '\' filename]);

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
    
    %Baseline CPP at P2
    if par.baseline == 1
        p2Idx = round((1.5+thisGap)*512);
        baseline = (nanmean(zCPP(:,p2Idx-0.1*512:p2Idx,idx),2));
        zCPP(:,p2Idx:end,idx) = zCPP(:,p2Idx:end,idx) - baseline;
    end
    
    c = 0;
    %Time resolved
    for times = 1:length(timeRange)
        c = c+1;        
        t = timeRange(times);
        
        
        if strcmp(par.modelName, 'PxCPP') 
               X = [nanzscore(EEG.sumData.p1coh(idx)).*sign(EEG.sumData.pdir(idx)+1),... %p1coh
                    nanzscore(EEG.sumData.p2coh(idx)).*sign(EEG.sumData.pdir(idx)+1),... %p2coh
                    nanzscore(EEG.sumData.p1coh(idx)).*sign(EEG.sumData.pdir(idx)+1).*nanzscore(squeeze(zCPP(:,t,idx))),... %p1coh x cpp amp
                    nanzscore(EEG.sumData.p2coh(idx)).*sign(EEG.sumData.pdir(idx)+1).*nanzscore(squeeze(zCPP(:,t,idx)))];   %p2coh x cpp amp
 
        end
        responses = EEG.sumData.keyCode(idx);
        responses(responses == 67) = 0; 
        responses(responses == 77) = 1; 
       
        Y = responses; 
        X = nanzscore(X);
      
        [B,dev,stats] = glmfit(X,Y,'binomial');

        if strcmp(par.modelName, 'PxCPP')
            GA.(gapLab).p1coh(id,:,c) = B(2);         
            GA.(gapLab).p2coh(id,:,c) =  B(3); 
            GA.(gapLab).p1cohxCPP(id,:,c) = B(4);          
            GA.(gapLab).p2cohxCPP(id,:,c) = B(5); 
        end
        GA.time = timeVector(timeRange);
    end
        
end

end