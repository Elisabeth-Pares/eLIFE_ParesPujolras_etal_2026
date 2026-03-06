function [EEG] = GT_addMetadata(sub,exp, EEG) 
% Add stimulus details to EEG dataset for later sorting 

%Load behavioural file
metadata = readtable([exp.behpath '/' exp.thisExpCode '_ExpResults_P' num2str(sub) '.csv']);
%Drop first row
metadata = metadata(metadata.ID ~= 0,:);
metadata = sortrows(metadata, 'trial');
     
%Get one row per trial  -take only 1st pulses
summaryData = metadata(metadata.event == 21,:);

%Deal with special subjects 
% Sub 3: Missing EEG data for first trials in block 6, 9 and 1
% Sub 4: Missing 1st trial in blocks 11 and 13
if isfield(EEG, 'missingEEG')
    summaryData(EEG.missingEEG,:) = []; %If any EEG trials are missing (recording started late), remove from summary data 
end
    
%Number of blocks in each experiment. All blocks had 70 trials. 
if exp.thisExp == 1
    nBlock = 10;
elseif exp.thisExp == 2
    nBlock = 14;
end  

%Add block data
%Initialise
summaryData.block = repmat(999,height(summaryData),1);
for b = 1:nBlock
    idx = (1:70)+70*(b-1);
    [summaryData(summaryData.trial <= idx(end) & summaryData.trial >= idx(1),:).block] = repmat(b,sum(summaryData.trial <= idx(end) & summaryData.trial >= idx(1)),1);
end

%Add high-low coherence data 

%Initialise fields
summaryData(:,'lowCoh') = deal(table(999)); %{'null'};
summaryData(:,'highCoh') = deal(table(999)); %{'null'};
summaryData(:,'p1hlcoh') = {'null'};
summaryData(:,'p2hlcoh') = {'null'};
summaryData(:,'avCoh') = {'null'};
summaryData(:,'bothCoh') = {'null'};

trial.p1lowCoh = []; trial.p1highCoh= []; trial.p2lowCoh = []; trial.p2highCoh = []; trial.lowCoh = []; trial.highCoh = []; trial.midCoh = [];
trial.low_low = []; trial.high_high = []; trial.low_high = []; trial.high_low = [];


for b = 1:nBlock
   idx = summaryData.block == b;
   lowCoh = min(summaryData(idx,:).p1coh); 
   highCoh = max(summaryData(idx,:).p1coh); 
   
   summaryData{idx,'lowCoh'} =  lowCoh; 
   summaryData{idx, 'highCoh'}= highCoh; %,70,1);
   
   %Pulse 1
   trial.p1lowCoh=  [trial.p1lowCoh; summaryData(idx,:).p1coh == lowCoh];
   trial.p1highCoh = [trial.p1highCoh; summaryData(idx,:).p1coh == highCoh];
   %Pulse 2
   trial.p2lowCoh = [trial.p2lowCoh; summaryData(idx,:).p2coh == lowCoh];
   trial.p2highCoh = [trial.p2highCoh; summaryData(idx,:).p2coh == highCoh];
    
   clear lowCoh; clear highCoh;
end

trial.lowCoh = [trial.p1lowCoh & trial.p2lowCoh];
trial.highCoh = [trial.p1highCoh & trial.p2highCoh];
trial.midCoh = [(trial.p1highCoh & trial.p2lowCoh) | (trial.p1lowCoh & trial.p2highCoh)];

trial.low_low = [trial.p1lowCoh & trial.p2lowCoh];
trial.high_high = [trial.p1highCoh & trial.p2highCoh];
trial.low_high = [(trial.p1lowCoh & trial.p2highCoh)];
trial.high_low = [(trial.p1highCoh & trial.p2lowCoh)];

   
summaryData(find(trial.p1lowCoh),'p1hlcoh') = deal({'low'});
summaryData(find(trial.p1highCoh),'p1hlcoh') = deal({'high'});
summaryData(find(trial.p2lowCoh),'p2hlcoh') = deal({'low'});
summaryData(find(trial.p2highCoh),'p2hlcoh') = deal({'high'});

summaryData(find(trial.lowCoh), 'avCoh') = deal({'lowCoh'});
summaryData(find(trial.highCoh), 'avCoh') = deal({'highCoh'});
summaryData(find(trial.midCoh), 'avCoh') = deal({'midCoh'});

summaryData(find(trial.low_low), 'bothCoh') = deal({'low_low'});
summaryData(find(trial.high_high), 'bothCoh') = deal({'high_high'});
summaryData(find(trial.low_high), 'bothCoh') = deal({'low_high'});
summaryData(find(trial.high_low), 'bothCoh') = deal({'high_low'});

EEG.metadata = metadata;
EEG.sumData = summaryData;

end
