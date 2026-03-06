function [EEG] = GT_alignmentCheck(sub,exp, EEG) 

%Do number of pulses & pulse events match?
t = 1;
for trigger = [4,5,8,9,12]; %leadin, p1, p2, gap, leadout
    if ~isnan(trigger)
        count(t,1) = sum([EEG.event(:).type] == trigger);
        if count(t,1) == 0 
            trigger = num2str(trigger)
            count(t,1) = sum(cell2mat(cellfun(@(c)strcmp(c,trigger), {EEG.event(:).type},'UniformOutput', false)));
        end            
            
    end
    t = t+1;
end

count(:,2) = [999];
t = 1;
for ev = [1,21,22,3,4];%leadin, p1, p2, gap, leadout
    count(t,2) = sum(EEG.metadata.event == ev);
    t = t+1;
end

count(:,3) = count(:,1) - count(:,2);
if sum(count(:,3) == 0)
    EEG.triggMatch = 'True';
else
    EEG.triggMatch = 'False';
    warning('Triggers do not match!')
end

EEG.eventCount = count;

clear count; 
end



