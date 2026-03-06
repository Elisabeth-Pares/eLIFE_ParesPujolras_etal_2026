function [EEG] = GT_removeEyeblinkEpochs(sub,exp, EEG)

for epoch = exp.epochs(1) %Doing just for SL for now
    EEG.eyeBlink_rej = [];
    EEG.eyeBlink_rej_perc=0;
    
    %Low-pass filter data just to find EB
    EEG1  = pop_basicfilter(EEG,  1:exp.nEEGchans , 'Boundary', 'boundary', 'Cutoff', [10], 'Design', 'butter', 'Filter', 'lowpass', 'order', 4);

    i = 0;
    for e = 1:length(EEG.epoch)
        
        %Reject trials with eyeblinks between 1st pulse onset and 2nd pulse end
        t.p1_idx = str2double([EEG.epoch(e).eventtype(:)]) == 5;  %str2double, cell2mat
        if sum(t.p1_idx) == 0; t.p1_idx = cell2mat([EEG.epoch(e).eventtype(:)]) == 5; end
        assert(sum(t.p1_idx ~= 0)); 
        t.p1_time= cell2mat([EEG.epoch(e).eventlatency(t.p1_idx)]);
        t.p1_times = find([EEG.times] == t.p1_time(1)); %cell2mat([EEG.epoch(e).eventlatency(t.p1_idx)]);
        
        t.lo_idx= str2double([EEG.epoch(e).eventtype(:)]) == 12;
        if sum(t.lo_idx) == 0; t.lo_idx = cell2mat([EEG.epoch(e).eventtype(:)]) == 12; end
        assert(sum(t.lo_idx ~= 0));
        t.lo_time = cell2mat([EEG.epoch(e).eventlatency(t.lo_idx)]);
        t.lo_times = find([EEG.times] == t.lo_time(1)); %cell2mat([EEG.epoch(e).eventlatency(t.p1_idx)]);
        
        if EEG.sumData.npulses(e) == 1 | (EEG.sumData.npulses(e) == 2 & EEG.sumData.gapdur(e) ==0) % if one pulse, exclude all trials with blinks during pulses
            t.minMax_eog = abs(min(EEG1.data(exp.eogChans(1:2), t.p1_times:t.lo_times,e),[],2) - max(EEG1.data(exp.eogChans(1:2), t.p1_times:t.lo_times,e),[],2));
            t.minMax_eeg = abs(min(EEG1.data(exp.FrontChans(:), t.p1_times:t.lo_times,e),[],2) - max(EEG1.data(exp.FrontChans(:), t.p1_times:t.lo_times,e),[],2));

        elseif (EEG.sumData.npulses(e) == 2 & EEG.sumData.gapdur(e) ~=0); % if two pulses, exclude between pulses
            
            t.gap_idx = str2double([EEG.epoch(e).eventtype(:)]) == 9; 
            if sum(t.gap_idx) == 0; t.gap_idx = cell2mat([EEG.epoch(e).eventtype(:)]) == 9; end
            assert(sum(t.gap_idx ~= 0));
           
            t.gap_time= cell2mat([EEG.epoch(e).eventlatency(t.gap_idx)]);
           
            if isempty(t.gap_time)
                warning('No gap trigger found!')
            end
            
            t.gap_times = find([EEG.times] == t.gap_time(1)); %cell2mat([EEG.epoch(e).eventlatency(t.p1_idx)]);
            
            t.p2_idx= str2double([EEG.epoch(e).eventtype(:)]) == 8;
            if sum(t.p2_idx) == 0; t.p2_idx = cell2mat([EEG.epoch(e).eventtype(:)]) == 8; end
            assert(sum(t.p2_idx ~= 0));
            t.p2_time = cell2mat([EEG.epoch(e).eventlatency(t.p2_idx)]);
            t.p2_times = find([EEG.times] == t.p2_time(1)); %cell2mat([EEG.epoch(e).eventlatency(t.p1_idx)]);
        
            t.minMax_p1 = abs(min(EEG1.data(exp.eogChans(1:2), t.p1_times:t.gap_times,e),[],2) - max(EEG1.data(exp.eogChans(1:2), t.p1_times:t.gap_times,e),[],2));
            t.minMax_p2 = abs(min(EEG1.data(exp.eogChans(1:2), t.p2_times:t.lo_times,e),[],2) - max(EEG1.data(exp.eogChans(1:2), t.p2_times:t.lo_times,e),[],2));
            
            t.minMax_p1_eeg = abs(min(EEG1.data(exp.FrontChans(:), t.p1_times:t.gap_times,e),[],2) - max(EEG1.data(exp.FrontChans(:), t.p1_times:t.gap_times,e),[],2));
            t.minMax_p2_eeg = abs(min(EEG1.data(exp.FrontChans(:), t.p2_times:t.lo_times,e),[],2) - max(EEG1.data(exp.FrontChans(:), t.p2_times:t.lo_times,e),[],2));
            
            t.minMax_eog = [t.minMax_p1, t.minMax_p2]; 
            t.minMax_eeg = [t.minMax_p1_eeg,t.minMax_p2_eeg];
        end

        if ~isempty(mean(t.minMax_eog > 100)) & ~isempty(find(t.minMax_eeg > 70)) %~isempty(find(t.minMax_eeg > 70)) % if both eye channels and frontal eeg channels show artefacts, reject
            EEG.eyeBlink_rej = [EEG.eyeBlink_rej; e];
            EEG.eyeBlink_rej_perc= length(EEG.eyeBlink_rej)/length(EEG.epoch)*100;
            i = i+1;
            disp(['Artefact #' num2str(i) ' in epoch #' num2str(e) 'rejected!']);
        else
            EEG.eyeBlink_rej = EEG.eyeBlink_rej;
        end
        clear t;
    end
    
end

disp(['Percentage trials w/eyeblink during pulses: ' num2str(EEG.eyeBlink_rej_perc)])

end