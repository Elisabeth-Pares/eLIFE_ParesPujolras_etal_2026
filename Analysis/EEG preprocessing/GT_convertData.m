function [] = GT_convertData(sub,exp, thisExp)

allEEG = [];

fullPath = ([exp.database, '/P' num2str(sub) '/']);

files = dir(fullPath);
files = files(contains({files.name}, '.bdf'))
files = {files.name}

full_EEG = [];
missingEEG = []; extraEEG = [];

% Exp 1: one single file per participant including 10 blocks
% Exp 2: 14 files (blocks) per participant that need to be concatenated into a single file. 

for f = 1:length(files)
    %Convert data from .bdf (BioSemi) format to .set (EEGlab) format
    if exp.thisExp == 1
        thisFile = ['EP1_P' num2str(sub) '.bdf'];
    elseif exp.thisExp == 2
        thisFile = ['EP2_P' num2str(sub) '_B' num2str(f) '.bdf'];
    end
    cd(fullPath)
    
    %Deal with special participants 
    if exp.thisExp == 2 & sub == 2 && f == 10 %Deal with special subjects; 
            EEG = pop_fileio([fullPath 'EP2_P' num2str(sub) '_B15.bdf';], 'dataformat', 'auto');
    else
        EEG = pop_fileio([fullPath thisFile], 'dataformat', 'auto');
    end
    EEG = eeg_checkset(EEG);
    
    %Detect missing trials to align w/behavioural data. This could happen
    %if block recording started after participants initiated a block.
    if exp.thisExp == 1 
        if sub == 2, extraEEG = 1;   %One extra trial before experiment onset.
            %Remove all events before trigger 29
            idx = find([EEG.event.type] == 29,1,'first')
            [EEG.event(1:idx)] = []; %Delete up to that.
            [EEG.urevent(1:idx)] = []; %Delete up to that.
        end
        %In sub 5, recording started late and there are 28 missing trials
        %at the beginning of the file (block 1)
        if sub == 5, missingEEG = [1:28]; end
    
    elseif exp.thisExp == 2 %Sanity check: make sure all blocks have 70 trials 
        trialCount(f) = sum([EEG.event.type] == 1);
        if trialCount(f) < 70 & f~=14 %Becausae the last block has fewer trials 
            missingEEG = [missingEEG, (1:(70-trialCount(f)))+(exp.blockLength*(f-1))];
        elseif trialCount(f) > 70
            extraEEG = [extraEEG, (1:trialCount(f)-70)+(exp.blockLength*(f-1))];
        end
    end
    
    
    if isempty(full_EEG)
        full_EEG = EEG;
    else
        full_EEG = pop_mergeset(full_EEG, EEG);
    end
end

if exp.thisExp == 2 & any(trialCount(1:13) ~= 70), warning('wrong trial numbers!'); end %Check that all blocks have 70 trials
full_EEG.missingEEG = missingEEG; 
full_EEG.extraEEG = extraEEG; 

%Add metadata 
full_EEG = GT_addMetadata(sub, exp, full_EEG)
 
%Check trigger alignment 
[full_EEG] = GT_alignmentCheck(sub,exp, full_EEG) 

%Make sure event field is numeric 
if isa(EEG.event(2).type, 'char')
    [EEG.event.type] = str2double({EEG.event.type})';
end

%Save file 
if ~isempty(full_EEG)
    filename = ['fc_' exp.thisExpCode '_P' num2str(sub)];
    EEG = pop_saveset(full_EEG, filename, exp.convpath);
end

end
