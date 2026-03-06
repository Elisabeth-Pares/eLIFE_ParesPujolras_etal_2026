function [GA, allP] = GT_Fig3_fitCPP_coh(exp,par,id, GA, allP)

%Regress CPP-P1 and CPP-P2 against P1coh and P2coh
%Is CPP build-up modulated by pulse coherence?
%CPP-P1/2(t) ~ p1coh + p2coh 

sub = exp.sub_id{exp.thisExp}(id);

%% Load CPP data 
   
% Which channels to use?
theseChans = exp.parietalFocus;

% Load EEG data & behaviour
filename = ['csd_preproc_P1_' exp.thisExpCode '_P' num2str(sub) '.set'];
EEG = pop_loadset([exp.filepath '\' filename]);

[bfilt,afilt] = butter(3, 6*2/EEG.srate, 'low');   % lo-pass - makes regressions more robust
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

%Load EEG data & behaviour
allP.sumData{id} = EEG.sumData;

timeVector = -1.5:1/200:3;
timeRange = find(timeVector>-0.7 & timeVector < 2.2);


%Align all gap conditions to P1: CPP-P1
CPP_P1(:,:,:) = zCPP(:,1.3*EEG.srate:2.5*EEG.srate,:);
  
CPP_P2= CPP_P1;
CPP_P2(:) = NaN;

if strcmp(par.modelName, 'CPP_P1'), gStart = 1;             %Include all conditions
elseif contains(par.modelName, 'CPP_P2'), gStart = 2; end   %Exclude 0 gap condition
%% Align data to P2, take each gap separately
for g = gStart:length(exp.gapDur{exp.thisExp})
    
    thisGap = exp.gapDur{exp.thisExp}(g);
    idx = EEG.sumData.gapdur == thisGap & EEG.sumData.npulses == 2;
    gapLab = exp.gapLabs{exp.thisExp}{g};
    
    %Baseline CPP before P2
    if par.baseline == 1
        p2Idx = round((1.5+0.2+thisGap)*EEG.srate);
        baseline = (nanmean(zCPP(:,p2Idx-0.1*EEG.srate:p2Idx,idx),2)); %retried with 200 ms baseline...changed from -0.15*200
        zCPP(:,p2Idx-0.2*EEG.srate:end,idx) = zCPP(:,p2Idx-0.2*EEG.srate:end,idx) - baseline;
    else
        p2Idx = round((1.5+0.2+thisGap)*EEG.srate);
        baseline = 0;
        zCPP(:,p2Idx:end,idx) = zCPP(:,p2Idx:end,idx) - baseline;
    end
    CPP_P2(:,:,idx) = zCPP(:,p2Idx-0.2*EEG.srate:(p2Idx+1*EEG.srate),idx);
end
    
c = 0;
%Time resolved
for t = 1:size(CPP_P1,2)
    c = c+1;
    
    X = [nanzscore(EEG.sumData.p1coh(:)),... %p1coh
        nanzscore(EEG.sumData.p2coh(:)),... %p2coh
        ];

    regNames = {'p1coh', 'p2coh'};

    if strcmp(par.modelName, 'CPP_P1')
        Y = squeeze(nanmean(CPP_P1(:,t,:),1));
    elseif strcmp(par.modelName, 'CPP_P2')
        Y = squeeze(nanmean(CPP_P2(:,t,:),1));
    elseif strcmp(par.modelName, 'CPP_P2_xCPP1amp') %R1: Consider interaction of p1coh*CPP-P1amp interaction
        Y = squeeze(nanmean(CPP_P2(:,t,:),1));
        
        peakP1 = round(0.6*512:0.7*512);
        CPP_P1amp = nanzscore(squeeze(nanmean(CPP_P1(:,peakP1,:),2)));

        X = [nanzscore(EEG.sumData.p1coh(:)),... 
        nanzscore(EEG.sumData.p2coh(:)),... 
        nanzscore(CPP_P1amp),...
        nanzscore(EEG.sumData.p1coh(:)).*CPP_P1amp,... 
        nanzscore(EEG.sumData.p2coh(:)).*CPP_P1amp,... 
        ];
        regNames = {'p1coh', 'p2coh', 'cppAmp','p1cohxcpp', 'p2cohxcpp'};
    end

    X = nanzscore(X);
    
    m = regstats(nanzscore(Y),nanzscore(X),'linear',{'beta','tstat','rsquare', 'r'}); 
    
    % Save beta coefficients and residuals
    for r = 1:length(regNames)
        thisReg = [regNames{r}];
        GA.(thisReg).data(id,:,c) = m.beta(r+1); %Beta coefficients
    end
    if size(m.r,1) == size(CPP_P1,3)
        GA.(thisReg).residuals{id}(:,c) = m.r; %Residuals
    else
        GA.(thisReg).residuals{id}(:,c)= NaN; 
    end
    GA.time = timeVector(timeRange);
      
end



end