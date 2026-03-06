%% GT_Figure2
% Plots grand-averaged traces of MBL, CPP & occipital alpha in Figure 2. 
% Also plots bottom panel in Fig. S9


%% MBL - Figure 2A
for ee = 1:2 %For both experiments
    exp.thisExp = ee; %Experiment to analyse
    [exp, plt]  = GT_setup(exp);
    exp.sub_id{1}            = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; % P8 excluded because only 4 blocks were recorded.
    exp.sub_id{2}            = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23]; % P1 & P21 excluded because they could not perform task well

    %Load dummy EEG to have scalp locations
    EEG = pop_loadset([exp.filepath '/csd_preproc_P1_' exp.thisExpCode '_P2.set']);
    %Load ERP average data
    load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_TFA_average_beta_np_gap_dir_acc_.mat'])
    
    plt.epoch = 'P1'; %Locked to stimulus onset
    plt.font = 11;
    par.db = 1; %DB transform power data
    ds = (200/512); %Downsample (MBL data are 200Hz)
    thisGrouping = 'np_gap_dir_acc_' %Data are sorted by gap duration, dot direction and accuracy
    if par.db == 0, ylims = [-2.75 2.75]; else ylims = [-2.75,2.75]; end %Plotting parameters
    e=1;
    %%
    f = figure;
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 11 9];
    factor = 2; %Beta data were extracted as amplitude; need to convert to power by squaring values
    chanSel = [3];
    for acc = 1:2
        for g = 1:length(exp.gapDur{exp.thisExp}) %Longest gap only!
            for d = 1:2
                thisGap = exp.gapLabs{exp.thisExp}{g};
                thisAcc = exp.accLabs{acc};
                thisDir = exp.dirLabs{d};
                
                thisCol = plt.satPalettes{2,g};
                
                if acc == 1; lineType = '--'; else lineType = '-'; end
                
                %Compute lateralisation based on cluster of left/right hemisphere electrodes
                dataTag = ['twoPulse' thisGap  thisDir thisAcc];
                
                if par.db == 1; unitTag = '(dB)';
                    baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:),2).^factor;
                    data = 10*log10(allP.(plt.epoch).(thisGrouping).(dataTag).sl.^factor./baselines);
                else; unitTag = '(\muV/m^{2})';
                    baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:),2).^factor;
                    data =(allP.(plt.epoch).(thisGrouping).(dataTag).sl).^factor;%-baselines);
                end
                newData(d,:,:) = squeeze(nanmean(data([exp.movCluster(1,chanSel)],:,:),1))-squeeze(nanmean(data([exp.movCluster(2,chanSel)],:,:),1)); %L-R
                dataTopo(d,g,:,:,:) = squeeze(nanmean(data(:,1.5*200:2.6*200,:),2));
                
                endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.70+1.5]*200); %Plot up to the earlies time the lead-out dots might have occurred 
            end
            %Subtract (Right - left) based on dot motion
            if g == length(exp.gapDur{exp.thisExp}) ; newdataTopo = squeeze(nanmean(dataTopo(1,:,:,:,:),2) - nanmean(dataTopo(2,:,:,:,:),2)); end
            
            %Low pass filter at 6Hz for viz. 
            [bfilt,afilt] = butter(3, 6*2/200, 'low');   % lo-pass
            plotData = squeeze(nanmean(newData(1,1.3*200:endPoint,:)-newData(2,1.3*200:endPoint,:),1))';
            for s = 1:size(plotData,1) %Subject
                plotData(s,:) = filtfilt(bfilt,afilt, plotData(s,:));
            end
            
            if g == 1; fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.2); end; %Plot pulses in time
            if acc == 1; fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2); end%Plot pulses in time
            hold on
            stdshade(plotData, plt.alpha, thisCol, lineType,[],plt.smooth) %Plot ERP traces
            hold on;
            
            permData = [];
            
            % Add cluster statistics for difference from 0
            for s = 1:size(plotData,1)
                permData.mbl{s} = plotData(s,:);
            end
            exp.condDur = size(plotData,2)/200;
            cluster = GT_clusterMBL(exp,permData);
            
            clusterStats(g,acc,:) = repmat(0,1,3000);
            if isfield(cluster.stat1, 'posclusters')
                for i = 1:size(cluster.stat1.posclusters,2)
                    if i >0
                        if cluster.stat1.posclusters(i).prob < 0.05
                            clusterStats(g,acc,find(cluster.stat1.posclusterslabelmat == i)) = 1;
                        end
                    end
                end
            end
            if isfield(cluster.stat1, 'negclusters')
                for i = 1:size(cluster.stat1.negclusters,2)
                    if i >0
                        if cluster.stat1.negclusters(i).prob < 0.05
                            clusterStats(g,acc,find(cluster.stat1.negclusterslabelmat == i)) = 1;
                        end
                    end
                end
            end
            
            xticks(round(plt.xticks(e,:)*ds));
            xticklabels(plt.xlabs(e,:)); set(gca, 'fontsize', plt.font)
            yline(0, 'k'); xline(0.2*200, 'k')
            clear newData;
        end
    end
    % Plot cluster stats
    for acc = 2
        if acc == 1; fac = -1; else fac = 1; end %Flip sign to plot sig. stars if error
        for g = 1:length(exp.gapDur{exp.thisExp}) 
            thisCol = plt.satPalettes{2,g};
            
            ypos = -2.3*fac-0.1*g*fac;
            x = find(clusterStats(g,acc,:)); hold on;
            plot(x, repmat(ypos,1,length(x)), '*', 'color', thisCol, 'markersize',1.5); hold on;
        end
    end
    ylim(ylims)
    ylabel(['MBL ' unitTag ])
    xlabel('Time (s)')
    
    %Insert topography
    axes('position',[.1 .65 .25 .25]);
    
    GT_topoplot(nanmean(newdataTopo,2),EEG.chanlocs, ...
        'numcontour', 0, ...
        'plotrad',0.65, ...
        'headrad', 0.5,...
        'emarker2', {[cat(2,exp.movCluster(1,chanSel),exp.movCluster(2,chanSel))],'*','k',3,1}); % To plot location of the selected channel/s
    caxis([ylim]);
    
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 11 9];
    
    exportgraphics(f,[exp.figPath, 'GT_Fig2A_Exp' num2str(ee) '.tiff'], 'Resolution', 600)
    clear all;
end

%% CPP
% Plot CPP across 4 conditions, for ALL delays
for ee = 2 %For both experiments
    exp.thisExp = ee; %Experiment to analyse
    [exp, plt]  = GT_setup(exp);
    exp.sub_id{1}            = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; % P8 excluded because only 4 blocks were recorded.
    exp.sub_id{2}            = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23]; % P1 & P21 excluded because they could not perform task well

    %% Ploting parameters 
    ds = 1; sr = 512; %No downsampling
    xlims = [1.3,4.5]*200; plt.epoch = 'P1'; e = 1;     plt.font = 11;
    plt.baseline = 0; ylims = [-5 20];
    f = figure
    data_allC = [];
    plt.RIDE = 0;    %Possibility to clean colour-evoked VEPs using RIDE (not used in paper)
    plt.minusP1 = 0; %Remove P1 activity to compute CPP-P2/P1 amplitude ratio by subtracting P1 activity from longest gap condition (not used for this figure)
    
    %%Use if loading the file with ALL channels
    chans = exp.cppChannels; %CPP channels: [4,32,31,30,21,17,18,5,19,20];
%     chans = 1; %Because we load the file with selected channels only 
    
    %% Plot
    %Load dummy EEG to have scalp locations
    EEG = pop_loadset([exp.filepath '/csd_preproc_P1_' exp.thisExpCode '_P2.set']);
    
    %Load sample-wise regression to remove colour decisions if plt.correction
    %[not used in paper]
%     RIDE = load([exp.filepath '/' exp.thisExpCode '_allP_ride_allPdata.mat']);
%     RIDE = RIDE.allP;
    
    %Load ERP average data
    load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_ERP_averageData_cpp.mat'])
    
    %%
    peak_cpp1 = []; peak_cpp2 = [];
    for g = 1:length(exp.gapDur{exp.thisExp}) %Longest gap only!
        p1Idx = []; p2Idx = [];
        thisGrouping= 'np_gap_'; %Get correctly sorted data
        thisGap = exp.gapLabs{exp.thisExp}{g};
        gapDur = exp.gapDur{exp.thisExp}(g);
        
        longGap = exp.gapLabs{exp.thisExp}(length(exp.gapDur{exp.thisExp}));
        
        thisCol = plt.satPalettes{2,g};
        
        dataTag = ['twoPulse' thisGap];
        longTag = ['twoPulse' longGap{1}];
        
        rideTag = thisGap;
        baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,round(1.3*sr):round(1.5*sr),:),2);
        data =((allP.(plt.epoch).(thisGrouping).(dataTag).sl)-baselines);
%         baselines = nanmean(cppAll.(plt.epoch).(thisGrouping).(dataTag).sl(:,round(1.3*sr):round(1.5*sr),:),2);
%         data =((cppAll.(plt.epoch).(thisGrouping).(dataTag).sl)-baselines);
        
        %To plot topography at peak of gap = 0;
        if g == 1;
            dataTopo = squeeze(nanmean(data(:,round(2*sr):round(2.1*sr),:),2)-nanmean(data(:,round(1.4*sr):round(1.5*sr),:),2));
        end
        
        data = squeeze(nanmean(data([chans],:,:),1));
        dataLong = ((allP.(plt.epoch).(thisGrouping).(longTag).sl)-baselines);
        dataLong = squeeze(nanmean(dataLong([chans],:,:),1));
%         dataLong = ((cppAll.(plt.epoch).(thisGrouping).(longTag).sl)-baselines);
%         dataLong = squeeze(nanmean(dataLong([chans],:,:),1));
        
        if plt.RIDE ==1
            indices = round(1.3*512:4.499*512);
            for id = 1:size(RIDE.(rideTag).s1,1)
                thisStim(id,:) = nanmean(RIDE.(rideTag).s1{id}(:,[chans]),2);
            end
            data(indices,:) = data(indices,:) - thisStim';
        end
        
        p1Idx = round([1.9*sr:2.1*sr ]); %For peak amp, 0.4-0.6s after P1 onset
        p2Idx = round([(1.9 + 0.2 +gapDur)*sr:(2.1+0.2+gapDur)*sr]); %For peak amp, 0.4s-0.6 after P2 onset
        
        %Low pass filter at 6Hz for viz.
        [bfilt,afilt] = butter(3, 6*2/512, 'low');   % lo-pass
        for s = 1:size(data,2) %Subject
            data(:,s) = filtfilt(bfilt,afilt, data(:,s));
            dataP2(:,s) = data(:,s); %dataP1 will be the data after removing P1 activity based on longest gap condition
            
            dataLong(:,s) = filtfilt(bfilt,afilt, dataLong(:,s));
            dataP2(1:(1.5+1.4)*sr,s) = data(1:(1.5+1.4)*sr,s)-dataLong(1:(1.5+1.4)*sr,s);
            peak_cpp1(g,s) = nanmean(data(p1Idx,s));
            peak_cpp2(g,s) = nanmean(data(p2Idx,s));
            
            baseline1 = 0; %data are baselined already
            baseline2 = nanmean(dataP2(p2Idx(1)-0.5*sr:p2Idx(1)-0.4*sr,s));
            baseline2_uncorr(s,:) = nanmean(data(p2Idx(1)-0.45*sr:p2Idx(1)-0.4*sr,s));
            
            baseIdx = (baseline2(1));
            % %         xline(p2Idx(1)-(0.45+1.3)*sr, '--r');   hold on; xline(p2Idx(1)-1.3*sr, '--b');
            exc_cpp1(g,s) = nanmean(data(p1Idx(:),s))-baseline1;
            exc_cpp2(g,s) = nanmean(dataP2(p2Idx(:),s))-baseline2;
            exc_cpp2_unc(g,s) = nanmean(data(p2Idx(:),s))-baseline2_uncorr(s,:); %Compute P2 amplitude on raw potentials
        end
        
        %Get grand-average amplitudes 
        baseline2 = nanmean(nanmean(dataP2(p2Idx(1)-0.5*sr:p2Idx(1)-0.4*sr,:),1),2);
        
        ga_exc_cpp1(g,:) = nanmean(nanmean(data(p1Idx(:),:),2),1);
        ga_exc_cpp2(g,:) = nanmean(nanmean(dataP2(p2Idx(:),:),2),1)-baseline2;
        ga_exc_cpp2_uncorr(g,:) = nanmean(nanmean(data(p2Idx(:),:),2),1)-nanmean(baseline2_uncorr);
        
        %Remove data after median RTs post dot offset,%conservative
        endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.75+1.5]*sr);
        
        fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2); %Plot pulses in time
        hold on
        if plt.minusP1 == 1 %Subtract CPP-P1 based on longest gap?
            stdshade(dataP2(1.3*sr:endPoint,:)', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
        else
            stdshade(data(1.3*sr:endPoint,:)', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
        end
        
        hold on;
        
        xlim([1,2.5*sr])
        permData = [];
        % Add cluster statistics for difference from 0
        for s = 1:size(data(1.3*sr:endPoint,:)',1)
            permData.cpp{s} = data(1.3*sr:endPoint,s)';
        end
        exp.condDur = size(data(1.3*sr:endPoint,s),1)/sr;
        cluster = GT_clusterCPP(exp,permData);
        
        clusterStats(g,:) = repmat(0,1,3000);
        if isfield(cluster.stat1, 'posclusters')
            for i = 1:size(cluster.stat1.posclusters,2)
                if i >0
                    
                    if cluster.stat1.posclusters(i).prob < 0.05
                        clusterStats(g,find(cluster.stat1.posclusterslabelmat == i)) = 1;
                    end
                end
            end
        end
        if isfield(cluster.stat1, 'negclusters')
            for i = 1:size(cluster.stat1.negclusters,2)
                if i >0
                    if cluster.stat1.negclusters(i).prob < 0.05
                        clusterStats(g,find(cluster.stat1.negclusterslabelmat == i)) = 1;
                    end
                end
            end
        end
        
        hold on
        thisCol = plt.colors(1)% plt.colors(7,2);
        if g == 1; fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], thisCol{1}, 'LineStyle','none','facealpha', 0.2); end %Plot pulses in time
        ylim(ylims);
        hold on;
        xline(plt.zero(e)*exp.srate*ds);
        yline(0, 'LineStyle', '--');
        
        xticks(round(plt.xticks(e,:)*ds));
        xticklabels(plt.xlabs(e,:)); set(gca, 'fontsize', plt.font)
        clear data;
        %     xlim([1, endPoint-1.3*sr])
    end
    ylabel('CPP (\muV/m^{2})')
    xlabel('Time (s)')
    
    xline(0.75*sr, '--k')
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 11 9];
    
    % Plot cluster stats
    for g = 1:length(exp.gapDur{exp.thisExp}) %Longest gap only!
        thisCol = plt.satPalettes{2,g};
        
        ypos = -2.5-0.5*g;
        x = find(clusterStats(g,:)); hold on;
        plot(x, repmat(ypos,1,length(x)), '*', 'color', thisCol, 'markersize',1.5); hold on;
        
    end
    
    %Insert topography
    axes('position',[.1 .65 .25 .25]);
    
    GT_topoplot(nanmean(dataTopo,2),EEG.chanlocs, ...
        'numcontour', 0, ...
        'plotrad',0.65, ...
        'headrad', 0.5,...
        'emarker2', {[chans],'*','k',3,1}); % To plot location of the selected channel/s
    caxis([-20 20]);
    
%     exportgraphics(f,[exp.figPath, 'GT_Fig2B_Exp' num2str(ee) '.tiff'], 'Resolution', 600)
    
    %t-test for CPP-P1 vs CPP-P2 max amplitudes, averaging across all non-zero gaps.
    gacpp1 = nanmean(exc_cpp1(2:end,:)); gacpp2 = nanmean(exc_cpp2(2:end,:)); gacpp2_unc = nanmean(exc_cpp2_unc(2:end,:))
    [H,P,CI,stats] = ttest(gacpp1,gacpp2)
    [H,P,CI,stats] = ttest(gacpp1,gacpp2_unc)
    
    %Compute what % of CPP-P1 is CPP-P2
    percCPP = ga_exc_cpp2/ga_exc_cpp1*100;
    percCPP_unc = ga_exc_cpp2_uncorr/ga_exc_cpp1*100;
    
    disp(['Exp. ratio:' num2str(nanmean(percCPP(2:end,exp.thisExp)))]) %Average over non-zero gaps, after removing P1 overlap
    disp(['Exp.  (unc) ratio:' num2str(nanmean(percCPP_unc(2:end,exp.thisExp)))]) %Average over non-zero gaps, on raw data
    
    clear all;

end

%% Occipital alpha
for ee = 1:2 %For both experiments
    exp.thisExp = ee; %Experiment to analyse
    [exp, plt]  = GT_setup(exp);
    exp.sub_id{1}            = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; % P8 excluded because only 4 blocks were recorded.
    exp.sub_id{2}            = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23]; % P1 & P21 excluded because they could not perform task well

%Load TFA average data
load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_TFA_average_alpha.mat'])

ds = (200/512);
clear diffData;
xlims = [1.3,4.5]*200;
plt.epoch = 'P1'; e = 1;
plt.baseline = 0;
par.db = 1;
factor = 2; %To convert beta amplitude values to power
sr = 200;
plt.font = 11;
data_allC = [];
%%
occChans = [9,10,17,16,21,30,29,38,39,22]
ylims = [-2.5,3];
f = figure
for g = 1:length(exp.gapDur{exp.thisExp}) %Longest gap only!
    
    thisGrouping= 'np_gap_';
    thisGap = exp.gapLabs{exp.thisExp}{g};
    
    thisCol = plt.satPalettes{2,g};
    %Compute lateralisation based on cluster of left/right hemisphere electrodes
    dataTag = ['twoPulse' thisGap];
    
    if par.db == 1
        baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:).^2,2);
        data = 10*log10((allP.(plt.epoch).(thisGrouping).(dataTag).sl).^2./baselines);
    else
        baselines = 0;%nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:),2);
        data =((allP.(plt.epoch).(thisGrouping).(dataTag).sl).^2-baselines);
    end
    data_allC(1,:,:) = squeeze(nanmean(data(occChans,:,:)));
    if g == 1;
        dataTopo = squeeze(nanmean(data(:,round(2*200):round(2.1*200),:),2)-nanmean(data(:,round(1.4*200):round(1.5*200),:),2));
    end %To plot topography at peak of gap = 0;
    
    data = squeeze(nanmean(data(occChans,:,:),1));
    %Remove data after median RTs post dot offset,conservative
    endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.75+1.5]*200);
    
    
    permData = [];
    % Add cluster statistics for difference from 0
    for s = 1:size(data(1.3*sr:endPoint,:)',1)
        permData.mbl{s} = data(1.3*sr:endPoint,s)';
    end
    exp.condDur = size(data(1.3*sr:endPoint,s),1)/sr;
    cluster = GT_clusterMBL(exp,permData);
    
    %save
    clusterStats(g,:) = repmat(0,1,3000);
    if isfield(cluster.stat1, 'posclusters')
        for i = 1:size(cluster.stat1.posclusters,2)
            if i >0
                
                if cluster.stat1.posclusters(i).prob < 0.05
                    clusterStats(g,find(cluster.stat1.posclusterslabelmat == i)) = 1;
                end
            end
        end
    end
    if isfield(cluster.stat1, 'negclusters')
        for i = 1:size(cluster.stat1.negclusters,2)
            if i >0
                if cluster.stat1.negclusters(i).prob < 0.05
                    clusterStats(g,find(cluster.stat1.negclusterslabelmat == i)) = 1;
                end
            end
        end
    end
        
    fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2); %Plot pulses in time
    hold on
    stdshade(data(1.3*200:endPoint,:)', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
    hold on;
    
    ylabel('Occipital alpha')
    xlabel('Time (s)')
    
    hold on
    thisCol = plt.colors(1)% plt.colors(7,2);
    if g == 1; fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], thisCol{1}, 'LineStyle','none','facealpha', 0.2); end %Plot pulses in time
    ylim(ylims);
    hold on;
    xline(plt.zero(e)*exp.srate*ds);
    yline(0, 'LineStyle', '--');
    
    xticks(round(plt.xticks(e,:)*ds));
    xticklabels(plt.xlabs(e,:)); set(gca, 'fontsize', plt.font)
%     set(gca, 'ydir', 'reverse')
    
    allAlpha{exp.thisExp}(g,:) = nanmean(data,2);
    clear data;
end

% Plot cluster stats
for g = 1:length(exp.gapDur{exp.thisExp}) %Longest gap only!
    thisCol = plt.satPalettes{2,g};
    
    ypos = -2-0.1*g;
    x = find(clusterStats(g,:)); hold on;
    plot(x, repmat(ypos,1,length(x)), '*', 'color', thisCol, 'markersize',1.5); hold on;
    
end

xlabel(plt.xlabel); ylabel(['Occipital alpha (dB)']);
xline(0.75*200, '--k')
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 11 9];

%Insert topography
axes('position',[.1 .65 .25 .25]);

GT_topoplot(nanmean(dataTopo,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {occChans,'*','k',3,1}); % To plot location of the selected channel
caxis([-3 3]); %colorbar;


filename = ['Fig2_' exp.thisExpCode '_alpha.tif'];
saveas(f, [exp.figurepath, '/ga/' filename]);

exportgraphics(f,[exp.figPath, 'GT_Fig2C_Exp' num2str(ee) '.tiff'], 'Resolution', 600)

end
%% Occipital alpha E1 vs E2 (Sup. Fig.S9)
%Load TFA average data
occChans = [9,10,17,16,21,30,29,38,39,22]

plt.font=11;
% f = figure
factor = 2; 
clear dataTopo; 
for e = 2
    hold on;
    exp.thisExp = e; %Experiment to analyse
    [exp, plt]  = GT_setup(exp);
     
    load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_TFA_average_alpha.mat'])
    
    ds = (200/512);
    clear diffData;
    xlims = [1.3,4.5]*200;
    plt.epoch = 'P1'; 
    plt.baseline = 0; ep = 1;
    par.db = 1;
    ylims = [-2,3];
    data_allC = [];
    plt.font = 11;
    
    for g = length(exp.gapDur{exp.thisExp}) %Longest gap only!
        
        thisGrouping= 'np_gap_';
        thisGap = exp.gapLabs{exp.thisExp}{g};
        
        thisCol = plt.satPalettes{5,exp.thisExp*2};
        %Compute lateralisation based on cluster of left/right hemisphere electrodes
        dataTag = ['twoPulse' thisGap];
        
        if par.db == 1
            baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:).^factor,2);
            data = 10*log10((allP.(plt.epoch).(thisGrouping).(dataTag).sl).^factor./baselines);
        else
            baselines = 0;%nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:),2);
            data =((allP.(plt.epoch).(thisGrouping).(dataTag).sl).^factor-baselines);
        end
        data_allC(1,:,:) = squeeze(nanmean(data(occChans,:,:)));
 
        dataTopo{e} = squeeze(data(:,:,:)); %Keep all channels to plot correct traces at the right channels 
        
        data = squeeze(nanmean(data(occChans,:,:),1));
        endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.75+1.5]*200);
        
        hold on;
        yy = xline([exp.gapDur{exp.thisExp}(g)+0.4]*200);
        yy.Color = thisCol;
        
        %     fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2); %Plot pulses in time
        hold on
        stdshade(data(1.3*200:endPoint,:)', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
        hold on;
        
        %Extract pre-P2 power
        preP2 = round([(exp.gapDur{exp.thisExp}(g)+1.5-0.1)*200: (exp.gapDur{exp.thisExp}(g) +1.5)*200]);
        p2Pow_m{e} = squeeze(nanmean(data(preP2,:)));
        p2Pow_sem{e} = squeeze(nanstd(nanmean(data(preP2,:)))./sqrt(size(data,2)));
        
        ylabel('Occipital alpha')
        xlabel('Time (s)')
        
        hold on
        thisCol = plt.colors(1)% plt.colors(7,2);
        ylim(ylims);
        hold on;
        xline(plt.zero(ep)*exp.srate*ds);
        
        xticks(round(plt.xticks(ep,:)*ds));
        xticklabels(plt.xlabs(ep,:)); set(gca, 'fontsize', plt.font)
%         set(gca, 'ydir', 'reverse')
        
        allAlpha{exp.thisExp}(g,:) = nanmean(data,2);
        clear data;
    end

    xlabel(plt.xlabel); ylabel(['Occipital alpha (dB)']);
    % xline(0.75*200, '--k')
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 11 9];
    
end
box on
    save(['GT_Alpha_E1_vs_E2.mat'], 'dataTopo'); 

filename = ['Fig2_' exp.thisExpCode '_alpha_E1E2.tif'];
saveas(f, [exp.figurepath, '/ga/' filename]);

%%
%% Plot
f=figure;
for e = 1:2
        thisCol = plt.satPalettes{5,e*2};
     [exp, plt]  = GT_setup(exp);
    g = length(exp.gapDur{e});
 
    ds = (200/512);
    xlims = [1.3,4.5]*200;
    plt.epoch = 'P1'; 
    plt.baseline = 0; ep = 1;
    par.db = 1;
    ylims = [-3,2];
    data_allC = [];
    plt.font = 11;
    
     if e == 1; 
        rectangle('Position',[preP2(1)-1.3*200 ylims(1) 0.1*200 ylims(end)-ylims(1)], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineStyle', '--'); hold on; 
%         fill(,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.2); 
    end %Plot pulses in time; end
    
    stdshade(squeeze(nanmean(dataTopo{e}(occChans,1.3*200:endPoint,:),1))', plt.alpha, thisCol, '-',[],plt.smooth) %Plot ERP traces
    data = squeeze(nanmean(dataTopo{e}(occChans,:,:),1));
    hold on; 
    xticks(round(plt.xticks(ep,:)*ds));
    xticklabels(plt.xlabs(ep,:)); set(gca, 'fontsize', plt.font)
    set(gca, 'ydir', 'reverse')
    
    xlabel(plt.xlabel); ylabel(['Occipital alpha (dB)']);
    ds = (200/512);
    clear diffData;
    
    %Extract pre-P2 power
    preP2 = round([(exp.gapDur{1}(3)+1.5+0.2-0.1)*200: (exp.gapDur{1}(3) +1.5+0.2)*200]);
    p2Pow_m{e} = squeeze(nanmean(data(preP2,:)));
    p2Pow_sem{e} = squeeze(nanstd(nanmean(data(preP2,:)))./sqrt(size(data,2)));
%     xline([preP2(1),preP2(end)]-1.3*200, 'color', thisCol);
     
    xline(plt.gapShades{e}(g+1,1)*ds, 'color', thisCol, 'linewidth', 2);
   

    xline(0.2*200, 'k');
    f.Units = 'centimeters';
    f.OuterPosition = [0.25 0.35 11 9];
end
[H,P,CI,stats] = ttest2(p2Pow_m{1}, p2Pow_m{2})
