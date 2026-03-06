%% Occipital alpha
%Load TFA average data
load(['C:\Users\elisa\Desktop\GapsTask\' exp.thisExpName '\Data\plot data\' exp.thisExpCode '_P1_allP_TFA_average_alpha.mat'])

%Load dummy EEG to have scalp locations
EEG = pop_loadset([exp.filepath '/csd_preproc_P1_' exp.thisExpCode '_P2.set']);

ds = (200/512);
clear diffData;
xlims = [1.3,4.5]*200;
plt.epoch = 'P1'; e = 1;
plt.baseline = 0; plt.font = 11;
par.db = 1;
factor = 2; %To convert beta amplitude values to power
sr = 200;
data_allC = [];

ylims = [-2.5,3]; dataTopo = [];
f = figure

occChans = [9,10,17,16,21,30,29,38,39,22]

clusterStats = [];
for g = length(exp.gapDur{exp.thisExp}) %Longest gap only!
        permData = [];

    for p1 = 1:2
    thisGrouping= 'np_gap_p1coh_';
    thisGap = exp.gapLabs{exp.thisExp}{g};
    thisP1 = exp.cohLabs{p1};
    
    thisCol = plt.satPalettes{1,p1*p1};
    %Compute lateralisation based on cluster of left/right hemisphere electrodes
    dataTag = ['twoPulse' thisGap thisP1];
    
    if par.db == 1
        baselines = nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:).^2,2);
        data = 10*log10((allP.(plt.epoch).(thisGrouping).(dataTag).sl).^2./baselines);
    else
        baselines = 0;%nanmean(allP.(plt.epoch).(thisGrouping).(dataTag).sl(:,[1.3*200:1.5*200],:),2);
        data =((allP.(plt.epoch).(thisGrouping).(dataTag).sl).^2-baselines);
    end
    data_allC(1,:,:) = squeeze(nanmean(data(occChans,:,:)));
    
    p2Idx = (0.2 + exp.gapDur{exp.thisExp}(g) + 1.5)*200;
    dataTopo(p1,:,:) = squeeze(nanmean(data(:,p2Idx-0.2*200:p2Idx,:),2));
    
    if p1 == 2; ddataTopo = dataTopo(1,:,:) - dataTopo(2,:,:); dataTopo = []; dataTopo = squeeze(ddataTopo); clear ddataTopo; end; 
    data = squeeze(nanmean(data(occChans,:,:),1));
    %Remove data after median RTs post dot offset,conservative
    endPoint = round([0.4 + exp.gapDur{exp.thisExp}(g) + 0.75+1.5]*200);
    
    if p1 == 2; fill(plt.gapShades{exp.thisExp}(1,:)*ds,[ylims, flip(ylims)], 'k', 'LineStyle','none','facealpha', 0.2);  hold on;
        fill(plt.gapShades{exp.thisExp}(g+1,:)*ds,[ylims, flip(ylims)], thisCol, 'LineStyle','none','facealpha', 0.2); %Plot pulses in time
    end %Plot pulses in time

     
    % Add cluster statistics for difference from 0
    for s = 1:size(data(1.3*sr:endPoint,:)',1)
        permData.(thisP1){s} = data(1.3*sr:endPoint,s)';
    end
    exp.condDur = size(data(1.3*sr:endPoint,s),1)/sr;
    
    if p1 == 2
        cluster = GT_clusterAlpha(exp,permData);
        
        %save
        clusterStats = repmat(0,1,3000);
        if isfield(cluster.stat1, 'posclusters')
            for i = 1:size(cluster.stat1.posclusters,2)
                if i >0
                    
                    if cluster.stat1.posclusters(i).prob < 0.05
                        clusterStats(find(cluster.stat1.posclusterslabelmat == i)) = 1;
                    end
                end
            end
        end
        if isfield(cluster.stat1, 'negclusters')
            for i = 1:size(cluster.stat1.negclusters,2)
                if i >0
                    if cluster.stat1.negclusters(i).prob < 0.05
                        clusterStats(find(cluster.stat1.negclusterslabelmat == i)) = 1;
                    end
                end
            end
        end
    end
    
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
%     clear data;
    end
end
box on 
xlabel(plt.xlabel); ylabel(['Occipital alpha (dB)']);
xline(0.75*200, '--k')
f.Units = 'centimeters';
f.OuterPosition = [0.25 0.35 11 9];

% Plot cluster stats
thisCol = plt.satPalettes{1,3};
ypos = -2;
x = find(clusterStats(:)); hold on;
plot(x, repmat(ypos,1,length(x)), '*', 'color', thisCol, 'markersize',1.5); hold on;



%Insert topography
axes('position',[.1 .65 .25 .25]);

GT_topoplot(nanmean(dataTopo,2),EEG.chanlocs, ...
    'numcontour', 0, ...
    'plotrad',0.65, ...
    'headrad', 0.5,...
    'emarker2', {occChans,'*','k',3,1}); % To plot location of the selected channel

filename = ['FigS9A_' exp.thisExpCode '_alpha_p1coh.tif'];
saveas(f, [exp.figurepath, '/ga/' filename]);
