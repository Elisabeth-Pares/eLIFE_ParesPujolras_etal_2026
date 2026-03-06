%% GT_Figure1
% Reproduces panels behavioural results in Figure 1

%% BEHAVIOURAL RESULTS
% ===========================================================================
%Load real data
f = figure;
for e = 1:2 %For both experiments
    
    subplot(1,2,e);
    
    exp.thisExp = e; %Experiment to analyse
    [exp, plt, pars]  = GT_setup(exp); %Get correct paths/values
    
    data = readtable([exp.behpath, '/' exp.thisExpCode '_procBehData.csv']); %Generated with R processing script
    
    data = table2array(data(:,[2,11,19,34,7])); %Get only relevant columns (subject, gap, accuracy, condition, npulses)
    
    if e ==2; data(data(:,1) == 1 | data(:,1) == 21,:) = []; end %Exclude bad subjects
    [ga, datsum] = GT_makeSummaryStructure_ave(data) % Retrieve grand-averages
    
    flippedCols = (plt.Rcolors);
    colors = (flippedCols);
    %Plot single-pulse accuracies
    rr1 = plot([-0.3:0.01:-0.2],repmat(squeeze(ga.p(1,1,1)),1,11),  'LineStyle','-','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    rr2 = plot([-0.2:0.01:-0.1],repmat(squeeze(ga.p(1,1,2)),1,11),  'LineStyle', '-','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    hold on;
    
    %Averaged over gaps
    gga{e}.mean = squeeze(nanmean(ga.p,2));
    for c = 1:6
        gga{e}.error(c) = nanstd(squeeze(nanmean(datsum.p(1,:,c,:),2))')./sqrt(size(datsum.p,4));
    end
    
    xpos = [-0.25,-0.15]
    for c = 1:2 %Add error bars
        x = xpos(c);
        y = ga.p(1,1,c);
        error = nanstd(squeeze(datsum.p(1,1,c,:))')./sqrt(size(datsum.p,4));
        errorbar(x,y,error, '.', 'Color', hex2rgb(colors{c*c}), 'MarkerSize', 10); hold on;
    end
    
    ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,3:end)), 'LineWidth', 2, 'LineStyle', '-'); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    %Average trials by P1/P2 coherence
    p1high = [5,6]; p1low = [3,4]; p2high = [3,5]; p2low = [4,6];
    for p = 1:2
        for coh = 1:2
            cohLab = exp.cohLabs{coh};
            idx = ['p' num2str(p) cohLab];
            idx = eval(idx);
            p1p2coh{e}.mean(p,coh,:,:) = squeeze(nanmean(ga.p(1,:,idx),3))
            p1p2coh{e}.error(p,coh,:,:) = nanstd(squeeze(nanmean(datsum.p(1,:,idx,:),3))')/sqrt(size(datsum.p,4));
        end
    end
    
    onePulse_hl{e}.mean(1,:,:) = squeeze(nanmean(ga.p(1,:,2),3));
    onePulse_hl{e}.error(1,:,:) = nanstd(squeeze(nanmean(datsum.p(1,:,2,:),3))')/sqrt(size(datsum.p,4));
    
    onePulse_hl{e}.mean(2,:,:) = squeeze(nanmean(ga.p(1,:,1),3));
    onePulse_hl{e}.error(2,:,:) = nanstd(squeeze(nanmean(datsum.p(1,:,1,:),3))')/sqrt(size(datsum.p,4));
    
    for c = 3:6
        y = ga.p(1,:,c);
        error = nanstd(squeeze(datsum.p(1,:,c,:))')./sqrt(size(datsum.p,4));
        x = exp.gapDur{e}; hold on;
        errorbar(x,y,error, '.', 'Color', hex2rgb(colors{c-2}), 'MarkerSize', 10); hold on;
    end
    
    for k = 1:4
        ll(k).Color = hex2rgb(flippedCols{k});
        hold on;
        
    end
    rr1.Color= hex2rgb(flippedCols{1});
    rr2.Color= hex2rgb(flippedCols{4});
    
    % Compute predicted p(correct) in two-pulse trials based on one-pulse trials - based on Kiani's eq. 4 to 6
    pcorrect = [0.6814,0.7869;0.6861,0.82]; %p(correct) in high/low, based on grand-average (Kiani did it based on log. regression, but said results based on proportions were similar)
    elow = norminv(pcorrect(e,1),0,1);
    ehigh = norminv(pcorrect(e,2),0,1);
    
    colors = flip(flippedCols);
    
    c = 0;
    for p1 = 1:2
        for p2 = 1:2
            c = c+1;
            if p1 == 1; thise1 = ehigh; else thise1 = elow; end
            if p2 == 1; thise2 = ehigh; else thise2 = elow; end
            pcorrect_perf(c) = 1-normcdf(0, thise1+thise2, sqrt(2));
            thisCol = hex2rgb(colors{c});
            xlims = [-0.05:2];
            line('XData',[xlims],'YData', repmat(pcorrect_perf(c), 1, length(xlims)), 'LineStyle', ':','Color', thisCol, 'LineWidth', 2); hold on;
            
        end
    end
    cLabs = {'HH', 'HL', 'LH', 'LL'};
    
    %Labels
    xticks([-0.25,-0.15,exp.gapDur{e}]);
    gapLabels = [];
    for g = 1:length(exp.gapDur{e});
        gapLabels = [gapLabels, {num2str(exp.gapDur{e}(g))}];
    end
    xticklabels({'L-', 'H-', gapLabels{:}})
    set(gca, 'fontsize', 12);
    ylim([0.6,0.95]);
    xlim([-0.4,1.18]);
    xlabel('Gap duration (s)', 'fontsize', 12);
    
    xline(-0.05, '-k')
    box off;
    if e == 1, ylabel(' P(correct)'); end
   
    datsum.p(1,:,:,:) = squeeze(datsum.p(2,:,:,:)); datsum.p(2,:,:,:) = [];
    datsum.n(1,:,:,:) = squeeze(datsum.n(2,:,:,:)); datsum.n(2,:,:,:) = [];
    
    save([exp.behpath, '/GT_ssDat_' exp.thisExpCode], 'datsum');
end

f.Units = 'centimeters';
f.OuterPosition = [2 2 25 10]

exportgraphics(f,[exp.figPath, 'GT_Fig1B.tiff'], 'Resolution', 600)


% %% Additional plots - not in paper
% %% Figure 1D - effects of P1/P2 coherence
% for p = 1:2
%     f = figure;
%     f.Units = 'centimeters';
%     f.OuterPosition = [2 2 8 15]
%     if p == 1; lineType = '-'; else lineType = ':'; end
%     for e = 1:2
%         subplot(2,1,e);
%         [exp, plt, pars]  = GT_setup(exp); %Get correct paths/values
%
%         for coh = 1:2
%             ll = plot(exp.gapDur{e},squeeze(p1p2coh{e}.mean(p,coh,:,:)), 'LineWidth', 2, 'LineStyle', lineType); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
%             x = exp.gapDur{e};
%             y = squeeze(p1p2coh{e}.mean(p,coh,:,:)); %hold on;
%             error = squeeze(p1p2coh{e}.error(p,coh,:,:))
%             errorbar(x,y,error, '.', 'Color',plt.satPalettes{1,coh*coh}, 'MarkerSize', 10); hold on;
%             ll(1).Color = plt.satPalettes{1,coh*coh}; hold on;
%         end
%          xticks([exp.gapDur{e}]);
%     box off;
%     if e == 2; ylabel('                                               P(correct)'); end
%     ylim([0.6,0.95]);
%     xlim([-0.1,1.18]);
%     set(gca, 'fontsize', 12);
%     if e == 2; xlabel('Gap duration (s)', 'fontsize', 12); else xlabel(''); end
%     end
% end
%
% %% Figure 1B - all conditions, averaged over gap
% colors = hex2rgb(flippedCols);
% f = figure;  f.Units = 'centimeters';  f.OuterPosition = [2 2 8 15]
% if p == 1; lineType = '-'; else lineType = ':'; end
% for e = 1:2
%     subplot(2,1,e);
%     [exp, plt, pars]  = GT_setup(exp); %Get correct paths/values
%
%     x = 1:6;
%     y = [gga{e}.mean([1,3,4,2,5,6])];
%     error = [gga{e}.error([1,3,4,2,5,6])];
%
%     colorIndices = [1,1,2,4,3,4];
%     for i = [1:6]
%         b = bar(x(i),y(i)); hold on; errorbar(x,y,error,'.k');
%         b(1).FaceColor = 'flat';
%         b(1).CData = [colors(colorIndices(i),:)];   hold on;
%         if i == 1 || i == 4;     hatchfill2(b,'single','HatchLineWidth', 1.5,'HatchAngle',45,'hatchcolor',[1 1 1]);hold on; end
%     end
%
%     box off;
%     if e == 2; xticks(1:6); xticklabels({'L-', 'LL', 'LH', 'H-', 'HL', 'HH'}); else xticks([]); end
%
%     if e == 2; ylabel('                                               P(correct)'); end
%     ylim([0.6,0.95]);
%     set(gca, 'fontsize', 12);
%     if e == 2; xlabel('Gap duration (s)', 'fontsize', 12); else xlabel(''); end
% end



