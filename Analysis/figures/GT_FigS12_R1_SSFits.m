%% Figure 4
%Load real data
% f = figure;
for e = 1:2 %For both experiments
    f = figure; t = tiledlayout(5,5);
    
    exp.thisExp = e; %Experiment to analyse
    [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values
    
    %Load ss data
    load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\SingleSubjects\GT_ssDat_' exp.code{e} '.mat']);
    
    for s = 1:length(exp.sub_id{e})
%         nexttile;
        ga = [];
        ga.p = datsum.p(:,:,:,s);
        ga.n = datsum.n(:,:,:,s);
        %Load this subjects simulation
        load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_S_ModelSims_ss\FixedBound_1DV_2DR_E' num2str(e) '_P' num2str(s) '_fitType0_N10000'])
        sims = allSims;
        sims.simdat(sims.simdat(:,2) == -1,2) = 0;
        
        pars = allPars;
        pars.accType = 'bin';
        
        acc = []; p1acc = []; accLab = {}; noBacc = [];
        
        for c =length(pars.cLab):-1:3
            thisCond = sims.simdat(:,1) == c;
            
            if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
            acc(c-2) = mean(sims.simdat(thisCond,thisAcc));
            accLab{c-2} = pars.cLab{c};
            %         acc_perf(c-2) = mean(sims.simdat(thisCond, 5)); %Perfect accumulation
        end
        
        acc = reshape(acc, [4,length(exp.gapDur{e})]);
        accLab = reshape(accLab, [4, length(exp.gapDur{e})]);
        box on
        
        for c = 1:2%length(pars.cLab)
            thisCond = sims.simdat(:,1) == c;
            if strcmp(pars.accType, 'bin'), thisAcc = 2; elseif strcmp(pars.accType, 'erf'), thisAcc = 3; end
            sims.simdat(sims.simdat(:,2) == -1,2) = 0;
            p1acc(c) = mean(sims.simdat(thisCond,thisAcc));
            p1accLab{c} = pars.cLab{c};
        end
        
        acc = [repmat(p1acc(2),1,length(exp.gapDur{e}));repmat(p1acc(1),1,length(exp.gapDur{e})); flip(acc)];
        
        % %         subplot(2,4,[1,2]+4*(e-1))
        %
        %         r1 = plot([-0.3:0.01:-0.2],repmat(acc(1,1)',1,11), 'LineStyle', '-','LineWidth', 2); hold on;
        %         rr1 = plot([-0.3:0.01:-0.2],repmat(squeeze(ga.p(1,1,1)),1,11),  'LineStyle',':','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
        %
        %         r2 = plot([-0.2:0.01:-0.1],repmat(acc(2,1)',1,11), 'LineStyle', '-','LineWidth', 2); hold on;
        %         rr2 = plot([-0.2:0.01:-0.1],repmat(squeeze(ga.p(1,1,2)),1,11),  'LineStyle', ':','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
        %
        allP_acc{e}(s,:,:) = acc;
        %
        %         hold on;
        %
        %         l = plot(exp.gapDur{e},acc(3:end,:)', 'LineStyle', '-','LineWidth', 2); hold on;
        %         ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,3:end)), 'LineWidth', 2, 'LineStyle', ':'); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
        %
        %         flippedCols = (plt.Rcolors);
        %
        %         for k = 1:4
        %             l(k).Color = hex2rgb(flippedCols{k});
        %             ll(k).Color = hex2rgb(flippedCols{k});
        %             hold on;
        %         end
        %         r1.Color = hex2rgb(flippedCols{1});
        %         rr1.Color= hex2rgb(flippedCols{1});
        %
        %         r2.Color = hex2rgb(flippedCols{4});
        %         rr2.Color= hex2rgb(flippedCols{4});
        %
        %         xticks([-0.25,-0.15,exp.gapDur{e}]);
        %         gapLabels = [];
        %         for g = 1:length(exp.gapDur{e});
        %             gapLabels = [gapLabels, {num2str(exp.gapDur{e}(g))}];
        %         end
        %         xticklabels({'L-', 'H-', gapLabels{:}})
        %         set(gca, 'fontsize', 11);
        %         ylim([0.6,0.9]);
        %         xlim([-0.4,1.08]);
        %         if e == 2; xlabel('Gap duration (s)', 'fontsize', 11); else xlabel(''); end
        %         %     ylabel('P(correct)', 'fontsize', 11);
        %
        %         xline(-0.05, '-k')
        %         box on;
        %         if e == 2; ylabel('                                               P(correct)'); end
        %
        %
        %% Early terminations
         for conds = 3:6%length(pars.cLab)/4
        thisLab = pars.cLab{conds}(6:end);
        indices{conds} = contains(pars.cLab, thisLab);
        condList{conds-2} = thisLab;
         end
        
        for c = 1:length(pars.cLab)
            thisP1prop_z1{e}(s,c) = nanmean(sims.bndTP1(c,:));
        end
        
        for conds = 3:6
            thisP1prop_z1xgap{e}(s,conds) = nanmean(thisP1prop_z1{e}(s,indices{conds}),2);
        end
        
        %% Plot DV at P1
        for c = 1:length(pars.cLab)
            idx = sims.bndTP1(c,:) == 0;
            thisP1_DV{e}(s,c) = nanmean(abs(sims.p1DV_end(c,idx)),2);
        end
        
        for conds = 3:6
            thisP1_DV_xgap{e}(s,conds) = nanmean(thisP1_DV{e}(s,indices{conds}),2);
        end
    end
end

%% GA DATA
acc = []
f= figure;
%     f = figure; t = tiledlayout(5,5);

for e = 1:2 %For both experiments
    acc_y = squeeze(nanmean(allP_acc{e}));
    acc_err = squeeze(nanstd(allP_acc{e})./sqrt(size(allP_acc{e},1)));
    
    exp.thisExp = e; %Experiment to analyse
    [exp,  pars]  = GT_S_setup(exp); %Get correct paths/values
    
    %Load ss data
    load(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode '.mat']);
    subplot(2,4,[1,2]+4*(e-1))
    r1 = plot([-0.3:0.01:-0.2],repmat(acc_y(1,1)',1,11), 'LineStyle', '-','LineWidth', 2); hold on;
    errorbar(-0.25, acc_y(1,1),acc_err(1,1), '.k'); hold on;
    
    rr1 = plot([-0.3:0.01:-0.2],repmat(squeeze(ga.p(1,1,1)),1,11),  'LineStyle',':','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    r2 = plot([-0.2:0.01:-0.1],repmat(acc_y(2,1)',1,11), 'LineStyle', '-','LineWidth', 2); hold on;
    errorbar(-0.15, acc_y(2,1),acc_err(2,1), '.k'); hold on;
    
    rr2 = plot([-0.2:0.01:-0.1],repmat(squeeze(ga.p(1,1,2)),1,11),  'LineStyle', ':','LineWidth', 2); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    
    hold on;
    
    l = plot(exp.gapDur{e},acc_y(3:end,:)', 'LineStyle', '-','LineWidth', 2); hold on;
    hold on;
    errorbar(repmat(exp.gapDur{e},4,1), acc_y(3:end,:),acc_err(3:end,:), '.k'); hold on;
    ll = plot(exp.gapDur{e},squeeze(ga.p(1,:,3:end)), 'LineWidth', 2, 'LineStyle', ':'); hold on; %', 'none', 'Marker', '*', 'MarkerSize', 5); hold on;
    
    flippedCols = (plt.Rcolors);
    
    for k = 1:4
        l(k).Color = hex2rgb(flippedCols{k});
        ll(k).Color = hex2rgb(flippedCols{k});
        hold on;
    end
    r1.Color = hex2rgb(flippedCols{1});
    rr1.Color= hex2rgb(flippedCols{1});
    
    r2.Color = hex2rgb(flippedCols{4});
    rr2.Color= hex2rgb(flippedCols{4});
    
    xticks([-0.25,-0.15,exp.gapDur{e}]);
    gapLabels = [];
    for g = 1:length(exp.gapDur{e});
        gapLabels = [gapLabels, {num2str(exp.gapDur{e}(g))}];
    end
    xticklabels({'L-', 'H-', gapLabels{:}})
    set(gca, 'fontsize', 11);
    ylim([0.6,0.9]);
    xlim([-0.4,1.08]);
    if e == 2; xlabel('Gap duration (s)', 'fontsize', 11); else xlabel(''); end
    %     ylabel('P(correct)', 'fontsize', 11);
    
    xline(-0.05, '-k')
    box on;
    if e == 2; ylabel('                                               P(correct)'); end
    
    %% Early terminations
    subplot(2,4,[3]+4*(e-1))
    % Panel C, %early terminations
    flippedCols = flip(plt.Rcolors);
    i = 0;
    for c = [6:-1:3]
        i = i+1;
        y = mean(thisP1prop_z1xgap{e}(:,c)); y_err = nanstd(thisP1prop_z1xgap{e}(:,c)./sqrt(size(thisP1prop_z1xgap{e}(:,c),1))); 
        l = bar(i,y); hold on; errorbar(i,y,y_err, '.k'); hold on; 
        l.FaceColor = hex2rgb(flippedCols{c-2});
        yline(0.5, '-k')
        %         ylabel({exp.name{e}});
        if e == 2; ylabel('                                               P(Early termination)'); end
        xticks([1,2,3,4]); xticklabels({'LL','LH','HL', 'HH'});
        
        ylim([0.2,0.9])
        set(gca, 'fontsize', 11);
        
        box on;
        
    end
    
    %% Plot DV at P1
    subplot(2,4,[4]+4*(e-1))
 
    flippedCols = flip(plt.Rcolors);
    i = 0;
    for c = [6:-1:3]
        i = i+1;
        y = mean(thisP1_DV_xgap{e}(:,c),1); y_err = std(thisP1_DV_xgap{e}(:,c))./sqrt(size(thisP1_DV_xgap{e},1));
        l = bar(i,y); hold on; errorbar(i,y,y_err, '.k');
        l.FaceColor = hex2rgb(flippedCols{c-2});
        %         yline(0.5, '--k')
        %         ylabel({exp.name{e}});
        if e == 2; ylabel('                                         |DV| at P1 if B not hit'); end
        xticks([1,2,3,4]); xticklabels({'LL','LH','HL', 'HH'});
        
                if e == 2; ylim([0.013,0.016]); elseif e == 1; ylim([0.015,0.02]); end
        set(gca, 'fontsize', 11);
        ylim([0.014,0.025]);
        box on;
        
    end
    
end
f.Units = 'centimeters';
f.OuterPosition = [2 2 21 14]
exportgraphics(f, [exp.figPath, 'GT_FigS13.tiff'], 'Resolution', 600)
