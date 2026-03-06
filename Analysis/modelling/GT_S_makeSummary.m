%Get grand-averaged data for each experiment so SONIC can directly load
%.mat files 
%Script does NOT need to be on SONIC

figure; 

for e = 1:2 %For both experiments
    subplot(1,2,e); 
    exp.thisExp = e; %Experiment to analyse
    [exp, plt, pars]  = GT_setup(exp); %Get correct paths/values
    
    
    data = readtable([exp.behpath, '/' exp.thisExpCode '_procBehData.csv']); %Generated with R processing script
    data = table2array(data(:,[2,11,19,34,7])); %Get only relevant columns (subject, gap, accuracy, condition, npulses)
    
    if e ==2; data(data(:,1) == 1 | data(:,1) == 21,:) = []; end
    [ga, datsum] = GT_makeSummaryStructure_ave(data) % Fit to grand-averages for simplicity
    
    colors = (plt.Rcolors);
    for c = 3:6;
        x = exp.gapDur{e};
        plot(x, squeeze(ga.p(1,:,c)), 'Color', hex2rgb(colors{c-2}), 'LineWidth', 2);
        hold on; ylabel('P(correct)'); xlabel('Gap dur (s)');
    end
    hold on; 
    plot(x, squeeze(ga.p(1,:,1)), 'Color', hex2rgb(colors{1}), 'LineWidth', 2, 'LineStyle', '--'); hold on; 
    plot(x, squeeze(ga.p(1,:,2)), 'Color', hex2rgb(colors{4}), 'LineWidth', 2, 'LineStyle', '--');
    ylim([0.65,0.9])
    title(exp.thisExpName)
    save(['C:\Users\elisa\Desktop\GapsTask\Modelling & simulations\GT_Sonic\GT_gaDat_' exp.thisExpCode], 'ga')
end
