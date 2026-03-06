% Run cluster-based permutation tests
% Permutation parameters 
function [cluster] = GT_clusterAlpha(exp,permData)

cfg_dat = struct;

chans = [1]; %Motor beta lateralisation
load('chanlocsBioSemi_128_ext8.mat') %Load channel locations

cfg_dat.label       = {chanlocs.labels};  %Channel locations; (exp.CPPCluster_small)

cfg = struct();
cfg.feedback = 'no';

cfg_dat.fsample = 200;  % sampling rate
cfg_dat.time = [-0.2:1/cfg_dat.fsample:exp.condDur-0.205]; 
cfg_dat.dimord = 'chan_time';

cfg.minnbchan   = [0];
cfg.latency     = 'all'; %Analysis period
cfg.channel     = 1;%%exp.centroParietal;%parietalFocus%exp.centroParietal;  %
cfg.avgoverchan = 'yes';
cfg.avgovertime = 'no'; %no
cfg.parameter   = 'values';
cfg.method      = 'montecarlo';
cfg.clusterthreshold    = 'parametric'%'nonparametric_common'; %EPP added

cfg.statistic   = 'depsamplesT';  
cfg.correctm    = 'cluster';
cfg.clusteralpha = 0.05;  % intial threhsold for forming clusters
cfg.clusterstatistic = 'maxsum';  % method for quantifying combined cluster test statistics

cfg.alpha       = 0.05;
cfg.tail        = 0; % two-sided test
cfg.correcttail = 'prob'; 

cfg.numrandomization = 1000;

%Design matrix 
Nsub = length(exp.sub_id{exp.thisExp}(1:end))
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

%Null condition 
for s = 1:length(exp.sub_id{exp.thisExp}(1:end))
    all_trl_high{s} = cfg_dat; all_trl_high{s}.values = permData.high{s}(chans,:);
    all_trl_low{s} = cfg_dat; all_trl_low{s}.values = permData.low{s}(chans,:);

end

stat1 = ft_timelockstatistics(cfg,all_trl_high{:},all_trl_low{:});

cluster.stat1 = stat1;