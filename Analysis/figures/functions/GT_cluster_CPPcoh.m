% Run cluster-based permutation tests
% Permutation parameters 
function [cluster] = GT_cluster_CPPcoh(exp,par,permData)

cfg_dat = struct;
chans = [1]; %Average of all CPP channels

cfg_dat.label       = {'CPP'}

cfg = struct();
cfg.feedback = 'no';

cfg_dat.fsample = 512;  % sampling rate
cfg_dat.time = [-0.2:1/cfg_dat.fsample:1]; 
cfg_dat.dimord = 'chan_time';

cfg.minnbchan   = [0];
cfg.latency     = par.testLatency;%'all'; %Analysis period
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

 
cfg.numrandomization = 10000;

%Design matrix 
Nsub = length(exp.sub_id{exp.thisExp}(1:end))
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

%Null condition 
for s = 1:length(exp.sub_id{exp.thisExp}(1:end))
    all_trl_p1{s} = cfg_dat; all_trl_p1{s}.values = permData.p1{s}(chans,:);
    all_trl_p2{s} = cfg_dat; all_trl_p2{s}.values = permData.p2{s}(chans,:);

    all_trl_null{s} = cfg_dat;
    all_trl_null{s}.values = zeros(size(all_trl_p1{1}.values,1),size(all_trl_p1{1}.values,2))
        
end

stat1 = ft_timelockstatistics(cfg,all_trl_p1{:},all_trl_null{:});
stat2 = ft_timelockstatistics(cfg,all_trl_p2{:},all_trl_null{:});

cluster.stat1 = stat1;
cluster.stat2 = stat2;