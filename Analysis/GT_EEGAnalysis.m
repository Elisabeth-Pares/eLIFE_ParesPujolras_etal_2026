%% Gaps Task EEG analysis
%%--------------------------------------------------------------------------------------------
% Analysis pipeline used in: 

% Parés-Pujolràs, E., Geuzebroek, A. C., O’Connell, R. G., & Kelly, S. P.
% (2025). Perceptual glimpses are locally accumulated and globally 
% maintained at distinct processing levels. eLife, 14.

% EEG analysis code for Experiments 1 & 2 
% Author: Elisabeth Parés-Pujolràs
% Date: 02/02/2025
%%--------------------------------------------------------------------------------------------

%% Set up experimental analysis parameters
clc
clear all

exp.thisExp = 2;        %Experiment to analyse: set to 1 or 2
[exp, plt]  = GT_setup(exp);

%Participant IDs for each experiment
exp.sub_id{1}            = [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]; % P8 excluded because only 4 blocks were recorded.
exp.sub_id{2}            = [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23];    % P1 & P21 excluded because they could not perform task well

% ==============================================================================================================================================
% PREPROCESSING 
% ==============================================================================================================================================
%% Convert data & add trial information to raw EEG data files 
% The output of this step are the raw data shared on the OSF repository
for sub = exp.sub_id{exp.thisExp}(1:end) %For each subject 
    GT_convertData(sub, exp)
end
%% Preprocess data (Part 1) - DONE
for sub = exp.sub_id{exp.thisExp}(1:end)%For each subject
    GT_preprocess_preICA(sub,exp)   
end
%% Manually remove eye movements using ICA 
%% Preprocess (Part 2) - DONE
for sub =  exp.sub_id{exp.thisExp}(1:end) %For each subject
    GT_preprocess_postICA(sub,exp)
end

%% CSD filter preprocessed data 
for sub = exp.sub_id{exp.thisExp}(1:end) %For each subject
    for epoch = 1:2
        GT_applyCSD(sub,exp, epoch)
    end
end

% ==============================================================================================================================================
% TF DECOMPOSITION 
% ==============================================================================================================================================
%% Morlet wavelet transform preprocessed data 
par.ds = 1; %Downsample data to 200Hz
for sub = exp.sub_id{exp.thisExp}(1:end) %For each subject
    for epoch = 1
        par.epoch = exp.epochs{epoch}  %Options: P1 (pulse-1 aligned), RL (response-locked)
        GT_morletWavelet(sub,exp,par,plt,0)
    end
end

% ==============================================================================================================================================
%% Average by condition
% ==============================================================================================================================================
par.average = 1;

allP = []; cppAll = [];
w = warning ('off','all');

for dataType = {'TFA'}; %'ERP' 'TFA',
    par.dataType = dataType;
    par.band = 'beta';
    for epoch = 1
        par.epoch = exp.epochs{epoch};
        id = 0;
        for sub = exp.sub_id{exp.thisExp}(1:end)
            disp(['Processing P' num2str(sub) '...'])
            id = id+1;
            if strcmp(par.dataType, 'ERP')
                %Tags specify which factors we group trials by.
                %NP = number of pulses; GAP = gap duration, etc.
                par.groupings = {'np_gap_', 'np_gap_p1coh_',...
                    'np_gap_p1coh_', 'np_gap_p1coh_acc_', ...
                    'np_gap_p1p2coh_', 'np_gap_p1p2coh_acc_','np_gap_dir_acc_',...
                    'np_gap_dir_acc_'}
            else
                par.groupings = {'np_gap_dir_acc_', 'np_gap_p1coh_dir_', 'np_gap_p1coh_dir_acc_', ...
                    'np_gap_p1p2coh_dir_', 'np_gap_p1p2coh_dir_acc_',...
                    'np_gap_dir_p1coh_'}; 
            end
            tic()
            [allP, cppAll] = GT_averageBygroup_tidy(id,exp, par, allP, cppAll)
            toc()
        end
        if strcmp(par.dataType, 'ERP')
            %This file contains all condition sorting, with all channels
            save([exp.plotpath '/' exp.thisExpCode, '_' par.epoch '_allP_' dataType{1}, '_averageData.mat'], 'allP', '-v7.3')
        else
            save([exp.plotpath '/' exp.thisExpCode, '_' par.epoch '_allP_' dataType{1}, '_averageData' par.band '.mat'], 'allP', '-v7.3')
        end

    end
end
















