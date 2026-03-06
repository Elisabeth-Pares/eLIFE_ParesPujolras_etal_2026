GT_fitRegressions 
%% ==============================================================================================================================================
% Single-trial regressions
% CPP_P1/CPP_P2: CPP_p1/p2 ~ p1coh + p2coh (Fig. 3)
% PxAlpha: Choice ~ p1coh + p2coh + p1coh*alpha + p2coh*alpha (Fig. S11)
% PxCPP: Choice ~ p1coh + p2coh + p1coh*CPP + p2coh*CPP (Fig. S11)
%% ==============================================================================================================================================

GA = [];

allP = [];

par.modelName = 'PxCPP'; %'PxCPP';%CPP_PxphAlpha'PxpreStimAlpha_PxAlpha'; %PxAlpha, PxAlphaxDelta, PxpreStimAlpha
par.baseline = 1; %Baseline before P2?


for id = 1:length(exp.sub_id{exp.thisExp})
    disp(['Processing ID ' num2str(id) '...'])
    
    if  strcmp(par.modelName, 'PxAlpha') %For Fig. S11
        [GA, allP] = GT_S11_fitRegression_Alpha(exp,par,id, GA, allP);
    elseif strcmp(par.modelName, 'PxCPP') %For Fig. S11
        [GA, allP] = GT_S11_fitBehRegression_CPP(exp,par,id, GA, allP);
    elseif  strcmp(par.modelName, 'CPP_P1') || strcmp(par.modelName, 'CPP_P2') %For regression results in Fig. 3
        [GA, allP] = GT_Fig3_fitCPP_coh(exp,par,id, GA, allP)
    end
    
end
% if  strcmp(par.modelName, 'PxAlpha') %For Fig. S11
%     save([exp.plotpath, '/' exp.thisExpCode 'AllP_beh_alpha.mat'], 'allP')
% elseif strcmp(par.modelName, 'PxCPP') %For Fig. S11
%     save([exp.plotpath, '/' exp.thisExpCode 'AllP_beh_CPP.mat'], 'allP')
save([exp.plotpath, '/' exp.thisExpCode '_' par.modelName '_behRegression_b' num2str(par.baseline) '.mat' ], 'GA')



