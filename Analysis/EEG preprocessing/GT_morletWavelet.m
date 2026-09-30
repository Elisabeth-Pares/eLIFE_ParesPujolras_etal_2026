%% Morlet wavelet transform
% Based on: http://www.mikexcohen.com/lecturelets/manytrials/manytrials.html
% Input:
%   -SUB: subject to be processed
%   -EXPER: structure with experimental parameters for wavelet
%   transformation (OPTIONAL - now everything you need is coded here.)
%   -PAR: parameters used in the files to be decomposed (filter, ica, nan
%   padding, csd...)
%   -PLOTTING: value indicating whether plots are required.Options:
%2- plot 3 sample trials for each individual subject.
% Output:
%   -DAT: cell array with 1 cell per participant, with 3 sub-structures for
%   each:
%       -tf (nfrex x time): average power across trials.
%       -st_tf (nfrex x time x ntrials): single trial power
%       -st_amp (time x ntrials): amplitude is the square root of power. it
%       is calculated as an average across all frequencies.

function [dat] = GT_morletWavelet(sub,exper,par,plt,plotting)

if plotting == 1; figure ; end

%Generate filename to be loaded
if contains(par.epoch,'P1') & ~contains(par.epoch, 'RL')
    filename = (['csd_preproc_P1_' exper.thisExpCode '_P' num2str(sub) '.set']);
elseif contains(par.epoch, 'RL')
    filename = (['csd_preproc_RL_' exper.thisExpCode '_P' num2str(sub) '.set']);
end

EEG = pop_loadset([exper.filepath, '/' filename]);

%Try downsampling for decomposition.
if par.ds == 1
    EEG = pop_resample(EEG, 200);
end

%Frequency decomposition parameters
min_freq = exper.tfa.min_freq;     % Minimum frequency
max_freq = exper.tfa.max_freq;     % Maximum frequency
num_frex = exper.tfa.num_frex;     % Specifies how many frequencies to compute. For a 1Hz resolution, do (max_freq-min_freq).
frex = exper.tfa.frex; % Define vector with centre frequencies to compute - set up in EP1_setup

 
% channel2use = {EEG.chanlocs(exper.movCluster(1,:)).labels, EEG.chanlocs(exper.movCluster(2,:)).labels};%{'B22', 'D19'};
channel2use = channel2use(1:exper.nEEGchans);

% other wavelet parameters
range_cycles = [7];

s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
wavtime = -2:1/EEG.srate:2;
half_wave = (length(wavtime)-1)/2;

% FFT parameters
nWave = length(wavtime);
nData = EEG.pnts * EEG.trials; %Because we are concatenating data 
nConv = nWave + nData - 1; % Number of times the wavelet is dot

% initialize output time-frequency data
tf = zeros(length(frex),EEG.pnts);

tic();
for c = 1:length(channel2use)
    ch = channel2use{c};
    disp(['Decomposing participant P' num2str(sub) ', channel ' ch '...']);
    % Concatenate all trials
    alldata = reshape(EEG.data(strcmpi(ch,{EEG.chanlocs.labels}),1:EEG.pnts,:) ,1,[]);
    %Fast fourier transform data before looping over frequencies for
    %convolution. This won't change over frequencies, so by running it
    %before convoluting we save time! 
    dataX   = fft(alldata, nConv); 
    
    % loop over frequencies
    for fi=1:length(frex)
        
        % create wavelet and get its FFT
        % the wavelet doesn't change on each trial, but it changes on each
        % frequency. Run convolution on concatenated data.
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX ./ max(waveletX); %Normalise 
        
        % now run convolution in one step
        as = ifft(waveletX.*dataX); %AS = analytic signal 
        as = as(half_wave+1:end-half_wave); %
        
        % and reshape back to time X trials
        % NOTE: check that your input dimensions were in this order (points
        % x trials; they should be if you are using eeglab).
        as = reshape(as,EEG.pnts,EEG.trials); %EEG.pnts
        %         disp(num2str(size(as)));
        
        % compute power for single trials
        dat.st_tf(fi,:,1:EEG.trials) = abs(as).^2; %
        
        % compute power and average over trials
        dat.tf(fi,:) = mean(abs(as).^2 ,2);
    end
    
    % compute amplitude for each trial, averaging across frequencies.
    dat.st_amp.alphaBand(:,1:EEG.trials) = squeeze(mean(sqrt(dat.st_tf(exper.tfa.alphaBand,:,1:EEG.trials))));
    dat.st_amp.betaBand(:,1:EEG.trials) = squeeze(mean(sqrt(dat.st_tf(exper.tfa.betaBand,:,1:EEG.trials))));
    dat.st_amp.alphaBetaBands(:,1:EEG.trials) = squeeze(mean(sqrt(dat.st_tf(exper.tfa.alphaBetaBands,:,1:EEG.trials))));
    
    %Save power data in a separate structure
    sumDat.alphaBand(c,:,:) = dat.st_amp.alphaBand;
    sumDat.betaBand(c,:,:) = dat.st_amp.betaBand;
    sumDat.alphaBetaBand(c,:,:) = dat.st_amp.alphaBetaBands;
    
    % save([exper.tfapath '/' exper.thisExpCode '_fullPower_P' num2str(sub) '_' ch '.mat'], 'dat', '-v7.3') %Save datafile with ALL power values for each chanel separately - TOO HEAVY OTHERWISE!
end
toc();

save([exper.tfapath '/Power_DS_' filename(1:end-3) '.mat'], 'sumDat', '-v7.3');


%Plot 3 sample trials (amplitude on the left, power spectrum on the
%right for 3 trials). This is a good way to check that data are
%computed well - the peaks in the amplitude plot should correspond to
%visible high power blobs in the time freq graph.
%%
if plotting == 1 %Plot sample single-trials for each participant
    figure
    for i = [1,3,5]
        subplot(3,2,i)
        plot(dat.st_amp.betaBand(plt.xlims(1,1):plt.xlims(1,2),i+10))
        xlabel('Time (s)'); ylabel('Beta amp');
        yLimits = get(gca,'YLim');  %# Get the range of the y axi
        %         line([500,500], [yLimits], 'color', 'k');
        %         xTick = [100,300,500];
        %         set(gca, 'XLim', [100,600],...
        %             'XTick', xTick,...
        %             'XTickLabel', {'-2', '-1', '0'})
        
        subplot(3,2,i+1)
        contourf(1:plt.xlims(1,2)-plt.xlims(1,1)+1,frex,dat.st_tf(:,plt.xlims(1,1):plt.xlims(1,2),i+10),40,'linecolor','none') %EEG.times(plt.xlims(1,1):plt.xlims(1,2)
        xlabel('Time (s)'); ylabel('Hz');
        yLimits = get(gca,'YLim');  %# Get the range of the y axi
        line([0,0], [yLimits], 'color', 'k' );
        yline([13,30], 'w')
        %         xTick = [-2000,-1000,0];
        %         set(gca, 'XLim', [-2000,500],...
        %             'XTick', xTick,...
        %             'XTickLabel', {'-2', '-1', '0'})
        colorbar
        colormap('jet');
    end
end

end