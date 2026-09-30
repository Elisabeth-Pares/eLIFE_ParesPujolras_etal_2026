function [sims, pars] = GT_S_simul_DDM_clean(pars, pvar, pfix, Sel)

pm = zeros(1,length(Sel)); % Parameters plus N and the seed
pm(Sel==1) = pvar;  %Parameters to fit
pm(Sel==0) = pfix;  %Fixed parameters

% Pull the parameter values out of the input parameter vector:
[parb, bc, dlow, dhigh,dp2, accT, accTvar, eta, sz, sd, N, seed] = deal(pm(1),pm(2),pm(3),pm(4),pm(5),pm(6),pm(7),pm(8), pm(9), pm(10), pm(11), pm(12));

%% Simulation parameters
s = 0.1;            % within-trial noise - scaling parameter
dt = 1/pars.sr;     % time increment based on sampling rate of simulations
t = [-.75:dt:3];    % Full time range to simulate.

wfDV = nan(length(dlow),N,length(t));   %waveform DV
bndTDist = nan(length(dlow), N);        %bound distribution

N = pars.nSim;      %Number of simulations
pulseDur = 0.201;   %All pulses of evidence last 200 ms

bndTP1 = []; bndTP2 =[];

%Explored variants of the model fitted on subsets of conditions; to fit to
%all conditions, set pars.simpleFit = 0. 
if pars.simpleFit == 0, condIdx = 1:length(pars.cLab);
elseif  pars.simpleFit== 1; condIdx = 1:6;
elseif pars.simpleFit == 2; condIdx = [1:6, length(pars.cLab)-3:length(pars.cLab)];
elseif pars.simpleFit == 3; condIdx = [1:2, length(pars.cLab)-3:length(pars.cLab)];
end

for c=condIdx(1:end)      
    
    if dhigh ~= 0, %If estimating two separate drift rates for high/low coherence
        ev1 = 1; ev2 = 1; 
    end
    
    % Initialise
    ch = zeros(N,1);
    pars.cLab{c};
    if contains(pars.cLab{c}, 'onePulse');
        isOnePulse = 1; else; isOnePulse = 0; end
    
    % Get times of evidence for accumulator
    thisGap = extractBefore(pars.cLab{c}, '_P1');
    gapDur = pars.gDur(c);
    p2Time = [gapDur + pulseDur];
    
    ch = repmat(nan, N, 1); ch_1p = repmat(nan, N, 1); ch_noB = repmat(nan, N, 1);
    for n=1:N % each single trial
        
        if dhigh == 0; d1 = dlow; d2 = dlow; %If estimating one single drift rate, assume d1 = d2
        else, %Use the name of conditions to set the appropriate drift rate for each pulse
            if contains(pars.cLab{c}, 'P1low') %Low drift rate
                d1 = dlow;      %Set drift rate for pulse 1 (d1) as dlow (free parameter for drift rate 1 = dlow, low coh dr)
            else d1 = dhigh;    %Set drift rate for pulse 1 (d1) as dhigh (free parameter for drift rate 2 = dhigh, high coh dr)
            end
            if contains(pars.drModulation, 'invDRP2')
                thisdp2 = 1/d1*dp2; %Inverse scaling with P1 coh
            else
                thisdp2 = dp2;
            end
            
            if contains(pars.cLab{c}, 'P2low')
                d2 = dlow*thisdp2;      %Set drift rate for pulse 2 (d2) as dlow (free parameter for drift rate 1 = dlow, low coh dr)
            else d2 = dhigh*thisdp2;    %Set drift rate for pulse 2 (d2) as dhigh (free parameter for drift rate 2 = dhigh, high coh dr)
            end
        end
        
        DVtt = [];                      %Initialise DV for this trial (DVtt)
        accDur = pulseDur*2 + gapDur;   %Set duration of accumulation process to the duration of both pulses, plus the gap duration
        
        % simulate trial - 'tt' stands for 'this trial'
        dtt_p1 = d1*ev1 + randn*eta;         % Simulate P1 ev
        dtt_p2 = d2*ev2 + randn*eta;         % Simulate P2 ev
        
        stt = (rand-0.5)*sz;              % starting point for this trial. If sz = 0, stt = 0.
        accTtt = accT+(rand-0.5)*accTvar; % accumulation onset for THIS trial. If accTvar == 0, start at time = 0
        accOffTtt = accTtt+accDur+(rand-0.5)*accTvar; % accumulation offset for this trial
        
        decayDelay = pars.decayDelay;     %How long does it take for the signal to start falling down after dots turn blue/a bound is hit?
        
        %Simulate each period of time in the trial: 
        % leadIn dots, p1, gap, p2, leadOut dots
        % leadIn, gap and leadOut are assumed to = 0, because accumulation
        % is not allowed to happen in these simulations anyway.
        
        leadIn = zeros(1,length(find(t<=0)))*s*sqrt(dt);

        p1Ev = dtt_p1*dt+randn(1,length(find(t>0 & t<=pulseDur)))*s*sqrt(dt);
        gapEv = zeros(1,length(find(t>pulseDur & t<p2Time)))*s*sqrt(dt);
        if ~contains(pars.cLab{c}, 'onePulse')
            p2Ev = dtt_p2*dt+randn(1,length(find(t>=p2Time & t<=p2Time+pulseDur)))*s*sqrt(dt);
        else
            p2Ev = [];
        end
        leadOut = zeros(1,length(find(t>p2Time+pulseDur & t <= t(end))))*s*sqrt(dt);
        noisyEv = [leadIn, p1Ev, gapEv, p2Ev, leadOut]; % This never needs to change - it is the physical reality of what was shown
        
        p1ERP = zeros(1,length(noisyEv)); p2ERP = zeros(1,length(noisyEv));
        %Initialise some output values
        p1DV = 0; p1DV_end = 0; p1DV_perf = 0; p2DV_perf = 0;
        
        % Now compute the accumulator waveform for this trial:
        DVtt = [leadIn(1:end-1),cumsum(noisyEv(find(t>= accT & t <= accOffTtt))),...
            ones(1,length(find(t>accOffTtt)))*sum(noisyEv(find(t>accT&t<=accOffTtt)))];
        newDVtt = DVtt; %Create copy of the DV; this corresponds to the CPP simulation, going down when 1) a bound is hit, or 2) dots turn blue. 
        
        %Make bound a vector instead of a point, so we can create
        %collapsing bound after a threshold is crossed
        b = repmat(parb, 1, length(DVtt));
        
        %To simulate MBL, take DVtt and set to Bound level after 1st bound crossing
        MBL = DVtt;
        bCross = find(MBL > b | MBL < -b, 1, 'first'); 
        if ~isempty(bCross)
            MBL(bCross:end) = MBL(bCross); 
        end
        clear bCross; 
        
        %Compute activities exclusively associated with P1 and P2
        p1Idx = (length(leadIn):length(leadIn)+length(p1Ev)+length(gapEv)-1);
        p2Idx = (length(leadIn)+length(p1Ev)+length(gapEv)):(length(leadIn)+length(p1Ev)+length(gapEv)+length(p2Ev)-1);             
        
        p1ERP(p1Idx) = [p1Ev, gapEv]; p1ERP = cumsum(p1ERP); decayP1 = p1ERP; 
        p2ERP(p2Idx) = p2Ev; p2ERP = cumsum(p2ERP); decayP2 = p2ERP; 
        
        %% Determine error/correct output of a given trial
        for p = 1:2
            if isOnePulse; p2Idx = length(p1Ev); end
            
            if p == 2 & isOnePulse || (p == 2 & p1DV ~= 0) % Added post || condition
                p2DV = 0;
                p2ERP = repmat(0,1,length(p1ERP)); decayP2 = p2ERP;
                break;
            end
            if pars.biased == 1 & p == 1 %[In the paper]
                %One single decision process that continues building up throughout
                %trial time, DV1 biases starting point of DV2, and final
                %choice = DV1 + DV2
                thisDV = DVtt;
            elseif pars.biased == 0 %[As in Kiani et al. 2013]
                %Two decision processes, final choice = DV1 + DV2
                thisDV = eval(['p' num2str(p) 'ERP']);
            end
            
            if pars.accGap == 1, upto = p2Idx(1);
            else upto = p1Idx(length(p1Ev)); end
            
            %Simulate perfect integration 
            if p == 1; p1cumev = cumsum(p1Ev); p1DV_perf = p1cumev(end);
            elseif p == 2; p2cumev = cumsum(p2Ev); p2DV_perf = p2cumev(end);end
        
            if p == 1 & ~isinf(b(1))                            %Bound set to inf in unbounded models
                bndT = find(abs(thisDV(1:upto))> b(1:upto),1);  % find the sample point where it crosses bound b
            elseif p == 2 & ~isinf(b(1))
                bndT = find(abs(thisDV(p2Idx(1):end))> b(p2Idx(1):end), 1); % find the sample point where it crosses bound b
                if ~isempty(bndT);
                    bndT = bndT + p2Idx(1)-1; end
            elseif isinf(b(1))
                bndT = [];
            end
            
            if ~isempty(bndT)
                if p == 1;
                    p1DV = thisDV(bndT); p1DVbias = thisDV(bndT);
                    p1DV_end = thisDV(bndT);
                elseif p == 2;
                    p2DV = thisDV(bndT);
                end
                
                bndTDist(c,n) = bndT;
                if pars.biased == 1;
                    if p == 1;
                        bndTP1(c,n) = 1;
                        bndTP2(c,n) = 0;
                    else
                        bndTP1(c,n) = 0;
                        bndTP2(c,n) = 1;
                    end
                else
                    if p == 1;
                        bndTP1(c,n) = 1;
                        bndTP2(c,n) = 0;
                    else
                        bndTP2(c,n) = 1;
                        bndTP1(c,n) = 0;
                        
                    end
                end
                
                pkDV = DVtt(bndT+decayDelay); % the amplitude of this "peak" where it has just crossed the bound
                pk_decayDV = newDVtt(bndT+decayDelay);
                pkDV_ERP = thisDV(bndT+decayDelay); % the amplitude of this "peak" where it has just crossed the bound
                pkP1 = decayP1(bndT+decayDelay);
                
                if (p == 1 & bndT <= upto) || (p == 2 & bndT >= p2Idx(1)) % p == 1& bndT <= p2Idx)
                    if decayDelay ~= 0; bndT = bndT + decayDelay; end %Test: introduce fixed delay after bound crossing before signal falls down

                    if p == 1
                        %Linear decay after hitting bound
                        DVtt(bndT+ [1:pars.stepAfterbndT]) = pkDV-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*pkDV;
                        newDVtt(bndT+ [1:pars.stepAfterbndT]) = pkDV-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*pkDV;
                        p1ERP(bndT+[1:pars.stepAfterbndT])=pkDV_ERP-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*pkDV_ERP;
                        decayP1(bndT+[1:pars.stepAfterbndT])=pkP1-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*pkP1;

                        DVtt((bndT+pars.stepAfterbndT+1):p2Idx(1))=0; 
                        newDVtt((bndT+pars.stepAfterbndT+1):p2Idx(1))=0; 
                        p1ERP((bndT+pars.stepAfterbndT+1):end)=0;
                        decayP1((bndT+pars.stepAfterbndT+1):end)=0;

                        if pars.biased == 0
                            DVtt(bndT+pars.stepAfterbndT+1:end) =0; 
                            newDVtt(bndT+pars.stepAfterbndT+1:end)=0; 
                        else
                            DVtt(bndT+pars.stepAfterbndT:end) = 0;
                            newDVtt(bndT+pars.stepAfterbndT+1:end)=0;
                        end
                        
                        b = b(1:length(DVtt));
                        DVtt = DVtt(1:length(t)); %cut if extra length added
                        newDVtt = newDVtt(1:length(t));

                        if pars.biased == 1 & p == 1
                            %One single decision process that continues building up throughout
                            %trial time, DV1 biases starting point of DV2, and final
                            %choice = DV1 + DV2
                            thisDV = DVtt;
                        elseif pars.biased == 0
                            %Two decision processes, final choice = DV1 + DV2
                            thisDV = eval(['p' num2str(p) 'ERP']);
                        end
                        
                    elseif p == 2
                        %Fast decay after hitting bound
                        pkDV_ERP = p2ERP(bndT); % the amplitude of this "peak" where it has just crossed the bound
                        
                        DVtt(bndT+[1:pars.stepAfterbndT]) = pkDV-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*pkDV;
                        newDVtt(bndT+[1:pars.stepAfterbndT]) = pk_decayDV-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*pk_decayDV;

                        p2ERP(bndT+[1:pars.stepAfterbndT])=(pkDV_ERP)-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*(pkDV_ERP-p2ERP(1));
                        decayP2(bndT+[1:pars.stepAfterbndT])=(pkDV_ERP)-[1:pars.stepAfterbndT]/(pars.stepAfterbndT)*(pkDV_ERP-p2ERP(1));

                        DVtt((bndT+pars.stepAfterbndT+1):end)=0; 
                        newDVtt((bndT+pars.stepAfterbndT+1):end)=0; 
                        decayP2((bndT+pars.stepAfterbndT+1):end)=0;

                        p2ERP((bndT+pars.stepAfterbndT+1):end)=p2ERP(bndT+pars.stepAfterbndT);

                        if pars.biased == 1
                            %One single decision process that continues building up throughout
                            %trial time, DV1 biases starting point of DV2, and final
                            %choice = DV1 + DV2
                            thisDV = DVtt;
                        elseif pars.biased == 0
                            %Two decision processes, final choice = DV1 + DV2
                            thisDV = eval(['p' num2str(p) 'ERP']);
                        end
                        
                        if pars.biased == 0 || p1DV == 0%If not biased OR no bound hit at all
                            p2DV = thisDV(bndT)  ; %0
                        elseif pars.biased == 1 & p1DV ~= 0 %If no bound hit, do not compute
                            p2DV = 0;
                            p2ERP = repmat(0, length(p2ERP), 1)'; %NO accum if bound already hit
                            decayP2 = repmat(0, length(p2ERP), 1)';
                        end
                        
                    end
                    
                end
            else
                if p == 1;
                    bndTP1(c,n) = 0;
                else
                    bndTP2(c,n) = 0;
                end
                
                if p == 1;
                    if pars.biased == 0 || isOnePulse
                        p1DV = thisDV(upto); p1DVbias = thisDV(upto);
                        p1DV_end = thisDV(upto);

                    elseif pars.biased == 1 
                        p1DV = 0; p1DVbias = thisDV(upto);
                        p1DV_end = thisDV(upto);
                    end
                    
      
                    if ~isOnePulse
                        p2ERP = zeros(1,length(noisyEv));
                        p2ERP(p2Idx) = p2Ev; 
                        p2ERP = cumsum(p2ERP);
                        decayP2(p2Idx) = p2ERP(p2Idx); 
                        
                        % Now compute the accumulator waveform for this trial:
                        DVtt = [leadIn(1:end-1),cumsum(noisyEv(find(t>= accT & t <= accOffTtt))),...
                            ones(1,length(find(t>accOffTtt)))*sum(noisyEv(find(t>accT&t<=accOffTtt)))];% for convenience, let's assume that when the accumulator "switches off" at a certain time, it stays at that level - this makes it easy to use the 'end' value below for the choice, when there is no bound
                        newDVtt = DVtt;
                        
                        %% Bring signal down for between end of P1 and onset of P2; does not apply to 0 gap trials
                        if gapDur ~= 0
                            decayEnd= min(upto+pars.stepAfterOffset, p2Idx(1));
                            temp = DVtt; 
                            temp(upto+[1:pars.stepAfterOffset]) = DVtt(upto+1)-[1:pars.stepAfterOffset]/(pars.stepAfterOffset)*DVtt(upto+1);%+randn(1,pars.stepAfterbndT)*s*sqrt(dt);
                            newDVtt(upto:decayEnd) = temp(upto:decayEnd);
                            if decayEnd<p2Idx(1)
                                newDVtt(decayEnd+1:p2Idx(1)-1) = 0;
                                 offset = max(decayEnd, p2Idx(1)-1);
                                 newDVtt(p2Idx(1):end) = DVtt(p2Idx(1):end)-DVtt(offset);
                            else
                                offset = DVtt(p2Idx(1)) - newDVtt(p2Idx(1));
                                newDVtt(p2Idx(1)+1:end) = newDVtt(p2Idx(1)+1:end)-offset;
                            end
                            clear offset; clear decayEnd; 
                        end
                        
                        decayP1(p1Idx(0.2*100):p1Idx(0.2*100)+pars.stepAfterOffset-1)=p1ERP(p1Idx(0.2*100))-[1:pars.stepAfterOffset]/(pars.stepAfterOffset)*p1ERP(p1Idx(0.2*100));
                        decayP1(p1Idx(0.2*100)+pars.stepAfterOffset:end) = 0;
                            
                        if pars.biased == 1, thisDV = DVtt; end
                    end
                elseif p == 2
                    decayP2(p2Idx(end):p2Idx(end)+pars.stepAfterOffset-1)=p2ERP(p2Idx(end))-[1:pars.stepAfterOffset]/(pars.stepAfterOffset)*p2ERP(p2Idx(end));
                    decayP2(p2Idx(end)+pars.stepAfterOffset:end) = 0;
                    
                    newDVtt(p2Idx(end):p2Idx(end)+pars.stepAfterOffset-1)=newDVtt(p2Idx(end))-[1:pars.stepAfterOffset]/(pars.stepAfterOffset)*newDVtt(p2Idx(end));
                    newDVtt(p2Idx(end)+pars.stepAfterOffset:end)=0; % 100 ms linear drop after bound crossing:

                    if pars.biased == 0 || p1DV == 0%If not biased OR no bound hit at all
                        p2DV = thisDV(p2Idx(end)); 
                    elseif pars.biased == 1 & p1DV ~= 0 %If no bound hit, do not compute
                        p2DV = 0;
                        p2ERP = repmat(0, length(p2ERP), 1)'; %NO accum if bound already hit
                    end
                    
                end
                
            end
        end
        
        %Ensure dimensions are right
        DVtt = DVtt(1:length(t));
        thisDV = thisDV(1:length(t));
        p1ERP = p1ERP(1:length(t)); 
        p2ERP = p2ERP(1:length(t)); 
        
        totalDV = p1DV + p2DV;
        
        %2-p & 1-p accuracies
        if totalDV > 0 %+stt> 0
            ch(n,1)=1; %Assume correct if DV>0 & no choice made
        elseif totalDV <= 0 %+stt<= 0
            ch(n,1)=0; %Assume error if DV<0 & no choice made
        end
        
        if p1DVbias > 0%+stt> 0
            ch_1p(n,1)=1; %Assume correct if DV>0 & no choice made
        elseif p1DVbias <=0%+stt<= 0
            ch_1p(n,1)=0; %Assume error if DV<0 & no choice made
        end
        
        %Compute perfect accuracy 
        if p1DV_perf+p2DV_perf > 0
        ch_perfect(n,1) = 1; 
        else
            ch_perfect(n,1) = 0; 
        end
        
        p1plusp2(c,n,:) = p1ERP  + p2ERP;
        wfDV(c,n,:)= DVtt; 
        allMBL(c,n,:) = MBL;
        decayDV(c,n,:) = newDVtt; 
        bounds(c,n,:) = b; %To plot with simulated ERPs
        simp1(c,n,:) = p1ERP;
        simp2(c,n,:) = p2ERP;
        simp1_decay(c,n,:) = decayP1;
        simp2_decay(c,n,:) = decayP2;
        sim_stt(c,n,:) = stt; %Random baseline shift - to use for plotting potential preStim alpha effects. Low bias --> bigger ERPs
        allp1DV(c,n,:) = p1DV;
        allp2DV(c,n,:) = p2DV;
        alltotalDV(c,n,:) = totalDV;
        allp1DV_end(c,n,:) = p1DV_end; %total ev at the end of P1
        allp1DV_perf(c,n,:) = p1DV_perf; %total ev at the end of P1
        allp2DV_perf(c,n,:) = p2DV_perf; %total ev at the end of P1
         
        clear DVtt; clear p1ERP; clear p2ERP; clear newDVtt;clear decayP1; clear decayP2; 
        
    end
    ch(ch == 0) = -1;
    ch_1p(ch_1p == 0) = -1;
    CH{c} =  ch; %Simulated choices
    CH_p{c} = 0.5+0.5.*erf(ch);
    CH_perf{c} = ch_perfect; %Simulated perfect accumulator choices
    CH_1p{c} = 0.5+0.5.*erf(ch_1p);
    clear ch; clear ch_1p; clear ch_noB;
    if ~(length(CH{c}) == length(CH_p{c}) & length(CH{c}) == length(CH_1p{c}))
        disp(['uh oh! Mismatch at ' pars.cLab{c}])
    end
end

% simulated data matrix to be returned:
simdat=[];
for c=condIdx(1:end)
    simdat = [simdat; c*ones(N,1) CH{c} CH_p{c} CH_1p{c}];
end

%% Single structure to output
sims.wfDV = wfDV;       	%Single trial DV, assuming no signal fall-down when dots turn blue
sims.decayDV = decayDV;     %Simulated CPP, assuming CPP turns off and falls down when 1) dots turn blue or 2) bound is hit
sims.MBL = allMBL;       	%Simulated MBL, which never decays, staying stable both when bound is hit and when dots turn blue
sims.decayP1 = simp1_decay; %Single trial P1 ERP, assuming CPP shuts down when dots turn blue
sims.decayP2 = simp2_decay; %Single trial P2 ERP, assuming CPP shuts down when dots turn blue

sims.p1pp2 = p1plusp2;
sims.p1DV = allp1DV;        %Single trial P1 ERP, assuming CPP does NOT shut down when dots turn blue, but only when a bound is hit
sims.p2DV = allp2DV;        %Single trial P2 ERP, assuming CPP does NOT shut down when dots turn blue, but only when a bound is hit

sims.p1DV_end = allp1DV_end;         
sims.p1DV_perf= allp1DV_perf;         
sims.p2DV_perf= allp2DV_perf;         
sims.totalDV = alltotalDV;

sims.bounds = bounds;       %Estimated bounds
sims.simp1 = simp1;         %P1 DV
sims.simp2 = simp2;         %P2 DV
sims.simdat = simdat;       %Simulated accuracies: condition, binary accuracy, prob. accuracy
sims.bndTDist = bndTDist;   %Bound crossing time
sims.bndTP1 = bndTP1;       %Was bound crossed after P1?
sims.bndTP2 = bndTP2;       %Was bound crossed after P2?
sims.stt = sim_stt;         %Random baseline shift (not used)

end