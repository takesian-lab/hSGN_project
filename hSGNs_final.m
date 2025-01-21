% Code for analyzing ephys data from hSGNs

%% User defined parameters:
clear;
close all;
%1) Choose data type and path:

%pathname = '\\apollo\research\ENT\Takesian Lab\Lucas\Manuscripts\4_hSGNs\for GitHub\Example data\Current clamp example'; 
pathname = '\\apollo\research\ENT\Takesian Lab\Lucas\Manuscripts\4_hSGNs\for GitHub\Example data\Voltage clamp example'; 


%2) Choose whole-cell recording mode:
%clamp_type = 'CC'; 
clamp_type = 'VC';

%3) Define specific parameters:
sampling_frequency = 10; % in kHz
first_step = -100; %first step of current or voltage injection to analyze(in pA or mV) 
%-80 pA and -100 mV in these exmaples
last_step = 40; %last step of current or voltage injection to analyze (in pA or mV) 
%400 pA and 40 mV in these examples
increase_step = 10; %current or voltage injection step increase (in pA or mV)
stimStart = 500; %in ms; time of injected current
stimLength = 1000; %in ms; duration of the injected pulse

% Define only for CC:
CC.obs_threshold = 12;%sweep in which first AP is observed
CC.plot_all = 1; %to plot all traces (0 or 1)
CC.plot_ind = 1; %to plot individual traces (0 or 1)

% Define only for VC:
VC.first_Na = 6; %sweep in which first Na current is observed
VC.preStim = stimStart-50; %50 pts before the initiation of the stimulus

% For RC calculation (only calculated in VC):
pulseAmp = -10; %in mV (VC)
pulseStart = 100; %in ms
pulseLength = 100; %in ms
baselineStart = 1; %in ms
baselineEnd = 99; %in ms


%% FROM HEREON IT RUNS AUTOMATICALLY


%% Load data, define sampling rate and remaining variables:

cd(pathname)
file = dir([pathname '/*.txt' ]);
data = load(fullfile(pathname,file.name));

samplingPoints = 1:size(data,1); %sampling points
pointspermillisecond = ((sampling_frequency)*10^3)/1000; %1000 ms in one s
samplingRate = 1/pointspermillisecond; %a point every (samplingRate) ms; sampling rate should be 0.1 if 10 kHz
%Sampling rate in ms; 0.1 ms = 100 us (10kHz - 10000 points per second --> 10 points per ms --> one point every 0.1 ms)
time = samplingPoints*samplingRate; 
stimEnd = stimStart + stimLength;

%% Create variables aacording to VC or CC and sampling rate:

datastimStart = stimStart/samplingRate; 
datastimEnd = stimEnd/samplingRate; 
stim_inject = first_step:increase_step:last_step;

if strcmp(clamp_type,'CC')    
    indexZero = find(stim_inject == 0); %sweep in which no current was injected
    threshold_current = stim_inject(CC.obs_threshold); % to report the current value of the observer threshold
    dataSAG = nan(1,size(data,2));
    dataminTime = nan(1,size(data,2));
    datasteadyTime = nan(1,size(data,2));
    dataminAmplitude = nan(1,size(data,2));
    datasteadyAmplitude = nan(1,size(data,2));
    datapeakAmplitude = cell(size(data,2),1);
    datapeakTime = cell(size(data,2),1);
    datathresholdTime = cell(size(data,2),1);
    datathresholdAmplitude = cell(size(data,2),1);
    datatroughAmplitude = cell(size(data,2),1);
    datatroughTime = cell(size(data,2),1);
    dataAPheight = cell(size(data,2),1);
    dataISI = cell(size(data,2),1);
    datahalfwidth = cell(size(data,2),1);
    datahalftime1 = cell(size(data,2),1);
    datahalftime2 = cell(size(data,2),1);
    spike_result = cell(length(stim_inject),16); %based on the number of output
    %variables that I'm currently using   
elseif strcmp(clamp_type,'VC')
    % %Load required filter:
    load 'Lowpass equiripple filter.mat' %should be in the same workspace
    datapreStim = VC.preStim/samplingRate;  
    leak_current_real = nan(1,3); % 3 points to do the linear fit
    leak_current = nan(1,size(stim_inject,2));
    current = nan(1,size(stim_inject,2));
    Na_current = nan(1,size(stim_inject,2));    
    RC_res = cell(length(stim_inject),4);
end

%% Main analysis

if strcmp(clamp_type,'CC')
    
    for i = 1:size(data,2) %CC.obs_threshold
        trace_data = smooth(data(:,i));
        baseline(i) = mean(trace_data(1:pulseStart-1)); %baseline before RC check
        
        if i < CC.obs_threshold
            spike_result{i,1} = stim_inject(i);
            spike_result{i,2} = baseline(i);
            spike_result{i,3} = NaN;
            if i < indexZero %for traces below I_inject == 0
                [minAmplitude, minTime] = min(trace_data(datastimStart:datastimStart+3000)); %min mV in the 1st ~third of the step
                [steadyAmplitude, steadyTime] = min(trace_data(datastimEnd-250:datastimEnd-1)); %min mV in the last 25 ms
                dataSAG(i) = steadyAmplitude-minAmplitude;
                minTime = datastimStart-1 + minTime;
                steadyTime = datastimEnd-250 + steadyTime;
                dataminTime(i) = minTime;
                datasteadyTime(i) = steadyTime;
                dataminAmplitude(i) = trace_data(minTime);
                datasteadyAmplitude(i) = trace_data(steadyTime);
                spike_result{i,3} = dataSAG(i);
            end
            for j = 4:size(spike_result,2)
                spike_result{i,j} = NaN;
            end
        elseif i >= CC.obs_threshold
            %Amplitude and time of peaks:
            [peakAmplitude, peakTime] = findpeaks(trace_data(datastimStart:datastimEnd-1),'MinPeakDistance',10,'MinPeakProminence',25);
            peaksVector = [peakAmplitude, peakTime]; %vector with the peaks found
            datapeakAmplitude{i} = peakAmplitude;
            peakTime = datastimStart-1 + peakTime; %when adding the peakTime in points to datapulseStart everything shitfs +1 in x, so -1 to shift it back to the right position
            datapeakTime{i} = peakTime;
            derivdata = diff(trace_data,1); %1st derivative to calculate threshold of each AP
            derivdata(20000)=derivdata(end); %extra point to make the array even
            a = derivdata(datastimStart+50:peakTime(1)); %region in which we'll look for first onset           
            limit = 0.075*max(a);
            [~,limit_index]=min(abs(a-limit));            
            limit2 = 0.075*max(derivdata);
            
            troughTime = [];
            troughAmplitude = [];
            thresholdTime = [];
            thresholdTimeSubseq = [];
            thresholdAmplitude = [];
            
            if size(peaksVector,1) > 1
                %Trough:
                %all except the last trough, since it could be after the
                %end of the pulse:
                for j = 1:length(peakTime)-1
                    [troughAmplitude(j), troughTime(j)] = min(trace_data(peakTime(j):peakTime(j+1)));
                    troughTime(j) = peakTime(j)-1 + troughTime(j);
                end
                %now I calculate the last trough
                [troughAmplitudeLast, troughTimeLast] = min(trace_data(peakTime(end):datastimEnd-1));
                troughTimeLast = peakTime(end) +troughTimeLast;
                %...and add it to the array only if it's before the end of
                %the pulse
                if troughTimeLast < datastimEnd
                    troughAmplitude = [troughAmplitude troughAmplitudeLast];
                    troughTime = [troughTime troughTimeLast];
                end
                
                datatroughTime{i} = troughTime;
                datatroughAmplitude{i} = troughAmplitude;
                
                %Onset:
                %first onset:
                thresholdTimeFirst = datastimStart + limit_index; %c(end);;% - 20; %-20 POINTS SEEMS TO WORK BEST
                
                %subsequent onsets:
                for j = 1:length(peakTime)                    
                    if j < length(peakTime) %to avoid going past the last peak
                        [~, thresholdTimeSubseq(j)] = min(abs(derivdata(troughTime(j):peakTime(j+1))-limit2));
                        thresholdTimeSubseq(j) = troughTime(j)-100 + thresholdTimeSubseq(j);                        
                    end                    
                end
                thresholdTime = [thresholdTimeFirst thresholdTimeSubseq-10];
                thresholdAmplitude = trace_data(thresholdTime);
                datathresholdAmplitude{i} = thresholdAmplitude;
                datathresholdTime{i} = thresholdTime;
                
            else %if only one AP
                %Trough:
                [troughAmplitude, troughTime] = min(trace_data(peakTime:datastimEnd-1));
                troughTime = peakTime-1+troughTime;
                datatroughTime{i} = troughTime;
                datatroughAmplitude{i} = troughAmplitude;
                
                %Onset:
                [~, thresholdTime] = min(abs(derivdata(datastimStart:peakTime)-limit));
                thresholdTime = datastimStart-1 + thresholdTime;
                thresholdAmplitude = trace_data(thresholdTime);
                datathresholdAmplitude{i} = thresholdAmplitude;
                datathresholdTime{i} = thresholdTime;
            end
            
            %Calculation of avg AP height (peak to trough):
            if size(peaksVector,1) > 1
                AP_height = [];
                ISI = [];
                half_width = []; %these three are defined as empty arrays to ensure that they reset in every cycle, since their size varies!
                halfTime1 = [];
                halfTime2 = [];
                %AP height:
                if length(troughTime) == length(peakTime) %if it was possible to calculate a trough for every peak (that is, if last trough doesn't fall after pulseEnd)
                    for j = 1:length(peakTime)
                        AP_height(j) = abs(peakAmplitude(j)-troughAmplitude(j));
                    end
                else
                    for j = 1:length(peakTime)-1
                        AP_height(j) = abs(peakAmplitude(j)-troughAmplitude(j));
                    end
                end
                
                %ISI:
                for j = 1:length(thresholdTime)
                    if j < length(thresholdTime)
                        ISI(j) = thresholdTime(j+1)-thresholdTime(j); %in points
                    end
                end
                if length(ISI) > 1
                    for k = 1:length(ISI)-1 %this will calculate the differences between ISIs across an AP
                        ISIindex(k) = (ISI(k+1)-ISI(k))/(ISI(k+1)+ISI(k)); 
                    end                    
                    adapt_index = (1/(length(ISI)-1)*sum(ISIindex)); %the closest to 0 the less adaptation in ISIs
                else
                    adapt_index = 0; %only applies if more than 2 AP
                end
                %Half-width:
                if length(troughTime) == length(peakTime)
                    for j = 1:length(peakTime)
                        t1_data(j) = 0.5*(abs(peakAmplitude(j)-thresholdAmplitude(j)));
                        t1_data(j) = t1_data(j) + thresholdAmplitude(j);
                        [~, halfTime1(j)] = min(abs(trace_data(thresholdTime(j):peakTime(j))-t1_data(j)));
                        halfTime1(j) = thresholdTime(j)-1 + halfTime1(j);                        
                        t2_data(j) = peakAmplitude(j)-(peakAmplitude(j)-t1_data(j));
                        [~, halfTime2(j)] = min(abs(trace_data(peakTime(j):troughTime(j))-t2_data(j)));
                        halfTime2(j) = peakTime(j)-1 + halfTime2(j);                        
                        half_width(j) = halfTime2(j)-halfTime1(j);
                    end                    
                else
                    for j = 1:length(peakTime)-1
                        t1_data(j) = 0.5*(abs(peakAmplitude(j)-thresholdAmplitude(j)));
                        t1_data(j) = t1_data(j) + thresholdAmplitude(j);
                        [~, halfTime1(j)] = min(abs(trace_data(thresholdTime(j):peakTime(j))-t1_data(j)));
                        halfTime1(j) = thresholdTime(j)-1 + halfTime1(j);                        
                        t2_data(j) = peakAmplitude(j)-(peakAmplitude(j)-t1_data(j));
                        [~, halfTime2(j)] = min(abs(trace_data(peakTime(j):troughTime(j))-t2_data(j)));
                        halfTime2(j) = peakTime(j)-1 + halfTime2(j);                        
                        half_width(j) = halfTime2(j)-halfTime1(j);
                    end
                    
                end
                
            else %if only one AP
                %AP height and ISI:
                if length(troughTime) == length(peakTime)
                    AP_height = abs(peakAmplitude-troughAmplitude);
                else
                    AP_height = NaN;
                end
                ISI = NaN;
                adapt_index = NaN;
                
                %Half-width:
                t1_data = 0.5*(abs(peakAmplitude-thresholdAmplitude));
                t1_data = t1_data + thresholdAmplitude; %target mV to look for t1                
                [~, halfTime1] = min(abs(trace_data(thresholdTime:peakTime)-t1_data));
                halfTime1 = thresholdTime-1 + halfTime1;                
                t2_data = peakAmplitude-(peakAmplitude-t1_data); %same as t1_data in fact, but nice to have it in case I want to plot it
                [~, halfTime2] = min(abs(trace_data(peakTime:troughTime)-t2_data));
                halfTime2 = peakTime-1 + halfTime2;                
                half_width = halfTime2-halfTime1; %in points
                
            end
            dataAPheight{i} = AP_height;
            dataISI{i} = ISI;
            datahalfwidth{i} = half_width;   
            datahalftime1{i} = halfTime1;
            datahalftime2{i} = halfTime2;           
            
            spike_result{i,1} = stim_inject(i);
            spike_result{i,2} = baseline(i);
            spike_result{i,3} = NaN;
            spike_result{i,4} = size(peaksVector,1); %stores number of detected APs in each trace
            spike_result{i,5} = size(peaksVector,1)/(length(datastimStart:datastimEnd-1)*samplingRate/1000); %spike freq in spikes/sec for f-i curve; denominator = 1 sec
            spike_result{i,6} = datathresholdAmplitude{i}(1); %onset amplitude (mV) for 1st spike
            spike_result{i,7} = (datathresholdTime{i}(1)-datastimStart)*samplingRate; %onset latency (ms) for 1st spike
            spike_result{i,8} = datapeakAmplitude{i}(1); %absolute peak amplitude (mV) for 1st spike
            spike_result{i,9} = (datapeakTime{i}(1)-datastimStart)*samplingRate; %peak latency (ms) for 1st spike
            spike_result{i,10} = dataAPheight{i}(1); %AP height (mV) for 1st AP in the step
            spike_result{i,11} = dataAPheight{i}(end); %AP height (mV) for the last AP, if there's more than one AP in the step
            spike_result{i,12} = mean(dataAPheight{i}); %avg AP height (mV)
            spike_result{i,13} = dataISI{i}(1)*samplingRate; %1st ISI (between 1st two APs) in ms
            spike_result{i,14} = mean(dataISI{i}*samplingRate); %avg ISI in ms
            spike_result{i,15} = adapt_index; %adaptation index: how ISI changes from AP to AP (is the train accelerating or not?)
            spike_result{i,16} = mean(datahalfwidth{i}*samplingRate); %avg HW in ms                    
        end       
    end
    
    spike_result_str = {'I_step', 'Vm', 'SAG_mV', 'Spike_count', 'Spikes/sec','ThresholdAmp_1st_AP_mV','ThresholdLatency_1stAP_ms','PeakAmp_1st_AP_mV','PeakLatency_1stAP_ms', 'AP_height_1stAP', 'AP_height_lastAP', 'avgAP_height_mV','1st_ISI_ms','avgISI_ms','Adapt_index','avgHW_ms'};
    spike_result = [spike_result_str; spike_result];
    CC_results.currentClamp = spike_result;
    CC_results.thresholdinpA = threshold_current;
    CC_results.AvgVm = mean(baseline');   
    
    if CC.plot_all == 1
        %Plot all the traces together:
        for i = 1:size(data,2)            
            subplot1 = subplot(5,10,i);
            trace_data = smooth(data(:,i));            
            if i >= CC.obs_threshold
                plot(time,trace_data,'k','LineWidth',0.5,'Parent',subplot1);
                hold on; plot(datapeakTime{i}*samplingRate, datapeakAmplitude{i}, 'bo','Parent',subplot1)
                hold on; plot(datatroughTime{i}*samplingRate, datatroughAmplitude{i}, 'bo','Parent',subplot1)
                hold on; plot(datathresholdTime{i}*samplingRate, datathresholdAmplitude{i}, 'ro','Parent',subplot1)
                hold on; plot(datahalftime1{i}*samplingRate, trace_data(datahalftime1{i}),'go','Parent',subplot1)
                hold on; plot(datahalftime2{i}*samplingRate, trace_data(datahalftime2{i}),'go','Parent',subplot1)
            elseif i < CC.obs_threshold
                plot(time,trace_data,'k','LineWidth',0.5,'Parent',subplot1);
                hold on;
                if i < indexZero
                    plot(dataminTime(i)*samplingRate, dataminAmplitude(i), 'ro', 'Parent',subplot1);
                    plot(datasteadyTime(i)*samplingRate, datasteadyAmplitude(i), 'bo', 'Parent',subplot1);
                end
            end
            xlabel('Time (msec)');
            ylabel('Membrane potential (mV)'); hold on;
            maxvalue = max(trace_data(4000:end));
            minvalue = min(trace_data(4000:end));
            ylim([minvalue-50 maxvalue+50]);
        end
    end
    
    if CC.plot_ind == 1
        %Plot traces individually:
        for i = 1:size(data,2)  %CC.obs_threshold
            trace_data = smooth(data(:,i));
            if i < indexZero
                figure;
                plot(time,trace_data,'k','LineWidth',0.5);
                hold on;
                plot(dataminTime(i)*samplingRate, dataminAmplitude(i), 'ro');
                plot(datasteadyTime(i)*samplingRate, datasteadyAmplitude(i), 'bo');
            elseif i >= CC.obs_threshold
                figure;
                plot(time,trace_data,'k','LineWidth',0.5);
                hold on;
                plot(datapeakTime{i}*samplingRate, datapeakAmplitude{i}, 'bo')
                plot(datatroughTime{i}*samplingRate, datatroughAmplitude{i}, 'bo')
                plot(datathresholdTime{i}*samplingRate, datathresholdAmplitude{i}, 'ro')
                plot(datahalftime1{i}*samplingRate, trace_data(datahalftime1{i}),'go');
                plot(datahalftime2{i}*samplingRate, trace_data(datahalftime2{i}),'go');
                plot(time,derivdata); %if you want to see the 1st derivative trace
                %from which AP parameters were estimated
            end
            xlabel('Time (msec)');
            ylabel('Membrane potential (mV)'); hold on;
            maxvalue = max(trace_data(4000:end));
            minvalue = min(trace_data(4000:end));
            ylim([minvalue-50 maxvalue+50]);                      
        end
    end
    
    clearvars -except CC_results pathname clamp_type
    
elseif strcmp(clamp_type,'VC')
    
    %Compute leak current:
    index_leak = [find(stim_inject==-80) find(stim_inject==-70) find(stim_inject==-60)];    
    for i = index_leak(1):index_leak(end)
        trace_data = smooth(data(:,i));
        avg_final_curr_leak = mean(trace_data(datastimEnd-10/samplingRate:datastimEnd)); %avg current in the last 10 ms of the step
        preStimBaseline_leak = mean(trace_data(datapreStim:datastimStart)); %baseline 50 ms before stimulus
        
        leak_current_real(i-2) = avg_final_curr_leak-preStimBaseline_leak; %amplitude between avg_final_curr and pre stimulus
        %I did i-2 because 1st i would otherwise be position 3, and I'd end up
        %with an array with two empty spots that wouldn't allow me to do a
        %proper fit
        %baseline(i) = mean(data(1:pulseStart));
    end    
    
    %Fit to obtain the equation for leak current:
    f = fit(stim_inject(index_leak(1):index_leak(end))',leak_current_real','poly1');    
    %Obtain the leak current for all the steps:
    for i = 1:size(stim_inject,2)
        leak_current(i) = f.p1*stim_inject(i)+f.p2;
    end
    
    %Calculate final currents (substracting the leak current) and RC parameters:   
    for i=1:size(stim_inject,2) %up to 40 mV until I figure out what happened with the rest
        trace_data2 = smooth(data(:,i));
        
        %compute baselines for each trace
        %baseline(i) = mean(data(1:pulseStart)); % baseline before Rm step
        %preStimBaseline(i)= mean(data(datapreStim:datastimStart)); %baseline 50 ms before stimulus
        
        %Compute Rs and Rm for each trace
        [Rin(i), Rs(i), Cm(i), error(i), tau(i)] = calcRs_VC_hSGNs(data(:,i), samplingRate, pulseAmp, pulseStart, pulseLength, baselineStart, baselineEnd);
        
        %Calculate current at the end of the step and substract leak:
        avg_final_curr = mean(trace_data2(datastimEnd-10/samplingRate:datastimEnd)); %avg current in the last 10 ms of the step
        preStimBaseline = mean(trace_data2(datapreStim:datastimStart)); %baseline 50 ms before stimulus
        curr = avg_final_curr-preStimBaseline;        
        current(i) = (avg_final_curr-preStimBaseline)-leak_current(i);
    end
    
    %Calculate Na currents:
    for i=VC.first_Na:size(stim_inject,2)
        trace_data3 = smooth(filter(Hd1,data(:,i)));
        [Na_peak, Na_peak_time] = min(trace_data3(datastimStart:datastimStart+100/samplingRate)); %look for the min between were the pulse starts and 100 ms later to find the peak of Na currents
        Na_peak_time = datastimStart+Na_peak_time-1;        
        [transient_max, transient_max_time] = max(trace_data3(datastimStart:Na_peak_time));
        transient_max_time = datastimStart+transient_max_time-1;        
        inflection = diff(trace_data3(transient_max_time:Na_peak_time),2);
        [~, real_inflection_time] = min(abs(0-inflection));
        real_inflection_time = transient_max_time+real_inflection_time;       
        Na_start = real_inflection_time;
        Na_current(i) = Na_peak-trace_data3(Na_start);        
    end
    
    RC_res{i,1} = Cm(i);
    RC_res{i,2} = Rin(i);
    RC_res{i,3} = Rs(i);
    RC_res{i,4} = tau(i);
    
    result_Array = [current; stim_inject; Na_current];
    result_Array = result_Array';
    colNum_result_Array = size(result_Array,2);
    
    for k=1:colNum_result_Array
        meanVar = mean(result_Array(:,k),'omitnan');
        stdVar = std(result_Array(:,k),'omitnan');
        semVar = stdVar/sqrt(sum(~isnan(result_Array(:,k))));
        cvVar = stdVar/abs(meanVar);
        Var(1,k) = meanVar;
        Var(2,k) = stdVar;
        Var(3,k) = semVar;
        Var(4,k) = cvVar;
    end
    
    newRes_Array = [result_Array; Var];    
    result_CellArray = num2cell(newRes_Array);
    t_str = {'Current_no_leak', 'Voltage_step_mV','Na_current'};    
    result_CellArray = [t_str; result_CellArray];
    result_CellArray = [repmat({'NaN'},size(result_CellArray,1),1) result_CellArray];    
    s_str = {'Mean','SD','SEM','CV'};
    pos = [3 2 1 0];
    
    for f = 1:size(pos,2)
        result_CellArray{end-pos(f),5} = {'NaN'};
        result_CellArray{end-pos(f),6} = {'NaN'};
        result_CellArray{end-pos(f),7} = {'NaN'};
        result_CellArray{end-pos(f),1} = s_str{f};
    end
    
    VC_results.voltageClamp = result_CellArray;
    VC_results.Cm = mean([RC_res{:,1}]);
    VC_results.Rm = mean([RC_res{:,2}]);
    VC_results.Rs = mean([RC_res{:,3}]);
    VC_results.Tau = mean([RC_res{:,4}]);
    clearvars -except VC_results pathname clamp_type
end
%% Store results:

if strcmp(clamp_type,'CC')
    results = 'CC_results';
    pathname = [pathname '\'];
    resultvar = [pathname 'CC_results.mat'];
    save(resultvar, results)
elseif strcmp(clamp_type,'VC')
    results = 'VC_results';
    pathname = [pathname '\'];
    resultvar = [pathname 'VC_results.mat'];
    save(resultvar, results)
    
end

