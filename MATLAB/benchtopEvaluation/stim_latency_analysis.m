clear
%% Add paths
addpath('.../violinplot'); % https://github.com/bastibe/Violinplot-Matlab.git

%% TTL pulse is on ch1 every StimulusCode==1, stim was to ch 3-5 when StimulusCode==1
path = 'path/to/data/';
fileNames = ['ttl_1.dat'; 'ttl_2.dat'; 'ttl_3.dat'];
latencies = cell(3,3);
stimModes = dictionary([0 1 2],["Volatile Commands" "Persistent Commands" "Persistent Functions"]);

%% Extract latencies
blockSize = nan(size(fileNames,1),1);
x_1 = 800;
for f = 1:size(fileNames,1)
  fprintf("Analyzing file " + fileNames(f,:) + "\n")
  [sig, st, parm] = load_bcidat([path fileNames(f,:)]);
  ind = 1:size(sig,1);
  starts = ind(diff(st.StimulusCode==1)==1);
  t = ind/parm.SamplingRate.NumericValue;
  ret = getPulseTimes(sig, t, starts, x_1, st, true);
  latencies{f,1} = ret(1,:);
  latencies{f,3} = ret(2,:);
  latencies{f,2} = stimModes(parm.StimulationMode.NumericValue);
  blockSize(f) = parm.SampleBlockSize.NumericValue;
end

if all(diff(blockSize) == 0)
    blockSize = blockSize(1);
else
    error('Recordings were acquird with a  varying sample block size')
end

acqLatency = reshape(cell2mat(latencies(:,3)),1,[]);
% Subtract the block size
acqLatency(acqLatency >= blockSize) = acqLatency(acqLatency >= blockSize) -blockSize;

fprintf("Done analyzing files\n")

%% visualize Combined stimulation and ttl (acquisition) latencies - Histogram
figure
hold on
bins = [6 6 50];
colors = [ "#456990","#ef767a", "#49beaa","#ffdd00"];
disp(newline)
% Acquisition latency 
histogram(acqLatency,6, 'FaceColor', colors(4), 'FaceAlpha', 0.75)
fprintf("Acquisition latency: " + mean(acqLatency) + " +/- " + std(acqLatency) + "\n")
% Stimulation latencies
for i = 1:size(latencies,1)
  histogram(latencies{i,1},bins(i), 'FaceColor', colors(i), 'FaceAlpha', 0.5)
  fprintf(latencies{i,2} + ": " + mean(latencies{i,1}) + " +/- " + std(latencies{i,1}) + "\n")
end
legend(["Acquisition latency" latencies{:,2}])
xlabel("Time (ms)")
ylabel("Occurences")
title("Acquisition and stimulation latencies")
% subtitle(sprintf("Acquisition sample block size: %d samples",blockSize))
% fontname("Myrad Pro")
fontsize( "increase" )

%% visualize Combined stimulation and ttl (acquisition) latencies - violin plots
% Data vector
allData = [ ...
    acqLatency(:); ...
    latencies{1,1}(:); ...
    latencies{2,1}(:); ...
    latencies{3,1}(:)];

% Matching group labels
allCats = [ ...
    repelem({'Acquisition latency'}, numel(acqLatency))'; ...
    repelem(latencies{1,2}, numel(latencies{1,1}))'; ...
    repelem(latencies{2,2}, numel(latencies{2,1}))'; ...
    repelem(latencies{3,2}, numel(latencies{3,1}))'];

color = [239 118 122;...
         73 190 170;...
         255 221 0;...
         69 105 144;...
        ]./255;

figure
vs = violinplot(allData, allCats, 'HalfViolin','right',...% left, full
    'QuartileStyle','shadow',... % boxplot, none
    'DataStyle', 'histogram',... % scatter, none
    'ShowNotches', false,...
    'ShowMean', false,...
    'ShowMedian', true,...
    'ViolinColor', color,...
    'Orientation', 'vertical');

ylabel('Latency [ms]')
grid on
fontsize(16,'points')

%% visualize
figure
hold on
bins = [6 6 50];
colors = ["#ef767a", "#456990", "#49beaa"];
disp(newline)
for i = 1:size(latencies,1)
  histogram(latencies{i,1},bins(i), 'FaceColor', colors(i), 'FaceAlpha', 0.5)
  fprintf(latencies{i,2} + ": " + mean(latencies{i,1}) + " +/- " + std(latencies{i,1}) + "\n")
end
legend(latencies{:,2})
xlabel("Time (ms)")
ylabel("Occurences")
title("Stimulation latencies for varying stimulation modes")
%fontname("Lucida Sans")
fontsize( "increase" )
%% visualize persistent modes
figure
hold on
bins = [6 6];
colors = ["#ef767a", "#456990"];
disp(newline)
for i = 1:size(bins,2)
  histogram(latencies{i,1},bins(i), 'FaceColor', colors(i), 'FaceAlpha', 0.5)
  fprintf(latencies{i,2} + ": " + mean(latencies{i,1}) + " +/- " + std(latencies{i,1}) + "\n")
end
legend(latencies{:,2})
xlabel("Time (ms)")
ylabel("Occurences")
title("Stimulation latencies for persistent stimulation modes")
%fontname("Lucida Sans")
fontsize( "increase" )
xlim([0 25])

%% visualize ttl latencies (acquisition latency)
figure
hold on
bins = [5 50 50];
colors = ["#ef767a", "#456990", "#49beaa"];
disp(newline)
for i = 1:size(latencies,1)
  histogram(latencies{i,3},bins(i), 'FaceColor', colors(i), 'FaceAlpha', 0.5)
end
legend(latencies{:,2})
xlabel("Time (ms)")
ylabel("Occurences")
title("TTL latencies for varying stimulation modes")
fontname("Lucida Sans")
fontsize( "increase" )



%% get peaks
function ret = getPulseTimes(sig, t, starts, x_1, st, draw)
% Inputs:
% sig: input signal
% t: time indices
% starts: indices of the TTL pulses
% x_1 : window (specified in samples)
% st : states struct
% draw: boolean value for plotting
%
% Ouputs:
% ret: returns Nx2 array where:
%   :(n,1) are latencies between the detected peaks and TTL pulses
%   :(n,2) are time indices of the detected TTL pulses

%assumes ttl pulses are on ch 1, stimulation is on ch 3
figure
thresh = 50; %doesn't need to be too high
ttlPulses = zeros(size(starts));
stimPulses = zeros(size(ttlPulses));
i = 1;
for s = starts
  n = s:s+x_1;
  x = t(n)-t(s);
  data = sig(n,1);
  ddata = diff(data);
  plot(x(1:end-1), ddata)
  allTTLPeaks = x(ddata>thresh);
  ttlPulses(i) = allTTLPeaks(1); %first peak
  dSig = diff(sig(n,3));
  allStimPeaks = x(dSig>thresh);
  hold off
  plot(x(1:end-1), ddata)
  hold on
  plot(x(1:end-1),dSig)
  if length(allStimPeaks) < 1
    allStimPeaks
  else
    stimPulses(i) = allStimPeaks(1); %first peak
  end
  if draw
    colors = ["#ef767a", "#456990", "#49beaa"];
    %add to existing figure
    hold off
    plot(x*1000, st.StimulusCode(n)*1800,'LineWidth', 2, 'Color', colors(1))
    hold on
    plot(x*1000, data,'LineWidth',2, 'Color', colors(2))
    plot(x*1000,sig(n,3),'LineWidth',2, 'Color', colors(3))
    legend("Trigger","TTL","Stimulation")
    xlabel("Time (ms)")
    ylabel("Amplitude")
    title("Latencies after BCI2000 trigger")
    box off

  end
  i = i + 1;
end
latencyMs = (stimPulses - ttlPulses)*1000;
acquisitionMs = ttlPulses*1000;
ret = [latencyMs; acquisitionMs];
end



