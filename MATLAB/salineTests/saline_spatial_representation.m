clear all; close all; clc
addpath('/Users/lampert.frederik/Documents/MATLAB/mex');
addpath('/Users/lampert.frederik/Documents/MyCodes/MATLAB/functions');
addpath('/Users/lampert.frederik/Documents/MATLAB/mnl_ieegBasics/functions');
%% ~~~~~~ select paths and load the data ~~~~~~~~~~~~~
fig_dir = '/Users/lampert.frederik/Documents/Mayo/CorTec/Figures/saline-test/implant04/stim/ses03/external';    % specify path where to save the figures
ddr = '/Users/lampert.frederik/Documents/Mayo/CorTec/data/other_recordings/Saline-tests'; % select data directory to open                          
sp = filesep;
%% ---------- Select file ----------------------
[fl, fl_pth] = uigetfile([ddr sp '*.dat'], 'Select a File');                 % select file
fl_id = strsplit(fl, '.');                                                 % get file name and extension
[signal, states, parameters, total_samples, file_samples ] ...      
    = load_bcidat([fl_pth, fl], '-calibrated');         % load using MEX package 
data = interpolate_lost_samples(signal, states);        % interpolate lost values
% ~~~~~~~~~~ Specify constants & Get parameters ~~~~~~~~~~~~~~~~
FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
T_ax = linspace(0,size(signal,1)/FS, size(signal,1));   % Creating time axis for plotting signal
amps = {'57.5dB',  '51.5dB', '45.5dB',  '39.5dB'};      % Specify available amplification values
ref_ch = parameters.ReferenceCh.NumericValue;           % Get refeence channel


%% ##############
SAVE_FIG = true;
chnls = [1:32];           % chanels selected
%##############
%% Plot the raw signal
displayChannels(data, FS, ref_ch, chnls, fl_id{1});
if SAVE_FIG
    % kjm_printfig([fig_dir fl_id{1} '_grid' num2str(grd)],[12 10]);
    saveas(gcf,[fig_dir sp fl_id{1} '_raw.fig']);
    saveas(gcf,[fig_dir sp fl_id{1} '_raw.png']);
end

%% ------------ Internal stim ------------------
detect_ch = 32;         % Specify a a good chnnel that wil be used to detect the onsets of the stimulation
bad_ch = [21];      % specify bad channels
stim_mode = 'bipolar'; %'monopolar'; %'bipolar'

stim = states.ImplantStimulation==1;       % Extract stim periods
bsline = mean(signal(states.StimulusCode==0,:),1);          % Calculate baseline from preiods without stimulation
stim_codes = unique(states.StimulusCode);  % Extract stimulus codes
stim_codes(stim_codes==0) = [];            % Get rid off baseline stimulus code

stim_duration = (5*parameters.StimulationPulses.NumericValue(3)+...        % 1x leading pulse + 4x charge balancing pulse
    2*parameters.StimulationPulses.NumericValue(4)+...                     % 2x deadzone 0 + 1x deadzone 1
    parameters.StimulationPulses.NumericValue(5))/1e3; 

results.DeclaredDur = zeros(parameters.NumberOfSequences.NumericValue*length(stim_codes),1);                   % initiate array for avg duration of the Implant stimulation state code (in s)
results.ActualDur = zeros(parameters.NumberOfSequences.NumericValue*length(stim_codes),1);                 % initiate array for avg actual stimuluation duration

for stm = 1:length(stim_codes)            % Itterate through the stimuli
    stm_code = stim_codes(stm);           % Variable for stimulus code
    stim_source = parameters.StimulationTriggers.NumericValue(3,stm_code);   % Get stimulation source channel
    stim_dest = parameters.StimulationTriggers.NumericValue(4,stm_code);     % Get stimulation destination channel

    results.(['StimCode' num2str(stm,'%02d')]).PeakAmpl = zeros(size(signal,2),parameters.NumberOfSequences.NumericValue);       % initiate array for avg positive peak amplitude
    % results.(['StimCode' num2str(stm,'%02d')]).NegPeakAmpl = zeros(size(signal,2),parameters.NumberOfSequences.NumericValue);       % initiate array for avg negative peak amplitude

    period_start = find(diff(stim & states.StimulusCode==stm_code)==1)+1; % Find the start of the stimulation period
    period_end = find(diff(stim & states.StimulusCode==stm_code)==-1);    % Find the end of the stimulation period

    if length(period_start)>parameters.NumberOfSequences.NumericValue || length(period_start)~=length(period_end)
        error('StimulusCode %d has inchorent number of stimulation onsets', stm_code)
    end

    for rep = 1:parameters.NumberOfSequences.NumericValue
        stim_period = signal(period_start(rep):period_end(rep),:);        % Extract periods of stimulation during stimulus code
        results.DeclaredDur(stm*rep) = length(stim_period)/FS;           % Get the duration of the Implant stimulation code
        % [start,bipol,stop] = GetStimulusOnsets(stim_period(:,detect_ch),stim_mode);  % Get stimulation onsets
        % results.ActualDur(stm*rep) = (stop-start)/FS;               % Store the duration of the stimulation

        % figure
        % plot(stim_period(:,6)); hold on; plot([start,bipol,stop],stim_period([start,bipol,stop],6),'o','MarkerSize',12)

        for ch = 1:size(stim_period,2)
            [m, idx] = max(abs(stim_period(:,ch)));        % Find the biggest peak in the data
            results.(['StimCode' num2str(stm,'%02d')]).PeakAmpl(ch,rep) = stim_period(idx,ch)-bsline(ch);
        end
    end
    
    avg_Peak = median(results.(['StimCode' num2str(stm,'%02d')]).PeakAmpl,2);

    fg1 = OverlayResults(avg_Peak, ref_ch, stim_source, stim_dest, bad_ch,[fl_id{1},' Peak amplitude distribution']);
    if SAVE_FIG
        % saveas(gcf,[fig_dir sp fl_id{1} '_ch' num2str(source_ch) '-' num2str(dest_ch) '.fig']);
        saveas(fg1,[fig_dir sp fl_id{1} '_StimCode' num2str(stm,'%02d') '.png']);
    end
end

save([fl_pth 'results'],"results");


%% ------ External stim ------
bad_ch = 21; % specify bad channels

% ~~~~~ SPECIFY PARAMETERS FOR CALCULATION OF SPECTRA ~~~~~~~
win_fxn = hann(1*FS); % windowing function
noverlap = floor(length(win_fxn)/2); %overlap between calculations
freqs = 1:500;       % frequencies included in spectra calculation
% ~~~~~ Calculate spectra ~~~~~~~~
[mps,spect_f] = pwelch(data,win_fxn,noverlap,freqs,FS);  % mean power spectrum across whole experiment

pk_80Hz = max(mps(75:85,:),[],1);  % find the higest value of the spectra between 75 and 85 Hz

fg2 = OverlayResults(pk_80Hz, ref_ch, [], [], bad_ch,[fl_id{1},' 80Hz Peak amplitude distribution']);
if SAVE_FIG
    % saveas(fg2,[fig_dir sp fl_id{1} '_ch' num2str(source_ch) '-' num2str(dest_ch) '.fig']);
    saveas(fg2,[fig_dir sp fl_id{1} '.png']);
end



%% Functions

function f = OverlayResults(results,reference,source,destination,badChannels,title_str)
radius = 45;
alph = 0.5;
color_range = 256;
cmap = jet(color_range);

setup = imread('/Users/lampert.frederik/Documents/Mayo/CorTec/Figures/saline-test/implant04/setup.png');
[height,width,depth] = size(setup);

xcords = [2073; 2056; 2022; 2017; 1983; 1964;... % Channels 1-6
    2280; 2306; 2321; 2340; 2357; 2379; ... % Channels 7-12
    1160; 1141; 1121; 1097; 1075; ... % Channels 13-17
    1364; 1383; 1405; 1425; 1451; ... % Channels 18-22
    315; 291; 269; 249; 232; ... % Channels 23-27
    517; 536; 558; 577; 606; ... % Channels 28-32
    2707    ];   % GND channel

ycords = [1406; 1224; 1044;861;684;500;...
    1398; 1221; 1049; 861; 681; 506; ...
    1526; 1354; 1177; 989; 806; ...
    1524; 1346; 1169; 989; 809; ...
    1406; 1226; 1044; 864; 686; ...
    1401; 1224; 1046; 861; 686; ...
    1164    ];

if destination==0   
    destination = 33;  % change the index of the destination to gnd channel
end

if ~isempty(badChannels)
    results(badChannels) = NaN;        % Exlude bad channnels from color ploting
end

results([source,destination])=NaN;      % Exclude the stimulation channels from the color ploting
results(length(xcords)) = NaN; % if the results vector is shorter make it he same size as xcords

X_norm = (results - min(results)) / (max(results) - min(results));
f = figure('Position', [100 100 1500 800])
imshow(setup); hold on
% Plot circles
for i = 1:length(results)
    if i == source
        circles = rectangle('Position', [xcords(i) - radius, ycords(i) - radius, 2*radius, 2*radius], ...
            'EdgeColor', [0.6350 0.0780 0.1840], 'FaceColor', [0 0 0 0.1], ...
            'LineWidth', 6,'Curvature', [1, 1]);
    elseif i == destination
        circles = rectangle('Position', [xcords(i) - radius, ycords(i) - radius, 2*radius, 2*radius], ...
            'EdgeColor', [0.6350 0.0780 0.1840], 'FaceColor', [0 0 0 0.1], ...
            'LineWidth', 6, 'LineStyle', '-.','Curvature', [1, 1]);
    elseif ismember(i,badChannels) || i==33
        circles = rectangle('Position', [xcords(i) - radius, ycords(i) - radius, 2*radius, 2*radius], ...
            'EdgeColor', [0 0 0], 'FaceColor', [0 0 0 0.1], ...
            'LineWidth', 2,'Curvature', [1, 1]);

    else
        % Map the normalized value to the colormap
        color_idx = round(X_norm(i) * (color_range-1)) + 1; % Get index for colormap
        color = cmap(color_idx, :); % Get the color
        % Plot the circle
        if i==reference
            circles = rectangle('Position', [xcords(i) - radius, ycords(i) - radius, 2*radius, 2*radius], ...
                'Curvature', [1, 1], ...
                'EdgeColor', '#77AC30', ...
                'FaceColor', [color alph], ...
                'LineWidth', 6);
        else
            circles = rectangle('Position', [xcords(i) - radius, ycords(i) - radius, 2*radius, 2*radius], ...
                'Curvature', [1, 1], ...
                'EdgeColor', color, ...
                'FaceColor', [color alph], ...
                'LineWidth', 2);
        end
    end
end
legend()
% Display colorbar for reference
colormap(cmap); colorbar;
caxis([min(results), max(results)]);
hold off;
title(title_str, 'Interpreter','none','FontSize',20)
end

function [start,center,stop] = GetStimulusOnsets(signal, mode)
    signal = abs(diff(signal));

    switch mode
        case 'monopolar'
            prmnc = 0.25*max(signal);
            control = true;
            hlp = 1;
            while control
                if hlp == 10
                    disp('Stimulation onsets could not be detected')
                    start = NaN; 
                    stop = NaN;
                    center = NaN;
                    return
                    % figure
                    % plot(signal);hold on
                    % plot(locs,pks,'o','MarkerSize',12)
                end
                [pks,locs] = findpeaks(signal,'MinPeakProminence', prmnc,'SortStr','descend');
                if length(locs)>=3 && (locs(3)>locs(1)&&locs(3)>locs(2)) % If enough peaks were found and they are oordered in a righ way return
                    locs = sort(locs,'ascend');
                    start = locs(1);
                    center = locs(2);
                    stop = locs(3);
                    control = false;
                elseif length(locs)>=3 && locs(2)>locs(1)
                    center = locs(1);
                    if length(signal(1:center-1)) < 3
                        start = 1;
                    else
                        start = findpeaks(signal(1:center-1),"NPeaks",1,"SortStr",'descend');
                    end
                    stop = findpeaks(signal(center+1:end),"NPeaks",1,"SortStr",'descend');
                    control = false;
                else
                    hlp = hlp+1;
                    prmnc = 0.75*prmnc;
                    continue
                end

            end

        case 'bipolar'
            [pk, center, w] = findpeaks(signal,'MinPeakProminence', 0.5*max(signal),'Npeaks',1);
            w = round(w);
            if center-w < 1
                start = 1;
            elseif center+w > length(signal)
                stop = length(signal);
            else
                bsln = mean(signal([1:center-w, center+w:end]));
                start = find(flip(signal(1:center))<=bsln,1,'first'); % find a first value below threshold to the left of the peak
                stop = find(signal(center:end)<=bsln,1,'first');    % find a first value below threshold to the right of the pak 
            end
            
            if isempty(start)
                start = 1; 
            elseif isempty(stop)
                stop = length(signal);
            else
                start = center-start;
                stop = center+stop;
            end
            figure; plot(signal); hold on; plot([start, center, stop],signal([start, center, stop]),'o','MarkerSize',12);
            line([1 length(signal)], [bsln bsln], 'Color', 'r', 'LineStyle', '--')
    end
end
