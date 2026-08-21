%% ================== USER PARAMETERS ==================
AMP_F = {'A57dB','A51dB','A45dB','A39dB'};   % amplification labels
insertedVpp = 2e-3;                     % inserted Vpp [V]
R = 989.2;                              % resistor value [Ohm]

transientTime = 6.999;                 % seconds to discard (steady-state)
minDuration   = 3;                     % minimum duration after trimming [s]

dataPath = 'path/to/data/TF_response';

%% ================== INITIALIZATION ==================
dataStruct = struct();

%% ================== FIND TF FOLDERS ==================
subdirs = dir(dataPath);
subdirs = subdirs([subdirs.isdir]);
subdirs = subdirs(~ismember({subdirs.name},{'.','..'}));

TF_folders = subdirs(startsWith({subdirs.name}, 'TF_response'));

%% ================== MAIN LOOP ==================
for i = 1:numel(TF_folders)

    fldrName = TF_folders(i).name;
    fldrPath = fullfile(dataPath, fldrName);

    % -------- Extract numeric frequency from folder name --------
    % Extract numeric frequency from folder name
    tokens = regexp(fldrName, '_([\d\.]+)Hz', 'tokens');
    
    if isempty(tokens)
        error('Could not extract frequency from folder name: %s', fldrName);
    end
    
    freqStr = tokens{1}{1};        % e.g. '100', '57.5'
    freqHz  = str2double(freqStr); % numeric value

    % -------- List BCI2000 files --------
    fls = dir(fullfile(fldrPath, '*.dat'));
    if isempty(fls)
        warning('No .dat files found in %s', fldrPath);
        continue
    end

    for f = 1:numel(fls)

        filePath = fullfile(fls(f).folder, fls(f).name);

        % -------- Load data --------
        [signal, states, parameters, total_samples ] ...      
        = load_bcidat(filePath, '-calibrated');         % load using MEX package   
        signal = double(signal);

        fs     = parameters.SamplingRate.NumericValue;
        ref_ch = parameters.ReferenceCh.NumericValue;
        ampIdx = parameters.AmplificationFactor.NumericValue + 1;
        amp_f  = AMP_F{ampIdx};

        % % -------- Reference subtraction --------
        % signal = signal - signal(:, ref_ch);
       
        % -------- Detrend entire signal --------
        % Removes slow drift / offset across channels
        % signal = detrend(signal);

        % -------- Remove transient --------
        n0 = round(transientTime * fs) + 1;

        if size(signal,1) - n0 < minDuration * fs
            warning('Signal too short after trimming: %s', fls(f).name);
            continue
        end

        % Remove ref channel and set up channel numbers
        nCh = 1:size(signal,2);                     % Channel numbers
        signal(:, ref_ch) = [];                     % Remove reference channel before analysis
        nCh(ref_ch) = [];                           % Remove also from channel Numbering
        
        % Clip signals
        signal = signal(n0:end, :);                 % 

        % -------- Peak-to-peak extraction --------
        % Per-channel Vpp
        chVpp = peak2peak(signal, 1);

        % Median across channels
        Vpp = median(chVpp, 'omitnan');

        % -------- Store results --------
        if ~isfield(dataStruct, amp_f)
            % Build variable names: Frequency, Median Vpp, Ch1 ... ChN
            varNames = [{'Frequency_Hz','Vpp_muV'}, ...
                arrayfun(@(ch) sprintf('Ch%d', ch), nCh, ...
                'UniformOutput', false)];

            dataStruct.(amp_f) = table( ...
                'Size', [0 numel(varNames)], ...
                'VariableTypes', repmat({'double'}, 1, numel(varNames)), ...
                'VariableNames', varNames);
        end

        % -------- Add one complete row at once --------
        n = height(dataStruct.(amp_f)) + 1;

        dataStruct.(amp_f){n, :} = [freqHz, Vpp, chVpp];

    end
end

dataStruct.MeasurementInfo.InsertedVpp = insertedVpp;
dataStruct.MeasurementInfo.R = R;
dataStruct.MeasurementInfo.TransientRemoved_s = transientTime;

%% ================== SORT RESULTS ==================
ampNames = fieldnames(dataStruct);
ampNames(strcmp(ampNames,'MeasurementInfo')) = [];

for a = 1:numel(ampNames)
    amp = ampNames{a};
    dataStruct.(amp) = sortrows(dataStruct.(amp), 'Frequency_Hz');
end

%% ================== SAVE RESULTS ==================
outFile = fullfile(dataPath, 'TF_response_Vpp.xlsx');

for a = 1:numel(ampNames)
        amp = ampNames{a};
    writetable(dataStruct.(amp), outFile, 'Sheet', amp);
end

infoT = struct2table(dataStruct.MeasurementInfo);
writetable(infoT, outFile, 'Sheet', 'MeasurementInfo');