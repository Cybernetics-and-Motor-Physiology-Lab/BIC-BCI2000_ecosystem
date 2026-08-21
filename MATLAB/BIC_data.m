classdef BIC_data
    properties
        % Required
        dataRoot               % Root folder path containing recordings.mat, etc. 

        % Dependency tracking
        DependencyReport = struct()   % Stores dependency check results

        % Optional subject-specific and general CorTec metadata
        subjectName = ''       % Subject identifier (string)
        previousDataRoot = ''  % Previous data root
        implantDate = []       % datetime or empty
        excludedDates = []     % Vector of datetimes
        leadAssignment = dictionary; % Dictionary with the name of the lead and struct containing the iformation
        rerefMethods = {'Native','CAR','BP','leadCAR','leadSpecific'};  % Default rereferencing methods
        rerefSpecs = struct('construct','','specification','','ReRefCh',[]);
        amplFactors = {'57.5dB','51.5dB','45.5dB','39.5dB'};      % Default Amplification factors
        LongtermStimFolders = {};
        SpikeDetectionParams = struct('LF',10,'HF',85, 'th_mult', 3.65, 'Amp_th_min',85);
        BadChannelThresh = 3;
        spectralParameters = struct('window',2 , ...
            'overlapCoeff',0.5, ...
            'windowFunction','hann', ...
            'NFFT',2000, ...
            'bands', struct('Delta', [1 4], 'Theta', [4 8], 'Alpha', [8 15], 'Beta',  [15 30], 'Gamma', [30 100]));
        knownBadChannels = [];

        % Working 
        newRecordings = struct([])   % Holds newly found .dat files
        fileStatus = cell2table(cell(0,4), 'VariableNames', {'FileName', 'InRecordings', 'Discarded', 'PSD'});
        subjectSummary = struct([])
        longtermStimFeatures = struct([])
        longtermStimSummary = struct([])
    end

    properties (Dependent) % File/folder paths (auto-initialized)      
        recordingsFile
        recordingsMATFolder
        impedanceFile
        longtermStimFile
        longtermStimFolder
        discardedFile
        spikeEventsFile
    end

    properties (SetAccess = private)
        recordings = struct([])            % Loaded from recordings.mat
        impedanceData = struct([]) 
        longtermStimData = struct([]) 
        discardedRecordings = struct([]);  % Struct array of skipped/corrupted recordings
        samplingRate = 1000    % Predefined sampling rate
    end

    methods
        %% ~~~~~~~~~~~~~~~ Class initializaion function ~~~~~~~~~~~~~~~~~~~
        function obj = BIC_data(dataRoot)
            % Constructor - requires only the root data directory
            arguments
                dataRoot (1, :) char
            end
            obj.dataRoot = dataRoot;

            % Check core dependencies
            [obj, ~] = checkDependencies(obj,'all', 'warning');

            % % Set default subpaths
            % obj.recordingsFile      = fullfile(dataRoot, 'recordings.mat');
            % obj.recordingsMATFolder = fullfile(dataRoot, 'recordings_MAT');
            % obj.impedanceFile       = fullfile(dataRoot, 'impedance.mat');
            % obj.longtermStimFile    = fullfile(dataRoot, 'longterm_stim.mat');
            % obj.longtermStimFolder  = fullfile(dataRoot, 'longterm_stim_MAT');
            % obj.discardedFile = fullfile(dataRoot, 'discarded.mat');

            % Create recordingsMATFolder if it doesn't exist
            if ~exist(obj.recordingsMATFolder, 'dir')
                mkdir(obj.recordingsMATFolder);
            end
            if ~exist(obj.longtermStimFolder, 'dir')
                mkdir(obj.longtermStimFolder);
            end

            % Auto-load metadata if files exist
            if isfile(obj.recordingsFile)
                tmp = load(obj.recordingsFile);
                obj.recordings = tmp.dataFiles;
            end
            if isfile(obj.impedanceFile)
                tmp = load(obj.impedanceFile);
                obj.impedanceData = tmp.impedance;
            end
            if isfile(obj.longtermStimFile)
                tmp = load(obj.longtermStimFile);
                obj.longtermStimData = tmp.dataFiles;
            end
            if isfile(obj.discardedFile)
                tmp = load(obj.discardedFile);
                obj.discardedRecordings = tmp.discarded;
            end
        end

        function val = get.recordingsFile(obj)
            val = fullfile(obj.dataRoot, 'recordings.mat');
        end
        
        function val = get.recordingsMATFolder(obj)
            val = fullfile(obj.dataRoot, 'recordings_MAT');
        end

        function val = get.impedanceFile(obj)
            val = fullfile(obj.dataRoot, 'impedance.mat');
        end

        function val = get.longtermStimFile(obj)
            val = fullfile(obj.dataRoot, 'longterm_stim.mat');
        end
        
        function val = get.longtermStimFolder(obj)
            val = fullfile(obj.dataRoot, 'longterm_stim_MAT');
        end

        function val = get.discardedFile(obj)
            val = fullfile(obj.dataRoot, 'discarded.mat');
        end

        function val = get.spikeEventsFile(obj)
            val = fullfile(obj.dataRoot, 'spike_events.mat');
        end

        %% ~~~~~~~~~~~~~~~ Data search and tracking functions ~~~~~~~~~~~~~
        function obj = synchronizeFileStatus(obj)
            % synchronizeFileStatus - updates fileStatus based on current recordings,
            % discarded files, and existing PSD .mat files

            disp('Synchronizing file status...');

            % Step 1: Find all .dat files in dataRoot
            allFiles = dir(fullfile(obj.dataRoot, '**', '*.dat'));
            allFiles = rmfield(allFiles, 'isdir');
            allNames = {allFiles.name}';

            % Step 2: Get reference names
            recNames = {};
            discNames = {};
            spikeDetected = false(size(allNames));

            if ~isempty(obj.recordings)
                recNames = {obj.recordings.name};
                for i = 1:numel(obj.recordings)
                    refs = fieldnames(obj.recordings(i));
                    refHasSpikes = any(contains(refs, obj.rerefMethods) & ...
                        structfun(@(r) isfield(r, 'SpikeDetect'), obj.recordings(i)));
                    if refHasSpikes
                        nameIdx = strcmp(allNames, obj.recordings(i).name);
                        spikeDetected(nameIdx) = true;
                    end
                end
            end
            if ~isempty(obj.discardedRecordings)
                discNames = {obj.discardedRecordings.name};
            end

            % Step 3: Check .mat files in recordings_MAT
            matFiles = dir(fullfile(obj.recordingsMATFolder, '*.mat'));
            matNames = erase({matFiles.name}, '.mat');
            matNames = strcat(matNames, '.dat');

            % Step 4: Build the table
            inRec = ismember(allNames, recNames);
            inDisc = ismember(allNames, discNames);
            hasPSD = ismember(allNames, matNames);
            hasSpike = spikeDetected;

            obj.fileStatus = table(allNames, inRec, inDisc, hasPSD, hasSpike, ...
                'VariableNames', {'FileName', 'InRecordings', 'Discarded', 'PSD','SpikeDetect'});

            disp(['Total files tracked: ', num2str(height(obj.fileStatus))]);
        end

        function obj = searchForNewData(obj)
            % searchForNewData - Scans dataRoot for new .dat files and updates obj.fileStatus
            % New, non-discarded files are stored in obj.newRecordings

            disp(['Searching for .dat files in ', obj.dataRoot]);

            % Always re-sync fileStatus
            obj = obj.synchronizeFileStatus();

            % Extract list of new, valid files
            newMask = ~obj.fileStatus.InRecordings & ...
                ~obj.fileStatus.Discarded & ...
                ~obj.fileStatus.PSD;

            newFileNames = obj.fileStatus.FileName(newMask);

            % Find file structs
            allFiles = dir(fullfile(obj.dataRoot, '**', '*.dat'));
            % Ignore hidden .dat files such as '.name.dat'
            allFiles = allFiles(~startsWith({allFiles.name}, '.'));
            % Remove isdir field for consistency
            allFiles = rmfield(allFiles, 'isdir');

            newFilesIdx = ismember({allFiles.name}, newFileNames);
            obj.newRecordings = allFiles(newFilesIdx);

            disp(['Found ', num2str(numel(obj.newRecordings)), ' new file(s) ready for processing.']);
        end

        function obj = searchForNewImp(obj)
            % searchForNewImp - Scans for impedance log files and parses them
            % Updates obj.impedanceData and optionally filters out corrupted entries.
            disp(['Searching for Impedance files in ', obj.dataRoot]);

            % Define filename patterns
            fls_pattern = {
                '*_impedances*.txt';
                '*_info*.txt';
                'Grd_*.log'
                };

            % Collect all impedance file paths
            ImpFiles = dir(fullfile(obj.dataRoot, '**', fls_pattern{1}));
            for pat = 2:length(fls_pattern)
                ImpFiles = [ImpFiles; dir(fullfile(obj.dataRoot, '**', fls_pattern{pat}))];
            end

            % Remove zero-byte files
            emptyIdx = [ImpFiles.bytes] == 0;
            ImpFiles(emptyIdx) = [];
            ImpFiles = rmfield(ImpFiles, 'isdir');

            % Compare with existing impedanceData (by filename)
            if ~isempty(obj.impedanceData)
                [~, new_idx] = setdiff({ImpFiles.name}, {obj.impedanceData.name});
                ImpFiles = ImpFiles(new_idx);
                disp(['Found ', num2str(length(ImpFiles)), ' new impedance file(s).']);
            else
                disp(['Found ', num2str(length(ImpFiles)), ' impedance file(s).']);
            end

            % Initialize removal index
            emptyIdx = [];

            % Process each file
            for i = 1:length(ImpFiles)
                progressbar(i, length(ImpFiles))  % Console progress
                fullPath = fullfile(ImpFiles(i).folder, ImpFiles(i).name);

                % Attempt to extract impedance
                impValues = extractImpedanceValues(fullPath);

                if isempty(impValues)
                    emptyIdx = [emptyIdx, i];  % Mark as unusable
                    continue;
                end

                % Convert to numeric array
                impArray = cell2mat(impValues(:,2));

                % Pad to 32 channels if short
                if numel(impArray) < 32
                    impArray(end+1:32) = nan;
                end

                % Assign impedance values as fields
                for j = 1:32
                    ch = sprintf('Ch%d', j);
                    ImpFiles(i).(ch) = impArray(j);
                end
            end

            % Remove empty entries
            ImpFiles(emptyIdx) = [];

            % Append to impedanceData
            if ~isempty(ImpFiles) % If new data were found
                if isempty(obj.impedanceData)
                    obj.impedanceData = ImpFiles;
                else
                    obj.impedanceData = [obj.impedanceData; ImpFiles];
                end
            end

        end

        %% ~~~~~~~~~~~~~~~ Data porcessing functions ~~~~~~~~~~~~~~~~~~~~~~
        function obj = processRecordings(obj, reftypes, OVERWRITE)
            % processRecordings - Process .dat files into structured metadata
            % Optional Inputs:
            %   reftypes  - cell array of rereferencing methods (default: obj.rerefMethods)
            %   OVERWRITE - logical (default: false)
            % Output:
            %   Updates obj.recordings

            % Handle optional inputs
            if nargin < 2 || isempty(reftypes)
                reftypes = obj.rerefMethods;
            end
            if nargin < 3
                OVERWRITE = false;
            end

            % Ensure reftypes is always a cell array of char
            if ischar(reftypes) || isstring(reftypes)
                reftypes = {char(reftypes)};  % convert to single-element cell
            end

            % Validate contents of reftypes
            validTypes = {'Native','CAR','BP','leadCAR','leadSpecific'};
            if any(~ismember(reftypes, validTypes))
                error('Invalid reftypes. Must be one or more of: %s', strjoin(validTypes, ', '));
            end

            % Select which files to process
            if OVERWRITE
                dataFiles = obj.recordings;
            elseif ~isempty(obj.newRecordings)
                dataFiles = obj.newRecordings;
            else
                disp('No new recordings found and OVERWRITE is false. Nothing to process.');
                return;
            end

            % Processing parameters 
            min_rec_time = 10;               % seconds
            ampl_f = obj.amplFactors;
            bCh_thr =  obj.BadChannelThresh;               % bad channel threshold
            knownBadChnls = obj.knownBadChannels;
            cnstrct = obj.rerefSpecs.construct;
            spec = obj.rerefSpecs.specification;           % rereferencing spec
            RR_chnls = obj.rerefSpecs.ReRefCh;  % example channel group

            % Initialize
            nFiles = numel(dataFiles);
            discarded = false(nFiles, 1);
            dataFiles = initializeDataStruct(dataFiles,reftypes);

            % Inform the user 
            fprintf('Processing %d files\n', nFiles);

            % Process files
            parfor i = 1:nFiles
                try
                    if contains(dataFiles(i).name, ' copy')
                        discarded(i) = true;
                        continue;
                    end

                    % Change filepaths if needed
                    if contains(dataFiles(i).folder, obj.dataRoot)
                        fl_pth = fullfile(dataFiles(i).folder, dataFiles(i).name);
                    else
                        fl_pth = fullfile(...
                            strrep(dataFiles(i).folder,obj.previousDataRoot, obj.dataRoot),...
                            dataFiles(i).name);
                    end

                    % Load recording, extract recording info & parameters
                    [signal, states, parameters, N_samp] = load_bcidat(fl_pth, '-calibrated');

                    % Discard corrupt/short recordings & Interpolate Lost samples
                    FS = parameters.SamplingRate.NumericValue;
                    if sum(isnan(signal),'all') > numel(signal)/2 || length(signal) < min_rec_time*FS
                        discarded(i) = true;
                        continue;
                    else
                        signal = interpolate_lost_samples(double(signal), states,'linear');
                    end

                    % Metadata
                    dataFiles(i).SamplingRate = FS;
                    dataFiles(i).NativeRef = parameters.ReferenceCh.NumericValue;
                    dataFiles(i).Duration = string(seconds(N_samp/FS),"hh:mm:ss.SSS");
                    dataFiles(i).PacketLoss = calculatePacketLoss(states);

                    try
                        dataFiles(i).AmplFactor = ampl_f{parameters.AmplificationFactor.NumericValue+1};
                        dataFiles(i).Ground = parameters.UseGround.NumericValue;
                    catch
                        dataFiles(i).AmplFactor = '57.5dB';
                        dataFiles(i).Ground = 0;
                    end

                    try
                        dataFiles(i).Impedance = extractImpedanceValues(strrep(fl_pth, '.dat', '_info0.txt'));
                    catch
                        dataFiles(i).Impedance = {};
                    end

                    % Stimulation info
                    sc_dur = [];
                    sc = [];
                    sc_onsets = [];
                    if parameters.EnableStimulation.NumericValue == 1
                        [sc_dur, sc] = runlength(states.ImplantStimulation, N_samp);
                        sc_onsets = [1; cumsum(sc_dur(1:end-1)) + 1];
                        dataFiles(i).StimInfo.StimPulses = parseBCI2000_parameters_values(parameters, 'StimulationPulses');
                        dataFiles(i).StimInfo.StimTriggers = parseBCI2000_parameters_values(parameters, 'StimulationTriggers');
                        dataFiles(i).StimInfo.StimSequence.Sequence = sc;
                        dataFiles(i).StimInfo.StimSequence.SequenceOnset = sc_onsets;
                        dataFiles(i).StimInfo.StimSequence.SequenceDuration = sc_dur;
                    end

                    % rest_indices = {};

                    % Detect Bad channels
                    if parameters.EnableStimulation.NumericValue == 1 && ~isempty(sc_dur)
                        rest_idx = [sc_onsets(sc==0), sc_onsets(sc==0)+sc_dur(sc==0)-1];
                        rest_indices = arrayfun(@(s, e) s:e, rest_idx(:,1), rest_idx(:,2), 'UniformOutput', false);
                        [bad_chnls, badChannelInfo] = identifyBadChannels(signal([rest_indices{:}], :), cnstrct, FS,bCh_thr, knownBadChnls, dataFiles(i).Impedance, 0);
                    else 
                        [bad_chnls, badChannelInfo] = identifyBadChannels(signal, cnstrct, FS,bCh_thr, knownBadChnls, dataFiles(i).Impedance, 0);
                    end
                    dataFiles(i).BadChannels.BadChannels = bad_chnls;
                    dataFiles(i).BadChannels.Info = badChannelInfo;

                    % Rereferencing + stats
                    for ridx = 1:numel(reftypes)
                        ref = reftypes{ridx};

                        try
                            % Rereference
                            switch ref
                                case 'Native'
                                    ref_signal = signal;
                                    ch_names = arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(signal,2), 'UniformOutput', false);
                                    ch_names(bad_chnls) = arrayfun(@(n) sprintf('BAD-Ch %d', n), bad_chnls, 'UniformOutput', false);
                                case 'BP'
                                    [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, spec, true);
                                    dataFiles(i).(ref).RefSpec = spec;
                                case 'leadSpecific'
                                    try
                                        [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, RR_chnls, true);
                                        dataFiles(i).(ref).RefSpec = RR_chnls;
                                    catch
                                        continue;
                                    end
                                otherwise
                                    [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, [], true);
                            end

                            dataFiles(i).(ref).ChannelNames = ch_names;
                            stats = struct();
                            stats.Average = mean(ref_signal, 1);
                            stats.Median = median(ref_signal, 1);
                            stats.STD = std(ref_signal, 1);
                            stats.Variance = var(ref_signal, 1);
                            stats.RMS = rms(ref_signal, 1);
                            stats.Skewness = skewness(ref_signal, 1);
                            stats.Kurtosis = kurtosis(ref_signal, 1);
                            stats.ZeroCrossRate = zerocrossrate(ref_signal);
                            stats.AutoCorr = corr(ref_signal);
                            stats.AutoCov = cov(ref_signal);
                            dataFiles(i).(ref).Stats = stats;

                        catch ME_ref
                            dataFiles(i).(ref).ProcessingErrorMsg = ME_ref.message;
                            dataFiles(i).(ref).ProcessingFailed = true;
                        end
                    end
                catch ME
                    discarded(i) = true;
                    dataFiles(i).ProcessingErrorMsg = ME;
                end
                progressbar(i, nFiles)  % Console progress
            end

            % Clean and assign results
            baseFieldnames = {'name', 'folder', 'date', 'bytes', 'datenum'};
            badFiles = dataFiles(discarded);
            badFiles = rmfield(badFiles,setdiff(fieldnames(dataFiles),baseFieldnames));
            obj.discardedRecordings = [obj.discardedRecordings; badFiles];
            % Append to object properties
            dataFiles(discarded) = [];
            obj.newRecordings = struct([]); % clean the new recordings
            if OVERWRITE
                obj.recordings = dataFiles;
            else
                obj.recordings = [obj.recordings; dataFiles];
            end
        end

        function obj = processLongtermStim(obj, rereftypes, OVERWRITE)
            % processLongtermStim - Processes long-term stim recordings and computes STFTs
            % Inputs (optional):
            %   rereftypes - rereferencing types to use (default = obj.rerefMethods)
            %   OVERWRITE  - whether to overwrite existing .mat files (default = false)

            if nargin < 2 || isempty(rereftypes)
                rereftypes = obj.rerefMethods;
            end
            if nargin < 3
                OVERWRITE = false;
            end

            % Step 1: Validate all folder paths
            resolvedFolders = cell(size(obj.LongtermStimFolders));
            for i = 1:numel(obj.LongtermStimFolders)
                folder = obj.LongtermStimFolders{i};
                if isfolder(folder)
                    resolvedFolders{i} = string(folder);
                elseif isfolder(fullfile(obj.dataRoot, folder))
                    resolvedFolders{i} = string(fullfile(obj.dataRoot, folder));
                else
                    error('Folder not found: "%s"', folder);
                end
            end
            resolvedFolders = string(resolvedFolders);  % Convert cell to string array

            % Step 2: Filter recordings in longterm folders
            % isLongterm = arrayfun(@(r) any(startsWith(r.folder, resolvedFolders)), obj.recordings);
            isLongterm = arrayfun(@(r) contains(r.folder, resolvedFolders), obj.recordings);
            dataFile = obj.recordings(isLongterm);

            % Step 3: Check existing .mat files in folders
            allMatFiles = [];
            for i = 1:numel(resolvedFolders)
                allMatFiles = [allMatFiles; dir(fullfile(resolvedFolders{i}, '*.mat'))];
            end
            existingMatNames = erase({allMatFiles.name}, '.mat');
            existingMatNames = strcat(existingMatNames, '.dat');

            % Step 4: Filter recordings based on overwrite
            if OVERWRITE
                toProcess = dataFile;
            else
                notProcessedMask = ~ismember({dataFile.name}, existingMatNames);
                toProcess = dataFile(notProcessedMask);
            end

            if isempty(toProcess)
                disp('No new long-term stim files to process.');
                return;
            end

            % Re-referencing parameters
            cnstrct = obj.rerefSpecs.construct;
            spec = obj.rerefSpecs.specification;
            RR_chnls = obj.rerefSpecs.ReRefCh;

            % Specify spike detecion parameters
            LF =  obj.SpikeDetectionParams.LF;
            HF = obj.SpikeDetectionParams.HF;
            th_mult = obj.SpikeDetectionParams.th_mult;
            minAmp = obj.SpikeDetectionParams.Amp_th_min;

            % Initialize STFT parameters
            FS = obj.samplingRate;
            wnd = obj.spectralParameters.window;
            noverlap_coeff = obj.spectralParameters.overlapCoeff;
            win_fxn = obj.spectralParameters.windowFunction;
            nfft = obj.spectralParameters.NFFT;

            switch lower(win_fxn)
                case 'hann'
                    wndw_fxn = hann(wnd * FS, 'periodic');
                case 'hamming'
                    wndw_fxn = hamming(wnd * FS, 'periodic');
                otherwise
                    error('Unsupported window function: %s', win_fxn);
            end
            noverlap = floor(length(wndw_fxn) * noverlap_coeff);

            %Specify outpudir
            outputDir = obj.longtermStimFolder;
            if ~isfolder(outputDir)
                mkdir(outputDir);
            end

            % Initialize structure for detected spike events
            spikeExamples = rmfield(toProcess, ...
                setdiff(fieldnames(toProcess), ...
                {'name', 'folder','date'}));

            % Start processing
            nFiles = numel(toProcess);
            fprintf('Processing %d long-term stim file(s)',nFiles);

            for f = 1:nFiles
                progressbar(f, nFiles);
                recIdx = find(strcmp({obj.recordings.name}, toProcess(f).name), 1);
                try
                    % Change filepaths if needed
                    if contains(toProcess(f).folder, obj.dataRoot)
                        fl_pth = fullfile(toProcess(f).folder, toProcess(f).name);
                    else
                        fl_pth = fullfile(...
                            strrep(toProcess(f).folder,obj.previousDataRoot, obj.dataRoot),...
                            toProcess(f).name);
                    end

                    [signal, states, ~, ~] = load_bcidat(fl_pth, '-calibrated');
                    signal = double(signal);
                    signal = interpolate_lost_samples(signal, states, 'linear');

                    bad_chnls = toProcess(f).BadChannels.BadChannels;
                    fileStruct = struct('sourceFile',fullfile(toProcess(f).folder, toProcess(f).name));

                    % Select resting periods
                    sc = toProcess(f).StimInfo.StimSequence.Sequence;
                    sc_onsets = toProcess(f).StimInfo.StimSequence.SequenceOnset;
                    sc_dur = toProcess(f).StimInfo.StimSequence.SequenceDuration;
                    rest_idx = [sc_onsets(sc==0), sc_onsets(sc==0)+sc_dur(sc==0)-1];
                    % Select stimulation periods
                    stim_idx = [sc_onsets(sc~=0), sc_onsets(sc~=0)+sc_dur(sc~=0)-1];

                    for r = 1:length(rereftypes)
                        ref = rereftypes{r};
                        try
                            switch ref
                                case 'Native'
                                    ref_signal = signal;
                                    ch_names = arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(signal,2), 'UniformOutput', false);
                                    ch_names(bad_chnls) = arrayfun(@(n) sprintf('BAD-Ch %d', n), bad_chnls, 'UniformOutput', false);
                                case 'BP'
                                    [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, spec, true);
                                case 'leadSpecific'
                                    try
                                        [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, RR_chnls, true);
                                    catch
                                        continue;
                                    end
                                otherwise
                                    [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, [], true);
                            end

                            % Initialize structs
                            spectral_features = struct();
                            temporal_features = struct();
                            rest_sig = [];
                            spect_f = [];   % store frequency vector once computed
                            % temp_struct = struct('stat',[],'timestamp',[],'Channelinformation',[],'Denosing',[],'RFOMP',[]);

                            for s = 1:size(rest_idx,1)
                                % --------------- Extract rest periods ---------
                                rest_sig = ref_signal(rest_idx(s,1):rest_idx(s,2),:);
                                trialName = sprintf('Trial%d', s);

                                % --------------- Temporal features ----
                                if size(rest_sig,1) < 2
                                    temporal_features.(trialName).Covariance = [];
                                    temporal_features.(trialName).CrossCorrelation = [];
                                    temporal_features.(trialName).MutualInformation = [];
                                else
                                    temporal_features.(trialName).Covariance = single(cov(rest_sig, 'omitrows'));
                                    temporal_features.(trialName).CrossCorrelation = single(corrcoef(rest_sig, 'Rows','pairwise'));
                                    temporal_features.(trialName).MutualInformation = [];
                                end

                                % --------------- Spectral features -----
                                if size(rest_sig,1) < length(wndw_fxn)
                                    spectral_features.(trialName).PSD = [];
                                    spectral_features.(trialName).Coherence = [];
                                else
                                    % Caluclate and extract spectral features
                                    spectralFeatures = extractSpectralFeatures(rest_sig, FS, wndw_fxn, noverlap, nfft, obj.spectralParameters.bands);
                                    % Convert to singles and save
                                    spectral_features.(trialName) = struct2single(spectralFeatures);
                                end

                                % -------------------- Detect spikes ------------
                                % Upsample signal to 2k Hz
                                q = 2;
                                rest_sig = resample(rest_sig,q,1,'Dimension',1);
                                %Create montage  for spike elimination
                                montage = struct('SampleRate', q*FS, 'ChannelNames', {ch_names});
                                % For spike detection and elimination functions please refer to 
                                % Ayyoubi, A. H., Besheli, B. F., Swamy, C. P., Okkabaz, J. L., Miller, K. J., Worrell, G. A., & Ince, N. F. (2025). 
                                % Spurious Spike Elimination using Sparse Signal Processing Improves Seizure Onset Zone Delineation in Brief Intraoperative iEEG Recordings. 
                                % The ... Midwest Symposium on Circuits and Systems conference proceedings : MWSCAS. Midwest Symposium on Circuits and Systems, 2025, 666–670. 
                                % https://doi.org/10.1109/mwscas53549.2025.11244470
                                %Detect spikes
                                IIS_Det = Spike_detection_function(rest_sig, montage, [], 512, [LF, HF], minAmp, th_mult);
                                % Filter spikes 
                                Spike_RFOMP =  pSpike_elimination_RFOMP(IIS_Det.spike,montage.SampleRate);

                                % --- Extract spike examples (up to 10 randomly selected that passed the denoising) ---
                                denoisedSpikesIdx = find(Spike_RFOMP.feature.Pred==1);  % [samples x nSpikes]
                                nSpikes = size(denoisedSpikesIdx, 1);
                                if nSpikes > 0
                                    nExamples = min(10, nSpikes);
                                    randIdx = randperm(nSpikes, nExamples);
                                    spikeEventSubset = IIS_Det.spike.events.Raw(randIdx,:);  % [samples x nSpikes]
                                    spikeEventSubsetChLabels = ch_names(IIS_Det.spike.Channelinformation(randIdx));
                                    spikeEventTimeStamp = IIS_Det.spike.timestamp(randIdx);
                                else
                                    spikeEventSubset = [];  % no spikes found
                                    spikeEventSubsetChLabels = [];
                                    spikeEventTimeStamp = [];
                                end

                                % Save spike examples and statistic
                                IIS_Det.spike = rmfield(IIS_Det.spike,'events');      % Save spike detection result back to master recording   % Remove the events to save memory

                                obj.recordings(recIdx).(ref).SpikeDetect = IIS_Det.spike;
                                obj.recordings(recIdx).(ref).SpikeDetect.Montage = montage; % !!!! changed from channel names
                                obj.recordings(recIdx).(ref).SpikeDetect.RFOMP = Spike_RFOMP; % Save the output of spike filtering

                                % Save to spikeExamples struct
                                spikeExamples(i).(ref).spikes = spikeEventSubset;
                                spikeExamples(i).(ref).channelLabels = spikeEventSubsetChLabels;
                                spikeExamples(i).(ref).timeStamps = spikeEventTimeStamp;
                            end

                            stim_ps = cell(1,size(stim_idx,1));
                            stim_sig = [];
                            for s = 1:size(stim_idx,1)
                                stim_sig = single(ref_signal(stim_idx(s,1):stim_idx(s,2),:));
                                if length(stim_sig)<length(wndw_fxn)
                                    stim_ps{s} = [];
                                else
                                    stim_ps{s} = pwelch(stim_sig,wndw_fxn,noverlap,nfft,FS); %nfft/n
                                end
                            end

                            % Chat save the adjusted parts into the struct
                            % Save to structured field
                            fileStruct.(ref).StimPSDs = stim_ps;
                            fileStruct.(ref).SpectralFeatures = spectral_features;
                            fileStruct.(ref).TemporalFeatures = temporal_features;
                            fileStruct.(ref).ChannelNames = ch_names;
                            fileStruct.(ref).Info.WindowFunc = win_fxn;
                            fileStruct.(ref).Info.WindowLength = length(wndw_fxn);
                            fileStruct.(ref).Info.Overlap = noverlap;
                            fileStruct.(ref).Info.Frequencies = spect_f;
                            fileStruct.(ref).Info.RestIndices = rest_idx;
                            fileStruct.(ref).Info.StimIndices = stim_idx;
                        catch ME_ref
                            fileStruct.(i).(ref).ProcessingErrorMsg = ME_ref.message;
                            fileStruct.(i).(ref).ProcessingFailed = true;
                        end
                    end

                    % Save to .mat file (one per recording)
                    outName = strrep(toProcess(f).name, '.dat', '.mat');
                    tempFilename = fullfile(outputDir, outName);
                    parsave(tempFilename, fileStruct);
                catch ME
                    obj.recordings(recIdx).ProcessingErrorMsg = ME;
                    % warning('Error processing file %s: %s', toProcess(f).name, ME.message);
                end
            end
            % Append to object properties
            if OVERWRITE
                obj.longtermStimData = toProcess;
            else
                obj.longtermStimData = [obj.longtermStimData; toProcess];
            end

            % Save or append spikeExamples
            if isfile(obj.spikeEventsFile) && ~OVERWRITE
                % Append new entries to existing
                loaded = load(obj.spikeEventsFile, 'spikeExamples');
                spikeExamples = [loaded.spikeExamples; spikeExamples];
                save(obj.spikeEventsFile, 'spikeExamples', '-v7.3');
            else
                save(obj.spikeEventsFile, 'spikeExamples', '-v7.3');
            end
            disp(['Spike event examples saved to: ', obj.spikeEventsFile]);
        end

        function obj = extractSpectralFeatures(obj, reftypes, OVERWRITE)
            % extractSpectralFeatures - Computes PSDs, CPSDs, power in bands and coherence for rereferenced signals
            % Optional Inputs:
            %   reftypes   - rereferencing methods (default: obj.rerefMethods)
            %   OVERWRITE  - whether to recompute for all files (default: false)

            % Handle optional args
            if nargin < 2 || isempty(reftypes)
                reftypes = obj.rerefMethods;
            end
            if nargin < 3
                OVERWRITE = false;
            end
            
            psdFiles = obj.fileStatus.FileName(obj.fileStatus.InRecordings & ...
                ~obj.fileStatus.Discarded & ...
                ~obj.fileStatus.PSD);

            % Select recordings to process
            if OVERWRITE
                dataFiles = obj.recordings;
            elseif ~isempty(psdFiles)
                dataFiles = obj.recordings(ismember({obj.recordings.name},psdFiles));
            else
                disp('No new recordings found and OVERWRITE is false. Nothing to process.');
                return;
            end

            % Ensure all files are in obj.recordings
            allNames = {dataFiles.name};
            existingNames = {obj.recordings.name};

            [~, missingIdx] = setdiff(allNames, existingNames);
            missingFiles = dataFiles(missingIdx);

            if ~isempty(missingFiles)
                disp('Processing new recordings missing from the recordings structe');
                obj.newRecordings = missingFiles;
                obj = obj.processRecordings(reftypes);  % Call own method
            end

            outDataFolder = obj.recordingsMATFolder;
            if ~isfolder(outDataFolder)
                mkdir(outDataFolder);
            end

            % Prepare full struct aligned with obj.recordings
            recNames = {obj.recordings.name};
            matchedIdx = cellfun(@(n) find(strcmp(recNames, n), 1, 'first'), allNames,'UniformOutput',false);
            tempData = obj.recordings(cell2mat(matchedIdx));

            % Rereference
            cnstrct = obj.rerefSpecs.construct;
            spec = obj.rerefSpecs.specification;
            RR_chnls = obj.rerefSpecs.ReRefCh;

            % Extract spectral parameters from object
            FS = obj.samplingRate;
            wnd = obj.spectralParameters.window;
            noverlap_coeff = obj.spectralParameters.overlapCoeff;
            win_fxn = obj.spectralParameters.windowFunction;
            nfft = obj.spectralParameters.NFFT;
            % Windowing setup
            switch lower(win_fxn)
                case 'hann'
                    wndw_fxn = hann(wnd*FS, 'periodic');
                case 'hamming'
                    wndw_fxn = hamming(wnd*FS, 'periodic');
                otherwise
                    error('Unsupported window function: %s', win_fxn);
            end
            noverlap = floor(length(wndw_fxn) * noverlap_coeff);
            
            % Inform the user 
            disp('Calculating spectra...');

            % Spectral calculation
            nFiles = numel(tempData);
            for f = 1:nFiles
                % Adjust filepath if needed
                if contains(tempData(f).folder, obj.dataRoot)
                    fl_pth = fullfile(tempData(f).folder,tempData(f).name);
                else
                    fl_pth = fullfile(strrep(tempData(f).folder, ...
                        obj.previousDataRoot, obj.dataRoot), ...
                        tempData(f).name);
                end

                bad_chnls = tempData(f).BadChannels.BadChannels;
                spectraStruct = struct('sourceFile',fl_pth);
                try
                    [signal, states, ~, ~] = load_bcidat(fl_pth, '-calibrated');
                    signal = double(signal);

                    % Interpolate
                    signal = interpolate_lost_samples(signal, states,'linear');

                    for ridx = 1:length(reftypes)
                        ref = reftypes{ridx};
                        try
                            switch ref
                                case 'Native'
                                    ref_signal = signal;
                                case 'BP'
                                    [ref_signal, ~] = custom_reref(signal, cnstrct, bad_chnls, ref, spec, true);
                                case 'leadSpecific'
                                    try
                                        [ref_signal, ~] = custom_reref(signal, cnstrct, bad_chnls, ref, RR_chnls, true);
                                    catch
                                        continue;
                                    end
                                otherwise
                                    [ref_signal, ~] = custom_reref(signal, cnstrct, bad_chnls, ref, [], true);
                            end

                            % Caluclate and extract spectral features
                            spectralFeatures = extractSpectralFeatures(ref_signal, FS, wndw_fxn, noverlap, nfft, obj.spectralParameters.bands);

                            % Convert to singles and save in spectra struct
                            spectralFeatures = struct2single(spectralFeatures);
                            spectraStruct.(ref) = spectralFeatures;
                        catch ME_ref
                            spectraStruct.(ref).ProcessingErrorMsg = ME_ref.message;
                        end
                    end
                catch ME
                    warning('Failed to calculate PSD for file %s: %s', tempData(f).name, ME.message);
                end
                % ~~~~~~~~ Save the current file structure ~~~~~~~~~~~~
                tempFilename = fullfile(outDataFolder,strrep(tempData(f).name,'.dat','.mat'));
                parsave(tempFilename, spectraStruct);
                progressbar(f, nFiles)  
            end        
        end

        function obj = detectSpikes(obj, reftypes, OVERWRITE)
            %detectSpikes - iterates over recordings and detects and
            %filteres spijkes
            % Optional Inputs:
            %   reftypes   - rereferencing methods (default: obj.rerefMethods)
            %   OVERWRITE  - whether to recompute for all files (default: false)

            % Handle optional inputs
            if nargin < 2 || isempty(reftypes)
                reftypes = obj.rerefMethods;
            end
            if nargin < 3
                OVERWRITE = false;
            end

            % Ensure reftypes is always a cell array of char
            if ischar(reftypes) || isstring(reftypes)
                reftypes = {char(reftypes)};  % convert to single-element cell
            end

            % Validate contents of reftypes
            validTypes = {'Native','CAR','BP','leadCAR','leadSpecific'};
            if any(~ismember(reftypes, validTypes))
                error('Invalid reftypes. Must be one or more of: %s', strjoin(validTypes, ', '));
            end

            % Specify spike detecion parameters
            LF =  obj.SpikeDetectionParams.LF;
            HF = obj.SpikeDetectionParams.HF;
            th_mult = obj.SpikeDetectionParams.th_mult;
            minAmp = obj.SpikeDetectionParams.Amp_th_min;
            cnstrct = obj.rerefSpecs.construct;

            % Check if spikeExamples file exists
            if isfile(obj.spikeEventsFile) && ~OVERWRITE
                % Identify which files still need spike detection
                pendingMask = obj.fileStatus.InRecordings & ~obj.fileStatus.SpikeDetect;
                pendingFiles = obj.fileStatus.FileName(pendingMask);

                % Filter recordings to process only pending ones
                recIdx = ismember({obj.recordings.name}, pendingFiles);
                recordingsToProcess = obj.recordings(recIdx);
            else
                recordingsToProcess = obj.recordings;
            end

            % Select only passive (non-stim) recordings
            nonStimMask = arrayfun(@(r) ~isfield(r, 'StimInfo') || isempty(r.StimInfo), recordingsToProcess);
            recordingsToProcess = recordingsToProcess(nonStimMask);  %  Filtered

            % Initialize structure for detected events
            spikeExamples = rmfield(recordingsToProcess, ...
                setdiff(fieldnames(recordingsToProcess), ...
                {'name', 'folder','date'}));

            % Inform the user
            nFiles = numel(recordingsToProcess);
            fprintf('Starting spike detection in %d files\n', nFiles);

            % Load recordings and detect spikes
            for i = 1:nFiles
                % Skip if already done and not overwriting
                if ~OVERWRITE && isfield(obj.recordings(i).(reftypes{1}), 'SpikeDetect')
                    continue;
                end

                try
                    % Adjust filepath if needed
                    if contains(recordingsToProcess(i).folder, obj.dataRoot)
                        fl_pth = fullfile(recordingsToProcess(i).folder,recordingsToProcess(i).name);
                    else
                        fl_pth = fullfile(strrep(recordingsToProcess(i).folder, ...
                            obj.previousDataRoot, obj.dataRoot), ...
                            recordingsToProcess(i).name);
                    end

                    % Load recording, extract recording info & parameters
                    [signal, states, parameters] = load_bcidat(fl_pth, '-calibrated');
                    signal = interpolate_lost_samples(double(signal), states); % interpolate over lost samples
                    FS = parameters.SamplingRate.NumericValue;
                    bad_chnls = recordingsToProcess(i).BadChannels.BadChannels;

                    % Re-reference the signal
                    for ridx = 1:numel(reftypes)
                        ref = reftypes{ridx};
                        try
                            switch ref
                                case 'Native'
                                    ref_signal = signal;
                                    ch_names = arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(signal,2), 'UniformOutput', false);
                                    ch_names(bad_chnls) = arrayfun(@(n) sprintf('BAD-Ch %d', n), bad_chnls, 'UniformOutput', false);
                                case 'BP'
                                    [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, obj.recordings(i).(ref).RefSpec, true);
                                case 'leadSpecific'
                                    try
                                        [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, obj.recordings(i).(ref).RefSpec, true);
                                    catch
                                        continue;
                                    end
                                otherwise
                                    [ref_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, ref, [], true);
                            end


                            % Remove Bad chanels
                            if ~strcmp(ref,'BP')        % Handle the BP rereferncing where bad channels are omited
                                ref_signal(:,bad_chnls) = [];
                                ch_names(bad_chnls) = [];
                            end


                            % -------------------- Spike detection and Filtering from HERE ------------
                            % Upsample signal to 2k Hz
                            q = 2;
                            ref_signal = resample(ref_signal,q,1,'Dimension',1);
                            montage = struct('SampleRate', q*FS, 'ChannelNames', {ch_names});

                            % Detect spikes
                            % For spike detection and elimination functions please refer to 
                                % Ayyoubi, A. H., Besheli, B. F., Swamy, C. P., Okkabaz, J. L., Miller, K. J., Worrell, G. A., & Ince, N. F. (2025). 
                                % Spurious Spike Elimination using Sparse Signal Processing Improves Seizure Onset Zone Delineation in Brief Intraoperative iEEG Recordings. 
                                % The ... Midwest Symposium on Circuits and Systems conference proceedings : MWSCAS. Midwest Symposium on Circuits and Systems, 2025, 666–670. 
                                % https://doi.org/10.1109/mwscas53549.2025.11244470
                            IIS_Det = Spike_detection_function(ref_signal, montage, [], 512, [LF, HF], minAmp, th_mult);

                            % Filter spikes
                            Spike_RFOMP = pSpike_elimination_RFOMP(IIS_Det.spike,montage.SampleRate);

                            % --- Extract spike examples (up to 10 randomly selected that passed the denoising) ---
                            denoisedSpikesIdx = find(Spike_RFOMP.feature.Pred==1);  % [samples x nSpikes]
                            nSpikes = size(denoisedSpikesIdx, 1);
                            if nSpikes > 0
                                nExamples = min(10, nSpikes);
                                randIdx = randperm(nSpikes, nExamples);
                                spikeEventSubset = IIS_Det.spike.events.Raw(randIdx,:);  % [samples x nSpikes]
                                spikeEventSubsetChLabels = ch_names(IIS_Det.spike.Channelinformation(randIdx));
                                spikeEventTimeStamp = IIS_Det.spike.timestamp(randIdx);
                            else
                                spikeEventSubset = [];  % no spikes found
                                spikeEventSubsetChLabels = [];
                                spikeEventTimeStamp = [];
                            end

                            % Save spike examples and statistic
                            recIdx = find(strcmp({obj.recordings.name}, recordingsToProcess(i).name)); % Save spike detection result back to master recording
                            IIS_Det.spike = rmfield(IIS_Det.spike,'events');        % Remove the events to save memory
                            obj.recordings(recIdx).(ref).SpikeDetect = IIS_Det.spike;
                            obj.recordings(recIdx).(ref).SpikeDetect.Montage = montage; % !!!! changed from channel names
                            obj.recordings(recIdx).(ref).SpikeDetect.RFOMP = Spike_RFOMP.feature; % Save the output of spike filtering
                            obj.recordings(recIdx).(ref).SpikeDetect.NumberOfDenoisedSpikes = sum((Spike_RFOMP.feature.Pred==1));
                            % -------------------- TO HERE --------------------

                            % Save to spikeExamples struct
                            spikeExamples(i).(ref).spikes = spikeEventSubset;
                            spikeExamples(i).(ref).channelLabels = spikeEventSubsetChLabels;
                            spikeExamples(i).(ref).timeStamps = spikeEventTimeStamp;
                        catch ME_Ref
                            obj.recordings(recIdx).ProcessingErrorMsg = ME_Ref;
                        end
                    end
                catch
                    warning('Spike detection failed on file %s: %s', obj.recordings(i).name);
                    continue;
                end

                % --- Update fileStatus ---
                fileIdx = find(strcmp(obj.fileStatus.FileName, recordingsToProcess(i).name));
                if any(fileIdx)
                    obj.fileStatus.SpikeDetect(fileIdx) = true;
                end
                progressbar(i, nFiles)
            end

            % Save or append spikeExamples
            if isfile(obj.spikeEventsFile) && ~OVERWRITE
                % Append new entries to existing
                loaded = load(obj.spikeEventsFile, 'spikeExamples');
                spikeExamples = [loaded.spikeExamples; spikeExamples];
                save(obj.spikeEventsFile, 'spikeExamples', '-v7.3');
            else
                save(obj.spikeEventsFile, 'spikeExamples', '-v7.3');
            end

            disp(['Spike event examples saved to: ', obj.spikeEventsFile]);
        end


        %% ~~~~~~~~~~~~~~~~ Load & Save the Data and Object ~~~~~~~~~~~~~~~
        function saveData(obj, whatToSave)
            % saveData - Save selected internal data structures to .mat files
            % Usage:
            %   obj.saveData();                  % Save all
            %   obj.saveData('recordings');      % Save only recordings
            %   obj.saveData({'impedance','object'}) % Save multiple targets

            arguments
                obj
                whatToSave (1,:) string {mustBeMember(whatToSave, ...
                    ["all", "recordings", "impedance", "longterm", "discarded", "object"])} = "all"
            end

            root = obj.dataRoot;

            % ----- Recordings -----
            if any(strcmp(whatToSave, 'all')) || any(strcmp(whatToSave, 'recordings'))
                fpath = fullfile(root, 'recordings.mat');
                if ~isempty(obj.recordings) && confirmOverwrite(fpath)
                    dataFiles = obj.recordings;
                    save(fpath, 'dataFiles', '-v7.3');
                    disp('Saved recordings.mat');
                end
            end

            % ----- Impedance -----
            if any(strcmp(whatToSave, 'all')) || any(strcmp(whatToSave, 'impedance'))
                fpath = fullfile(root, 'impedance.mat');
                if ~isempty(obj.impedanceData) && confirmOverwrite(fpath)
                    impedance = obj.impedanceData;
                    save(fpath, 'impedance', '-v7.3');
                    disp('Saved impedance.mat');
                end
            end

            % ----- Long-term stim -----
            if any(strcmp(whatToSave, 'all')) || any(strcmp(whatToSave, 'longterm'))
                fpath = fullfile(root, 'longterm_stim.mat');
                if ~isempty(obj.longtermStimData) && confirmOverwrite(fpath)
                    dataFiles = obj.longtermStimData;
                    save(fpath, 'dataFiles', '-v7.3');
                    disp('Saved longterm_stim.mat');
                end
            end

            % ----- Discarded -----
            if any(strcmp(whatToSave, 'all')) || any(strcmp(whatToSave, 'discarded'))
                fpath = fullfile(root, 'discarded.mat');
                if ~isempty(obj.discardedRecordings) && confirmOverwrite(fpath)
                    discarded = obj.discardedRecordings;
                    save(fpath, 'discarded', '-v7.3');
                    disp('Saved discarded.mat');
                end
            end

            % ----- Object (lightweight) -----
            if any(strcmp(whatToSave, 'all')) || any(strcmp(whatToSave, 'object'))
                % Determine filename
                if ~isempty(obj.subjectName)
                    objectFile = sprintf('%s.mat', obj.subjectName);
                else
                    [~, folderName] = fileparts(obj.dataRoot);
                    objectFile = sprintf('%s.mat', folderName);
                end
                fpath = fullfile(root, objectFile);

                if confirmOverwrite(fpath)
                    objCopy = obj;
                    objCopy.recordings = [];
                    objCopy.impedanceData = [];
                    objCopy.longtermStimData = [];
                    objCopy.discardedRecordings = [];

                    save(fpath, 'objCopy', '-v7.3');
                    disp(['Saved object as ', objectFile]);
                end
            end
        end

        function obj = loadData(obj, dirPath, whichVars)
            % loadData - Load .mat data files into the object from dataRoot or specified directory
            % Inputs:
            %   dirPath    - (optional) directory to load from (default = obj.dataRoot)
            %   whichVars  - (optional) string or string array of fields to load
            % Valid values for whichVars: 'recordings', 'impedance', 'longterm', 'discarded'

            arguments
                obj
                dirPath {mustBeFolder}
                whichVars (1, :) string {mustBeMember(whichVars, ["recordings", "impedance", "longterm", "discarded"])} = ...
                    ["recordings", "impedance", "longterm", "discarded"]
            end

            if ~isfolder(dirPath)
                error('Specified directory "%s" does not exist.', dirPath);
            end

            % --- Update dataRoot dynamically ---
            if ~strcmp(dirPath , obj.dataRoot)
                obj.previousDataRoot = obj.dataRoot;
                obj.dataRoot = dirPath;
                disp(['Updated obj.dataRoot to: ', obj.dataRoot]);
            end

            % Define mapping of logical names to file names and object properties
            varMap = struct( ...
                'recordings', struct('file', 'recordings.mat',      'var', 'dataFiles', 'prop', 'recordings'), ...
                'impedance',  struct('file', 'impedance.mat',       'var', 'impedance',  'prop', 'impedanceData'), ...
                'longterm',   struct('file', 'longterm_stim.mat',   'var', 'dataFiles',  'prop', 'longtermStimData'), ...
                'discarded',  struct('file', 'discarded.mat',       'var', 'discarded',  'prop', 'discardedRecordings') ...
                );

            for v = whichVars
                info = varMap.(v);
                filepath = fullfile(dirPath, info.file);

                if isfile(filepath)
                    loaded = load(filepath);
                    if isfield(loaded, info.var)
                        obj.(info.prop) = loaded.(info.var);
                        disp(['Loaded ', info.prop, ' from ', filepath]);
                    else
                        warning('Expected variable "%s" not found in "%s". Skipping.', info.var, info.file);
                    end
                else
                    warning('File "%s" not found in directory "%s". Skipping.', info.file, dirPath);
                end
            end
        end

        function [signal, recInfo] = loadRecording(obj, file)
            % Loads individual recording by specifying it's index
            % in recordings structure or by it's name with .dat extension
            % Inputs: -> file :index or filename with name with .dat extension
            % Outputs: -> signal : raw signal extracted from .dat file (if .dat file exists)
            %          -> recInfo : struct containing information and PSDs (if .mat file exists)

            % Handle inputs
            if nargin < 2 || isempty(file)
                error('No file was specified. Please specify file by it''s name or index');
            else % Handle file (index or name)
                if isnumeric(file)
                    idx = file;
                elseif ischar(file) || isstring(file)
                    idx = find(strcmp({obj.recordings.name}, char(file)), 1);
                    if isempty(idx)
                        error('Specified file "%s" not found in obj.recordings.', file);
                    end
                else
                    error('Invalid input: fileSpec must be index or filename.');
                end
                fname = obj.recordings(idx).name;
            end

            % Extract filenames
            datFile = fullfile(obj.recordings(idx).folder, fname);
            matFile = fullfile(obj.recordingsMATFolder, strrep(fname,'.dat','.mat'));
            recInfo = obj.recordings(idx);

            if isfile(datFile)
                signal = load_bcidat(datFile, '-calibrated');
            else
                warning('Recording file %s does not longer exist on specified path.', fname);
                signal = [];
            end

            if isfile(matFile)
                spectra = load(matFile);
                reftypes = intersect(fieldnames(recInfo),fieldnames(spectra));
                for rf=1:numel(reftypes)
                    flds = fieldnames(spectra.(reftypes{rf}));
                    for fn = 1:numel(flds)
                        recInfo.(reftypes{rf}).(flds{fn}) = spectra.(reftypes{rf}).(flds{fn});
                    end
                end
            else
                warning('Recording file %s does not longer exist in %s folder.', fname, obj.recordingsMATFolder);
            end
        end

        %% ~~~~~~~~~~~~~~~ Stats & Summaries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function obj = getSpectralStats(obj, rerefTypes, OVERWRITE)
            % getSpectralStats - Extracts spectral statistics (RMS power) from PSDs
            % Inputs:
            %   rerefTypes - rereferencing methods to analyze (default = obj.rerefMethods)
            %   OVERWRITE  - whether to recompute stats even if already present

            if nargin < 2 || isempty(rerefTypes)
                rerefTypes = obj.rerefMethods;
            end
            if nargin < 3
                OVERWRITE = false;
            end

            bands = obj.spectralParameters.bands;
            bandNames = fieldnames(bands);
            % Inform the user
            disp('Extracting statistics from spectral data');
            nFiles = numel(obj.recordings);

            for i = 1:nFiles
                for r = 1:length(rerefTypes)
                    ref = rerefTypes{r};

                    % Skip if stats already exist and overwrite is false
                    if isfield(obj.recordings(i), ref) && ...
                            isfield(obj.recordings(i).(ref), 'SpectralStats') && ...
                            ~OVERWRITE
                        continue;
                    end

                    % Try to load PSD from recordings_MAT folder
                    matFilename = fullfile(obj.recordingsMATFolder, strrep(obj.recordings(i).name, '.dat', '.mat'));
                    if ~isfile(matFilename)
                        % warning('Missing PSD .mat file for %s. Skipping.', obj.recordings(i).name);
                        continue;
                    end
                    try
                        data = load(matFilename);
                    catch
                        continue
                    end
                    if ~isfield(data, ref) || ~isfield(data.(ref), 'PSD')
                        % warning('PSD not found in %s for ref %s. Skipping.', matFilename, ref);
                        continue;
                    end

                    psd = data.(ref).PSD;
                    freqs = data.(ref).info.Frequencies;
                    binWidth = mean(diff(freqs));

                    % Initialize stats
                    stats = struct();
                    stats.TotalRMS = sqrt(sum(psd * binWidth));

                    for b = 1:numel(bandNames)
                        name = bandNames{b};
                        range = bands.(name);
                        idx = freqs >= range(1) & freqs <= range(2);
                        stats.([name 'RMS']) = sqrt(sum(psd(idx, :) * binWidth));
                    end

                    %  Avg, Med, Std PSD (per frequency across channels)
                    stats.AvgPSD = mean(psd, 2);     
                    stats.MedPSD = median(psd, 2);   
                    stats.StdPSD = std(psd, 0, 2); 

                    % Coherence 
                    stats.Frequencies = freqs;
                    % Store in-memory
                    obj.recordings(i).(ref).SpectralStats = stats;
                end
                progressbar(i,nFiles)
            end

        end

        function [obj,summaryTable] = summarizeSubject(obj)
            % summarizeSubject - Computes and stores subject summary.
            % If output is requested, also returns a table version.

            recs = obj.recordings;
            if isempty(recs)
                warning('No recordings available.');
                if nargout > 1, summaryTable = table(); end
                return;
            end
            % Total duration
            totalDuration = sum(cellfun(@duration,{recs.Duration}));

            % Filter non-stim recordings
            nonStim = recs(~arrayfun(@(r) isfield(r, 'StimInfo') && ~isempty(r.StimInfo), recs));

            % Total duration of non-stim only (for normalization)
            durationNonStim = sum(cellfun(@duration,{nonStim.Duration}));
            durationNonStim = minutes(durationNonStim);

            % Time since implant
            if isempty(obj.implantDate)
                daysSinceImplant = NaN;
            else
                daysSinceImplant = fix(days(datetime('now') - obj.implantDate));
            end

            % Unique recording days
            recDates = unique(datetime({recs.date}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'));
            nonStimRecDates = dateshift(datetime({nonStim.date}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'), 'start', 'day'); % For better bad channel by month estiamtion
            uniqueRecDays = numel(unique(dateshift(recDates, 'start', 'day')));

            % Unique impedance days
            impDates = datetime({obj.impedanceData.date}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
            uniqueImpDays = numel(unique(dateshift(impDates, 'start', 'day')));

            % Packet loss
            pktLoss = arrayfun(@(r) double(r.PacketLoss), recs);
            avgPacketLoss = mean(pktLoss, 'omitnan');

            % Bad channels
            badChCount = arrayfun(@(r) numel(r.BadChannels.BadChannels), nonStim);
            avgBadChCount = mean(badChCount, 'omitnan');
           
            % Bad channel occurrence 
            allBadChs = arrayfun(@(r) r.BadChannels.BadChannels, nonStim, 'UniformOutput', false);
            allBadChs_array = cat(2,allBadChs{:});
            chNames = arrayfun(@(ch) sprintf('Ch %d', ch), allBadChs_array, 'UniformOutput', false);
            badChTbl = tabulate(chNames);  
            OverallBadChOccurrence = table(badChTbl(:,1), cell2mat(badChTbl(:,2)),badChTbl(:,3), ...
                'VariableNames', {'Channel', 'Count','Percentage'});  

            % Convert to month-level grouping (only non stim recordings)
            monthDates = dateshift(nonStimRecDates, 'end', 'month');
            uniqueMonths = unique(monthDates);
            
            % Initialize result struct
            monthlyBadChStats = struct;
            for m = 1:numel(uniqueMonths)
                thisMonth = uniqueMonths(m);
                % Create valid struct field name for month
                monthField = datestr(thisMonth, 'mmm_yyyy');  % e.g. 'Sep_2023'
                monthField = matlab.lang.makeValidName(monthField);
                idx = monthDates == thisMonth;
                % Collect all bad channels for this month
                badChThisMonth = allBadChs(idx);
                % Concatenate all channel indices into one vector
                allBadChThisMonth = [badChThisMonth{:}];
                if isempty(allBadChThisMonth)
                    % Store in dictionary
                    monthlyBadChStats.(monthField).NumberOfRecordings = sum(idx);
                    monthlyBadChStats.(monthField).BadChannelOccurence = [];
                    continue
                end
                % Create table
                badChTbl = tabulate(allBadChThisMonth); 
                nonemptyIdx = badChTbl(:,2) ~= 0; % drop empty fields
                badChTbl = table(badChTbl(nonemptyIdx,1), badChTbl(nonemptyIdx,2), (badChTbl(nonemptyIdx,2)./sum(idx))*100,...
                 'VariableNames', {'Channel', 'Count', 'Percentage'});
                % Store table
                monthlyBadChStats.(monthField).NumberOfRecordings = sum(idx);
                monthlyBadChStats.(monthField).BadChannelOccurence = badChTbl;
            end

            % Average RMS across rereferencing types
            rerefTypes = obj.rerefMethods;
            avgRMS = struct();
            for r = 1:numel(rerefTypes)
                ref = rerefTypes{r};
                totalRMS = [];
                for i = 1:numel(nonStim)
                    if isfield(nonStim(i), ref) && isfield(nonStim(i).(ref), 'SpectralStats')
                        rms = nonStim(i).(ref).SpectralStats.TotalRMS;
                        if ~isempty(rms)
                            totalRMS(end+1) = mean(rms);  % average across channels
                        end
                    end
                end
                avgRMS.(ref) = mean(totalRMS, 'omitnan');
            end

            % Spike counts and normalized counts
            spikeCounts = struct();
            spikeCountsNorm = struct();
            totalSpikePerChannel = struct();

            for r = 1:numel(rerefTypes)
                ref = rerefTypes{r};
                totalSpikes = 0;
                norm_spikeCount = nan(numel(nonStim),1);
                spikeMap = containers.Map('KeyType', 'char', 'ValueType', 'double');  % aggregate spike counts

                for i = 1:numel(nonStim)
                    if isfield(nonStim(i), ref) && isfield(nonStim(i).(ref), 'SpikeDetect') && ...
                            isfield(nonStim(i).(ref).SpikeDetect, 'RFOMP') && ...
                            isfield(nonStim(i).(ref).SpikeDetect.RFOMP, 'Pred')

                        badChnls = nonStim(i).BadChannels.BadChannels;
                        chanNames = nonStim(i).(ref).SpikeDetect.Montage.ChannelNames;      % Retrieve the channel names
                        spikeChannels = nonStim(i).(ref).SpikeDetect.Channelinformation(nonStim(i).(ref).SpikeDetect.RFOMP.Pred==1); % Retrieve channel with spikes
                        % Remove spikes detected on bad channels
                        spikeChannels(ismember(spikeChannels,badChnls)) = [];
                        totalSpikes = numel(spikeChannels);

                        % Accumulate into map
                        for c = 1:numel(spikeChannels)
                            ch = chanNames{spikeChannels(c)};

                            if spikeMap.isKey(ch)
                                spikeMap(ch) = spikeMap(ch) + 1;
                            else
                                spikeMap(ch) = 1;
                            end
                        end
                        rec_dur = minutes(duration(nonStim(i).Duration, 'InputFormat', 'hh:mm:ss.SSS'));
                        norm_spikeCount(i) = totalSpikes/(length(chanNames)*rec_dur); % Normalized spike count by dividing through number of channels and duration of recording
                    end
                end

                % Convert map to table
                allKeys = keys(spikeMap);
                allValues = values(spikeMap);
                spikeTbl = table(allKeys', cell2mat(allValues'), ...
                    'VariableNames', {'Channel', 'SpikeCount'});

                % Save
                spikeCounts.(ref) = totalSpikes;
                spikeCountsNorm.(ref) = mean(norm_spikeCount,'omitnan');
                spikeCountsNorm.([ref '_STD']) = std(norm_spikeCount,'omitnan');
                totalSpikePerChannel.(ref) = spikeTbl;
            end

            % Store summary as struct
            summary = struct();
            summary.TotalDuration = totalDuration;
            summary.NonStimRecordingsDuration = durationNonStim;
            summary.DaysSinceImplant = daysSinceImplant;
            summary.UniqueRecordingDays = uniqueRecDays;
            summary.UniqueImpedanceDays = uniqueImpDays;
            summary.AvgPacketLoss = avgPacketLoss;
            summary.AvgBadChCount = avgBadChCount;
            summary.BadChannelOccurrence = OverallBadChOccurrence;
            summary.MonthlyBadChannelOccurrence = monthlyBadChStats;
            summary.AvgTotalRMS = avgRMS;
            summary.SpikeCounts = spikeCounts;
            summary.SpikeCountsNormalized = spikeCountsNorm;
            summary.SpikeCountPerChannel = totalSpikePerChannel;

            % Save to property
            obj.subjectSummary = summary;

            % Return output as table if requested
            if nargout > 1
                summaryTable = table( ...
                    summary.TotalDuration, ...
                    summary.DaysSinceImplant, ...
                    summary.UniqueRecordingDays, ...
                    summary.UniqueImpedanceDays, ...
                    summary.AvgPacketLoss, ...
                    summary.AvgBadChCount, ...
                    'VariableNames', { ...
                    'TotalDurationHrs', ...
                    'DaysSinceImplant', ...
                    'RecordingDays', ...
                    'ImpedanceDays', ...
                    'AvgPacketLoss', ...
                    'AvgBadChCount' ...
                    });
            end
        end

        function impTable = summarizeImpedance(obj)
            % summarizeImpedance - Creates a summary table of average impedances per day
            %
            % Returns:
            %   impTable - table with columns: Date, DaysSinceImplant, Ch1...Ch32

            if isempty(obj.impedanceData)
                warning('No impedance data found.');
                impTable = table();
                return;
            end

            if isempty(obj.implantDate)
                warning('implantDate is not defined.');
                impTable = table();
                return;
            end

            % Convert dates to datetime and align to day
            impDates = dateshift( ...
                datetime({obj.impedanceData.date}, 'InputFormat', 'd-MMM-y HH:mm:ss'), ...
                'start', 'day');

            uniqueDates = unique(impDates, 'sorted')';
            if ~isempty(obj.excludedDates)
                [~,excludeMask] = intersect(uniqueDates,obj.excludedDates);
                uniqueDates(excludeMask) = [];
            end
            replace_inf = 4e4;

            plotData = struct();
            for i = 1:length(uniqueDates)
                day = uniqueDates(i);
                dayMask = impDates == day;
                imp = obj.impedanceData(dayMask);  % All files for this day

                % Metadata
                plotData(i).Date = day;  
                plotData(i).DaysSinceImplant = days(day - obj.implantDate);

                % Channel loop
                for ch = 1:32
                    chName = sprintf('Ch%d', ch);
                    z_ch = {imp.(chName)};  % array of values from all files
                    z_ch = cell2mat(z_ch);
                    z_ch(isinf(z_ch)) = replace_inf;
                    plotData(i).(sprintf('Impedance_Ch%d', ch)) = mean(z_ch, 'omitnan');
                end
            end

            % Convert to table
            impTable = struct2table(plotData);
        end

        function [summaryTable, glmTable] = summarizeLongtermStim(obj, refType)
            % summarizeLongtermStim Summarize long-term stimulation features.
            % Outputs
            %   summaryTable - group-level summary table
            %   glmTable     - recording-level long-format table for GLM/GLMM
            % Example
            %   [summaryTable, glmTable] = obj.summarizeLongtermStim('Native');

            arguments
                obj
                refType (1,:) char {mustBeMember(refType, ...
                    {'Native','CAR','leadCAR','leadSpecific','BP'})} = 'Native'
            end

            % Ensure collected features exist
            if isempty(obj.longtermStimFeatures)
                [obj, stimFeatures] = obj.collectLongtermStimFeatures();
            else
                stimFeatures = obj.longtermStimFeatures;
            end

            groupNames = stimFeatures.groups;
            nGroups = numel(groupNames);

            % Group-level summary variables
            Condition = strings(nGroups,1);
            NumRecordings = zeros(nGroups,1);
            TotalDuration = duration.empty(nGroups,0);
            NumRecordingDays = zeros(nGroups,1);
            MeanDurationPerDay = duration.empty(nGroups,0);

            % Recording-level GLM table
            glmTable = table();
            implantDate = obj.implantDate;

            for g = 1:nGroups
                conditionName = groupNames{g};
                Condition(g) = string(conditionName);
                recStruct = stimFeatures.(conditionName);
                NumRecordings(g) = numel(recStruct);

                TotalDuration(g,1) = stimFeatures.totalDuration.(conditionName);
                dayTable = stimFeatures.RecDays.(conditionName);
                NumRecordingDays(g) = height(dayTable);

                if height(dayTable) > 0
                    MeanDurationPerDay(g,1) = mean(dayTable.Duration);
                else
                    MeanDurationPerDay(g,1) = duration(0,0,0);
                end

                % Stimulation frequency
                if strcmpi(conditionName, 'baseline')
                    stimFrequency_Hz = 0;
                    conditionLabel = "Baseline";
                else
                    tok = regexp(conditionName, '(\d+)Hz', 'tokens', 'once');
                    stimFrequency_Hz = str2double(tok{1});
                    conditionLabel = "Stim" + string(stimFrequency_Hz) + "Hz";
                end

                % Recording-level rows
                for f = 1:numel(recStruct)
                    rec = recStruct(f);
                    if ~isempty(rec.(refType))
                        F = rec.(refType);
                    else
                        continue
                    end

                    recID = erase(string(rec.name), ".dat");
                    recordingDate = rec.date;
                    recordingDuration_min = minutes(rec.duration);

                    % Day ID
                    dayID = find(dateshift(dayTable.Date, 'start', 'day') == ...
                        dateshift(recordingDate, 'start', 'day'), 1);
                    % Recording order within condition and day
                    sameDayIdx = arrayfun(@(x) dateshift(x.date, 'start', 'day') == ...
                        dateshift(recordingDate, 'start', 'day'), recStruct);
                    sameDayRecs = recStruct(sameDayIdx);
                    sameDayDates = [sameDayRecs.date];
                    [~, orderIdx] = sort(sameDayDates);
                    originalIdx = find(sameDayIdx);
                    recordingOrder = find(originalIdx(orderIdx) == f, 1);
                    daysSinceImplant = days(recordingDate - implantDate); % Days since implant

                    %---- Spike features -----
                    spikeCount = F.nSpikes;
                    spikeRate = F.normSpikeCount;

                    % --- Exposure ---
                    nValidChannels = numel(F.GoodChannels);
                    exposure_channel_min = recordingDuration_min .* nValidChannels;
                    logExposure = log(exposure_channel_min);

                    % ---- Base GLM row ----
                    newRow = table( ...
                        string(conditionName), ...
                        conditionLabel, ...
                        stimFrequency_Hz, ...
                        recID, ...
                        string(rec.name), ...
                        rec.sourceRecordingIdx, ...
                        recordingDate, ...
                        dayID, ...
                        recordingOrder, ...
                        daysSinceImplant, ...
                        recordingDuration_min, ...
                        nValidChannels, ...
                        exposure_channel_min, ...
                        logExposure, ...
                        spikeCount, ...
                        spikeRate, ...
                        string(refType), ...
                        'VariableNames', { ...
                        'stimGroup', ...
                        'condition', ...
                        'stimFrequency_Hz', ...
                        'recID', ...
                        'fileName', ...
                        'sourceRecordingIdx', ...
                        'recordingDate', ...
                        'dayID', ...
                        'recordingOrder', ...
                        'daysSinceImplant', ...
                        'recordingDuration_min', ...
                        'nValidChannels', ...
                        'exposure_channel_min', ...
                        'logExposure', ...
                        'spikeCount', ...
                        'spikeRate', ...
                        'refType'});

                    % ---- Spectral features ----
                    skipFields = {'nSpikes','spikeCount','numSpikes', ...
                        'normSpikeCount','spikeRate','normalizedSpikeRate', ...
                        'GoodCHannels','GoodChannels','goodChannels','validChannels'};

                    featureNames = fieldnames(F);
                    spectralRow = table();
                    for k = 1:numel(featureNames)
                        featureName = featureNames{k};
                        if any(strcmp(featureName, skipFields))
                            continue
                        end
                        value = F.(featureName);
                        varBase = matlab.lang.makeValidName(featureName);

                        if isnumeric(value) && isscalar(value)
                            spectralRow.(varBase) = value;
                        elseif isnumeric(value) && isvector(value) && ~isempty(value)
                            spectralRow.([varBase '_mean']) = mean(value, 'omitnan');
                            spectralRow.([varBase '_median']) = median(value, 'omitnan');
                            spectralRow.([varBase '_std']) = std(value, 'omitnan');
                        elseif isnumeric(value) && ismatrix(value) && ~isempty(value)
                            spectralRow.([varBase '_mean']) = mean(value(:), 'omitnan');
                            spectralRow.([varBase '_median']) = median(value(:), 'omitnan');
                            spectralRow.([varBase '_std']) = std(value(:), 'omitnan');
                        end
                    end
                    newRow = [newRow spectralRow];
                    glmTable = [glmTable; newRow];
                end
            end

            % Final group-level summary table
            summaryTable = table(Condition, NumRecordings, TotalDuration, ...
                NumRecordingDays, MeanDurationPerDay, ...
                'VariableNames', {'StimType', 'NumberOfRecordings', ...
                'TotalDuration', 'NumberOfRecDays', 'MeanDurationPerDay'});

            % Final formatting
            glmTable.condition = categorical(glmTable.condition);
            glmTable.stimGroup = categorical(glmTable.stimGroup);
            glmTable.refType = categorical(glmTable.refType);

            % Store summary
            stimSummary = struct();
            stimSummary.summaryTable = summaryTable;
            stimSummary.glmTable = glmTable;
            stimSummary.groups = groupNames;
            stimSummary.refType = refType;

            for g = 1:nGroups
                conditionName = groupNames{g};
                stimSummary.(conditionName).NumberOfRecordings = NumRecordings(g);
                stimSummary.(conditionName).TotalDuration = TotalDuration(g);
                stimSummary.(conditionName).NumberOfRecDays = NumRecordingDays(g);
                stimSummary.(conditionName).MeanDurationPerDay = MeanDurationPerDay(g);
            end
            obj.longtermStimSummary = stimSummary;
        end


        %% ~~~~~~~~~~~~ Feature & Data extraction  ~~~~~~~~~~~~~~~~~~~~~~~~
        function psdStructs = collectPSDs(obj, file, rerefTypes)
            % collectPSDs - Collects PSDs and metadata from recordings and .mat files
            %
            % Inputs (optional):
            %   file    - index into obj.recordings OR filename ('.dat')
            %   rerefTypes  - rereferencing types to load (default = obj.rerefMethods)
            %
            % Output:
            %   psdStructs  - struct array with fields:
            %                 'name', 'date', 'Ground', 'Stim', 'Frequencies', reref PSDs
            % --- Handle defaults ---
            PSDmatFile = fullfile(obj.recordingsMATFolder, 'PSDs.mat');
            if  nargin < 3 || isempty(rerefTypes)
                rerefTypes = obj.rerefMethods;
            end

            % --- If no fileSpec is provided ---
            if nargin < 2 || isempty(file)
                % Check if 'PSDs.mat' exists and offer to load
                if isfile(PSDmatFile)
                    choice = questdlg('PSDs.mat already exists. Load it or recompute?', ...
                        'Load or Recompute', 'Load', 'Recompute', 'Load');

                    if strcmp(choice, 'Load')
                        psdStructs = load(PSDmatFile, 'psdStructs');
                        psdStructs = psdStructs.psdStructs;
                        disp('Loaded PSDs.mat from disk.');
                        return;
                    end
                else
                    choice = 'Recompute'; % If no PSDmatFile exists, compute from scratch
                end

                % If recompute: process all recordings
                recs = obj.recordings;
                recIndices = 1:numel(recs);
            else
                % --- Handle fileSpec (index or name) ---
                if isnumeric(file)
                    recIndices = file;
                elseif ischar(file) || isstring(file)
                    idx = find(strcmp({obj.recordings.name}, char(file)), 1);
                    if isempty(idx)
                        error('Specified file "%s" not found in obj.recordings.', file);
                    end
                    recIndices = idx;
                else
                    error('Invalid input: fileSpec must be index or filename.');
                end
                recs = obj.recordings(recIndices);
            end

            n = numel(recIndices);
            % --- Build field template dynamically based on rerefTypes ---
            baseFields = struct( ...
                'name','', ...
                'date','', ...
                'Ground',[], ...
                'Stim',[], ...
                'Frequencies',[]);

            for r = 1:numel(rerefTypes)
                baseFields.(rerefTypes{r}) = [];
            end
            psdStructs = repmat(baseFields, n, 1);

            disp('Collecting PSDs and metadata from recordings.');
            for ii = 1:n
                i = recIndices(ii);
                try
                    psdEntry = struct();
                    psdEntry.name = recs(ii).name;
                    psdEntry.date = recs(ii).date;
                    psdEntry.Ground = recs(ii).Ground;

                    % --- Determine stim status ---
                    psdEntry.Stim = isfield(recs(ii), 'StimInfo') && ~isempty(recs(ii).StimInfo);

                    % --- Load .mat file ---
                    matFilename = fullfile(obj.recordingsMATFolder, strrep(recs(ii).name, '.dat', '.mat'));
                    if ~isfile(matFilename)
                        warning('Missing .mat file for %s. Skipping.', recs(ii).name);
                        psdStructs(ii) = psdEntry;
                        continue;
                    end

                    try
                        data = load(matFilename);
                    catch
                        warning('Unable to load .mat file for %s. Skipping PSD loading.', recs(i).name);
                        continue
                    end

                    for r = 1:numel(rerefTypes)
                        ref = rerefTypes{r};
                        if r == 1 && isfield(data, ref) && isfield(data.(ref), 'info')
                            psdEntry.Frequencies = data.(ref).info.Frequencies;
                        end
                        if isfield(data, ref) && isfield(data.(ref), 'PSD')
                            psdEntry.(ref) = data.(ref).PSD;
                        else
                            psdEntry.(ref) = [];
                        end
                    end

                    psdStructs(ii) = psdEntry;

                catch ME
                    warning('Error reading file %s: %s', recs(ii).name, ME.message);
                    continue;  % fallback to empty if failed
                end
                progressbar(i,n)
            end
            % --- Confirm before overwriting ---
            if nargin<2 && strcmp(choice,'Recompute') && confirmOverwrite(PSDmatFile) % full set
                save(PSDmatFile, 'psdStructs', '-v7.3');
            end
        end

        function psdStructs = collectLongtermStimPSDs(obj, fileSpec, rerefTypes)
            % collectLongtermStimPSDs - Collects rest/stim PSDs from long-term stim .mat files
            % Inputs (optional):
            %   fileSpec    - index or name of a specific .dat file
            %   rerefTypes  - rereferencing types to include (default = obj.rerefMethods)
            %
            % Output:
            %   psdStructs  - struct array with fields:
            %                'name', 'date', reref-specific RestPSDs and StimPSDs

            fpath = fullfile(obj.longtermStimFolder, 'LongtermStim_PSDs.mat');

            if  nargin < 3 || isempty(rerefTypes)
                rerefTypes = obj.rerefMethods;
            end

            % --- If no fileSpec is provided, offer to load existing file ---
            if nargin < 2 || isempty(fileSpec)
                if isfile(fpath)
                    choice = questdlg('LongtermStim_PSDs.mat exists. Load it or recompute?', ...
                        'Load or Recompute', 'Load', 'Recompute', 'Load');
                    if strcmp(choice, 'Load')
                        loaded = load(fpath);
                        psdStructs = loaded.psdStructs;
                        disp('Loaded LongtermStim_PSDs.mat');
                        return;
                    end
                end
                recs = obj.longtermStimData;
                recIndices = 1:numel(recs);
            else
                % --- Handle fileSpec (index or name) ---
                if isnumeric(fileSpec)
                    recIndices = fileSpec;
                elseif ischar(fileSpec) || isstring(fileSpec)
                    idx = find(strcmp({obj.longtermStimData.name}, char(fileSpec)), 1);
                    if isempty(idx)
                        error('Specified file "%s" not found in obj.longtermStimData.', fileSpec);
                    end
                    recIndices = idx;
                else
                    error('Invalid input: fileSpec must be index or filename.');
                end
                recs = obj.longtermStimData(recIndices);
            end

            n = numel(recIndices);

            % --- Build dynamic struct template ---
            baseFields =  struct('name','', 'date','');
            for r = 1:numel(rerefTypes)
                baseFields.([rerefTypes{r} '_Rest']) = [];
                baseFields.([rerefTypes{r} '_Stim']) = [];
            end
            psdStructs = repmat(baseFields, n, 1);

            for ii = 1:n
                i = recIndices(ii);
                try
                    psdEntry.name = recs(ii).name;
                    psdEntry.date = recs(ii).date;

                    matFile = fullfile(obj.longtermStimFolder, strrep(recs(ii).name, '.dat', '.mat'));
                    if ~isfile(matFile)
                        warning('Missing .mat file for %s. Skipping.', recs(ii).name);
                        psdStructs(ii) = psdEntry;
                        continue;
                    end

                    data = load(matFile);

                    for r = 1:numel(rerefTypes)
                        ref = rerefTypes{r};

                        if isfield(data, ref)
                            % --- Rest PSDs ---
                            if isfield(data.(ref), 'RestPSDs')
                                restCells = data.(ref).RestPSDs;
                                if iscell(restCells) && numel(restCells) > 1
                                    psdEntry.([ref '_Rest']) = mean(cat(3, restCells{:}), 3);
                                elseif iscell(restCells) && numel(restCells) == 1
                                    psdEntry.([ref '_Rest']) = restCells{1};
                                end
                            else
                                psdEntry.([ref '_Rest']) = [];
                            end

                            % --- Stim PSDs ---
                            if isfield(data.(ref), 'StimPSDs')
                                stimCells = data.(ref).StimPSDs;
                                if iscell(stimCells) && numel(stimCells) > 1
                                    psdEntry.([ref '_Stim']) = mean(cat(3, stimCells{:}), 3);
                                elseif iscell(stimCells) && numel(stimCells) == 1
                                    psdEntry.([ref '_Stim']) = stimCells{1};
                                end
                            else
                                psdEntry.([ref '_Stim']) = [];
                            end
                        else
                            psdEntry.([ref '_Rest']) = [];
                            psdEntry.([ref '_Stim']) = [];
                        end
                    end
                    psdStructs(ii) = psdEntry;
                catch ME
                    warning('Error processing file %s: %s', recs(ii).name, ME.message);
                    psdStructs(ii) = [];
                end
                progressbar(ii, n);
            end

            % --- Ask before overwriting ---
            if numel(recIndices) == numel(obj.longtermStimData)
                if isfile(fpath)
                    choice = questdlg('Overwrite existing LongtermStim_PSDs.mat?', ...
                        'Confirm Overwrite', 'Yes', 'No', 'No');
                    if ~strcmp(choice, 'Yes')
                        disp('User cancelled saving.');
                        return;
                    end
                end
                save(fpath, 'psdStructs', '-v7.3');
                disp(['Saved to ', fpath]);
            end
        end

        function [obj, stimFeatures] = collectLongtermStimFeatures(obj, fileSpec, rerefMethod, channelSpec)
            % collectLongtermStimFeatures Collect and summarize long-term stimulation features.
            %
            % Syntax:
            %   [obj, stimFeatures] = collectLongtermStimFeatures(obj, fileSpec)
            %
            % Aggregates stimulation/rest PSDs, band-power features, spike counts,
            % recording durations, and recording-day statistics from long-term
            % stimulation recordings. Results are grouped by stimulation protocol,
            % stored in obj.longtermStimFeatures, and optionally saved to
            % LongtermStim_features.mat.
            %
            % fileSpec (optional): [] = all recordings, index, or filename.

            % Ensure that the channel specification is string
            if nargin < 4 || isempty(channelSpec)
                channelSpec = "all";
            elseif ischar(channelSpec) || isstring(channelSpec)
                channelSpec = string(channelSpec);
            else
                error('Invalid input: channelSpec must be a channel name as char/string or "all".');
            end

            if strcmpi(channelSpec, "all")
                fpath = fullfile(obj.longtermStimFolder, 'LongtermStim_features.mat');
            else
                safeChannelSpec = matlab.lang.makeValidName(channelSpec);
                fpath = fullfile(obj.longtermStimFolder, sprintf('LongtermStim_features_%s.mat', safeChannelSpec));
            end

            % If no reref method was specified go through all of them
            if  nargin < 3 || isempty(rerefMethod)
                rerefMethod = obj.rerefMethods;
            else
                rerefMethod = {rerefMethod};
            end

            % Load existing or determine recordings to process
            if nargin < 2 || isempty(fileSpec)
                if isfile(fpath)
                    choice = questdlg( ...
                        'LongtermStim_features.mat exists. Load it or recompute?', ...
                        'Load or Recompute', ...
                        'Load', 'Recompute', 'Load');
                    if strcmp(choice, 'Load')
                        loaded = load(fpath, 'stimFeatures');
                        stimFeatures = loaded.stimFeatures;
                        obj.longtermStimFeatures = stimFeatures;
                        disp('Loaded LongtermStim_features.mat');
                        return;
                    end
                end
                recIndices = 1:numel(obj.longtermStimData);
            else
                if isnumeric(fileSpec)
                    recIndices = fileSpec;
                elseif ischar(fileSpec) || isstring(fileSpec)
                    fileSpec = char(fileSpec);
                    idx = find(strcmp({obj.longtermStimData.name}, fileSpec), 1);
                    if isempty(idx)
                        error('Specified file "%s" not found in obj.longtermStimData.', fileSpec);
                    end
                    recIndices = idx;
                else
                    error('Invalid input: fileSpec must be an index or filename.');
                end
            end

            stimRecs = obj.longtermStimData;

            if isempty(stimRecs)
                warning('No long-term stimulation recordings available.');
                stimFeatures = struct();
                return;
            end

            recs = obj.recordings;

            % Define groups
            selectedNames = {stimRecs(recIndices).name};
            rawGroups = unique(cellfun(@(x) regexprep(x, 'S\d+R\d+.*$', ''), ...
                selectedNames, 'UniformOutput', false));
            groupNames = matlab.lang.makeValidName(rawGroups);

            % Initialize output struct
            stimFeatures = struct();
            stimFeatures.groups = groupNames;
            stimFeatures.totalDuration = struct();
            stimFeatures.RecDays = struct();

            for g = 1:numel(groupNames) % Iterate over ech stim group
                group = groupNames{g};
                % Cell array avoids dissimilar-struct assignment errors
                stimFeatures.(group) = struct([]);
                stimFeatures.totalDuration.(group) = duration(0,0,0);
                stimFeatures.RecDays.(group) = table(datetime.empty(0,1), duration.empty(0,1), ...
                    'VariableNames', {'Date','Duration'});
            end

            % Iterate over recordings
            n = numel(recIndices);
            disp('Collecting rest/stim features from long-term stim .mat files...');
            for ii = 1:n
                i = recIndices(ii);
                fname = stimRecs(i).name;
                fIdx = find(strcmp({recs.name}, fname), 1);

                if isempty(fIdx)
                    warning('Recording "%s" not found in obj.recordings. Skipping.', fname);
                    continue;
                end

                FS = recs(fIdx).SamplingRate;

                % Determine group
                fGrp = regexprep(fname, 'S\d+R\d+.*$', '');
                grpIdx = find(strcmp(fGrp, rawGroups), 1);

                if isempty(grpIdx)
                    warning('Could not assign recording "%s" to any stimulation group. Skipping.', fname);
                    continue;
                end

                group = groupNames{grpIdx};

                % Duration/date bookkeeping
                recDuration = duration(stimRecs(i).Duration);
                recDate = dateshift(datetime(stimRecs(i).date), 'start', 'day');
                stimFeatures.totalDuration.(group) = stimFeatures.totalDuration.(group) + recDuration;

                dayTable = stimFeatures.RecDays.(group);
                dayIdx = find(dayTable.Date == recDate, 1);

                if isempty(dayIdx)
                    dayTable = [dayTable; table(recDate, recDuration, ...
                        'VariableNames', {'Date','Duration'})];
                else
                    dayTable.Duration(dayIdx) = dayTable.Duration(dayIdx) + recDuration;
                end

                stimFeatures.RecDays.(group) = dayTable;

                % Basic recording entry
                recEntry = struct();
                recEntry.name = fname;
                recEntry.date = recDate;
                recEntry.duration = recDuration;
                recEntry.sourceRecordingIdx = fIdx;

                for r =1:numel(rerefMethod)
                    rrM = rerefMethod{r};
                    recEntry.(rrM) = [];
                end

                % Determine where the corresponding MAT file is stored
                if contains(lower(fname), 'baseline') || contains(lower(fname), 'rest')
                    % Baseline/rest recordings are stored with the regular recordings
                    MAT_filePath = fullfile( obj.recordingsMATFolder, strrep(fname, '.dat', '.mat'));
                else
                    % Long-term stimulation recordings are stored separately
                    MAT_filePath = fullfile( obj.longtermStimFolder, strrep(fname, '.dat', '.mat'));
                end

                if ~isfile(MAT_filePath)
                    warning('MAT file not found for "%s": %s', fname, MAT_filePath);
                    continue;
                end

                S = load(MAT_filePath);
                % Bad channels
                badChnls = [];
                try
                    badChnls = recs(fIdx).BadChannels.BadChannels;
                catch
                    warning('Could not retrieve bad channels for "%s".', fname);
                end

                % Iterate over rereferencing methods
                for rr = 1:numel(rerefMethod)
                    rrM = rerefMethod{rr};
                    if ~isfield(S, rrM)
                        warning('Missing rereference field "%s" in %s.', rrM, MAT_filePath);
                        continue;
                    end
                    rrStruct = S.(rrM);

                    if isfield(rrStruct, 'ChannelNames') && ~isempty(rrStruct.ChannelNames)
                        ch_names = rrStruct.ChannelNames;
                    elseif isfield(recs(fIdx), rrM) && ...
                            isfield(recs(fIdx).(rrM), 'ChannelNames') && ...
                            ~isempty(recs(fIdx).(rrM).ChannelNames)
                        ch_names = recs(fIdx).(rrM).ChannelNames;
                    else
                        warning('Missing ChannelNames for "%s" in "%s".', rrM, fname);
                        continue;
                    end

                    % Check if channelSpec is specified, then select only that channel,
                    % otherwise select all
                    if ~strcmpi(channelSpec, "all")
                        chNamesString = string(ch_names);
                        cleanChNames = erase(chNamesString, "Bad-");
                        chIdx = find(strcmp(cleanChNames, channelSpec), 1);

                        if isempty(chIdx) && ~strcmp(rrM, 'BP') % If channel is missing and it's not for 'BP recording  throw an error 
                            error('Specified channel "%s" is not present in recording "%s" for rereference method "%s".', channelSpec, fname, rrM);
                        elseif  isempty(chIdx) && strcmp(rrM, 'BP') % Otherwise skip iteration
                            warning('Specified channel "%s" is not present in recording "%s" for rereference method "%s".', channelSpec, fname, rrM);
                            continue
                        end

                        if startsWith(chNamesString(chIdx), "Bad-", "IgnoreCase", true)
                            warning('Specified channel "%s" is marked as bad in recording "%s" / "%s". Skipping.', channelSpec, fname, rrM);
                            continue;
                        end
                        goodChnls = chIdx;
                    else
                        if strcmp(rrM, 'BP')
                            goodChnls = 1:numel(ch_names);
                        else
                            goodChnls = 1:numel(ch_names);
                            badChnlsInRange = badChnls(badChnls >= 1 & badChnls <= numel(ch_names));
                            goodChnls(badChnlsInRange) = [];
                        end
                    end

                    recEntry.(rrM).ChannelNames = ch_names;
                    recEntry.(rrM).GoodChannels = goodChnls;

                    % ------ 1. Stimulation-period PSDs --------
                    recEntry.(rrM).StimPSD = struct();
                    recEntry.(rrM).StimPSD.trialMeanAcrossGoodChannels = [];
                    recEntry.(rrM).StimPSD.meanAcrossStimTrials = [];
                    recEntry.(rrM).StimPSD.Frequencies = [];

                    if isfield(rrStruct, 'StimPSDs') && ~isempty(rrStruct.StimPSDs)
                        StimPSDs = rrStruct.StimPSDs;
                        stimPsdTrialAvg = [];
                        for sp = 1:numel(StimPSDs)
                            thisPSD = StimPSDs{sp}; % expected [freq x channels]
                            if isempty(thisPSD)
                                continue;
                            end
                            % Verify that we are selecting only the
                            % specified channel if channelSpec was selected
                            validGood = goodChnls(goodChnls <= size(thisPSD,2));
                            if isempty(validGood)
                                continue;
                            end
                            stimPsdTrialAvg(:,end+1) = mean(thisPSD(:,validGood), 2, 'omitnan'); %#ok<AGROW>
                        end
                        recEntry.(rrM).StimPSD.trialMeanAcrossGoodChannels = stimPsdTrialAvg;
                        if ~isempty(stimPsdTrialAvg)
                            recEntry.(rrM).StimPSD.meanAcrossStimTrials = mean(stimPsdTrialAvg, 2, 'omitnan');
                        end
                        if isfield(rrStruct, 'Info') && isfield(rrStruct.Info, 'Frequencies')
                            recEntry.(rrM).StimPSD.Frequencies = rrStruct.Info.Frequencies;
                        end
                    end

                    % ------- 2. Spike summary --------
                    recEntry.(rrM).nSpikes = NaN;
                    recEntry.(rrM).normSpikeCount = NaN;

                    try
                        if ~isempty(recs(fIdx).StimInfo)
                            nonStimDurSamples = sum( ...
                                recs(fIdx).StimInfo.StimSequence.SequenceDuration( ...
                                recs(fIdx).StimInfo.StimSequence.Sequence == 0));
                            nonStimDurMin = minutes(seconds(nonStimDurSamples / FS));
                        else
                            nonStimDurMin = minutes(duration(recs(fIdx).Duration));
                        end
                        spikeDetect = recs(fIdx).(rrM).SpikeDetect;

                        if isfield(spikeDetect, 'RFOMP') && isfield(spikeDetect, 'Channelinformation')
                            % Verify that we are selecting only the specified channel in channelSpec
                            try
                                spikeChnls = spikeDetect.Channelinformation(spikeDetect.RFOMP.feature.Pred == 1);
                            catch
                                spikeChnls = spikeDetect.Channelinformation(spikeDetect.RFOMP.Pred == 1);
                            end

                            if ~strcmpi(channelSpec, "all")
                                spikeChnls = spikeChnls(spikeChnls == goodChnls);
                                nCh = 1;
                            elseif ~strcmp(rrM, 'BP')
                                spikeChnls = spikeChnls(~ismember(spikeChnls, badChnls));
                                nCh = numel(spikeDetect.Montage.ChannelNames) - numel(badChnls);
                            else
                                nCh = numel(spikeDetect.Montage.ChannelNames);
                            end

                            recEntry.(rrM).nSpikes = numel(spikeChnls);

                            if nonStimDurMin > 0 && nCh > 0
                                recEntry.(rrM).normSpikeCount = ...
                                    numel(spikeChnls) / (nonStimDurMin * nCh);
                            end
                        end
                    catch
                        warning('Could not extract spike information for "%s" / "%s".', fname, rrM);
                    end

                    % 3. Rest / baseline spectral feature summary
                    recEntry.(rrM).PSD = struct();
                    recEntry.(rrM).BandPower = struct();
                    recEntry.(rrM).BandPowerSummary = struct();

                    hasTrialSpectralFeatures = isfield(rrStruct, 'SpectralFeatures') && isstruct(rrStruct.SpectralFeatures);
                    hasDirectPSD = isfield(rrStruct, 'PSD') && ~isempty(rrStruct.PSD);
                    hasDirectBandPower = isfield(rrStruct, 'BandPower') && ...
                        isstruct(rrStruct.BandPower);

                    % Case A: baseline/passive recording, direct PSD/BandPower fields
                    if ~hasTrialSpectralFeatures && (hasDirectPSD || hasDirectBandPower)
                        % Direct PSD
                        if hasDirectPSD
                            psd = rrStruct.PSD; % expected [freq x channels]
                                                        % Verify that we are selecting only the
                            % specified channel if channelSpec was selected
                            validGood = goodChnls(goodChnls <= size(psd,2));
                            if ~isempty(validGood)
                                psdAvg = mean(psd(:,validGood), 2, 'omitnan');
                            else
                                psdAvg = [];
                            end
                            recEntry.(rrM).PSD.trialMeanAcrossGoodChannels = psdAvg;
                            recEntry.(rrM).PSD.meanAcrossTrials = psdAvg;
                            if isfield(rrStruct, 'Info') && isfield(rrStruct.Info, 'Frequencies')
                                recEntry.(rrM).PSD.Frequencies = rrStruct.Info.Frequencies;
                            elseif isfield(rrStruct, 'info') && isfield(rrStruct.info, 'Frequencies')
                                recEntry.(rrM).PSD.Frequencies = rrStruct.info.Frequencies;
                            else
                                recEntry.(rrM).PSD.Frequencies = [];
                            end
                        else
                            recEntry.(rrM).PSD.trialMeanAcrossGoodChannels = [];
                            recEntry.(rrM).PSD.meanAcrossTrials = [];
                            recEntry.(rrM).PSD.Frequencies = [];
                        end

                        % Direct BandPower
                        if hasDirectBandPower
                            bandNames = fieldnames(rrStruct.BandPower);
                            for b = 1:numel(bandNames)
                                currBand = bandNames{b};
                                bandData = rrStruct.BandPower.(currBand); % expected [time x channels] or [1 x channels]
                                if isempty(bandData)
                                    val = NaN;
                                else
                                    validGood = goodChnls(goodChnls <= size(bandData,2));
                                    if isempty(validGood)
                                        val = NaN;
                                    else
                                        val = mean(bandData(:,validGood), 'all', 'omitnan');
                                    end
                                end
                                recEntry.(rrM).BandPower.(currBand) = val;
                                recEntry.(rrM).BandPowerSummary.(currBand) = val;
                            end
                        end
                        % Done with this reref method
                        continue
                    end

                    % Case B: trial-based SpectralFeatures
                    if ~hasTrialSpectralFeatures
                        warning('Missing SpectralFeatures or direct PSD/BandPower for "%s" / "%s".', fname, rrM);
                        continue;
                    end

                    trials = fieldnames(rrStruct.SpectralFeatures);
                    isTrialField = ~cellfun(@isempty, regexp(trials, '^Trial\d+$', 'once'));
                    trials = trials(isTrialField);

                    if isempty(trials)
                        warning('No Trial# fields found in SpectralFeatures for "%s" / "%s".', fname, rrM);
                        continue;
                    end

                    trialNums = cellfun(@(x) sscanf(x, 'Trial%d'), trials);
                    [~, sortIdx] = sort(trialNums);
                    trials = trials(sortIdx);

                    psdTrialAvg = [];
                    bandPowerTrialAvg = struct();

                    for tr = 1:numel(trials)
                        trial = trials{tr};
                        trialStruct = rrStruct.SpectralFeatures.(trial);

                        % BandPower
                        if isfield(trialStruct, 'BandPower') && isstruct(trialStruct.BandPower)
                            bandNames = fieldnames(trialStruct.BandPower);
                            for b = 1:numel(bandNames)
                                currBand = bandNames{b};
                                bandData = trialStruct.BandPower.(currBand);
                                if isempty(bandData)
                                    val = NaN;
                                else
                                   % Verify that we are selecting only the
                            % specified channel if channelSpec was selected
                                    validGood = goodChnls(goodChnls <= size(bandData,2));
                                    if isempty(validGood)
                                        val = NaN;
                                    else
                                        val = mean(bandData(:,validGood), 'all', 'omitnan');
                                    end
                                end
                                bandPowerTrialAvg.(currBand)(tr,1) = val;
                            end
                        end

                        % PSD
                        if isfield(trialStruct, 'PSD') && ~isempty(trialStruct.PSD)
                            psd = trialStruct.PSD; % expected [freq x channels]
                            % Verify that we are selecting only the
                            % specified channel if channelSpec was selected
                            validGood = goodChnls(goodChnls <= size(psd,2));
                            if ~isempty(validGood)
                                psdTrialAvg(:,tr) = mean(psd(:,validGood), 2, 'omitnan'); %#ok<AGROW>
                            end
                        end
                    end

                    recEntry.(rrM).PSD.trialMeanAcrossGoodChannels = psdTrialAvg;
                    if ~isempty(psdTrialAvg)
                        recEntry.(rrM).PSD.meanAcrossTrials = mean(psdTrialAvg, 2, 'omitnan');
                    else
                        recEntry.(rrM).PSD.meanAcrossTrials = [];
                    end

                    try
                        firstTrial = rrStruct.SpectralFeatures.(trials{1});
                        if isfield(firstTrial, 'info') && isfield(firstTrial.info, 'Frequencies')
                            recEntry.(rrM).PSD.Frequencies = firstTrial.info.Frequencies;
                        elseif isfield(firstTrial, 'Info') && isfield(firstTrial.Info, 'Frequencies')
                            recEntry.(rrM).PSD.Frequencies = firstTrial.Info.Frequencies;
                        else
                            recEntry.(rrM).PSD.Frequencies = [];
                        end
                    catch
                        recEntry.(rrM).PSD.Frequencies = [];
                    end

                    recEntry.(rrM).BandPower = bandPowerTrialAvg;
                    bandNames = fieldnames(bandPowerTrialAvg);
                    for b = 1:numel(bandNames)
                        currBand = bandNames{b};
                        recEntry.(rrM).BandPowerSummary.(currBand) = ...
                            mean(bandPowerTrialAvg.(currBand), 'omitnan');
                    end
                end
                % Store recording entry
                if isempty(stimFeatures.(group))
                    stimFeatures.(group) = recEntry;
                else
                    stimFeatures.(group)(end+1) = recEntry;
                end
                progressbar(ii, n)
            end


            % Save only when full dataset was processed
            if numel(recIndices) == numel(obj.longtermStimData)
                if isfile(fpath)
                    choice = questdlg('Overwrite existing LongtermStim_features.mat?', ...
                        'Confirm Overwrite', 'Yes', 'No', 'No');
                    if ~strcmp(choice, 'Yes')
                        disp('User cancelled saving.');
                        obj.longtermStimFeatures = stimFeatures;
                        return;
                    end
                end
                % If fileSpec or channel spec was specified add it to the output file name
                fname = 'stimFeatures';
                if strcmp(channelSpec,'all')
                    fname = [fname '_' channelSpec];
                end
                save(fpath, fname, '-v7.3');
                disp(['Saved to ', fpath]);
            end
            obj.longtermStimFeatures = stimFeatures;
        end


        %% ~~~~~~~~~~~~ Visualization functions  ~~~~~~~~~~~~~~~~~~~~~~~~~~
        function [plotData, plotDataStd] = plotAveragePSD(obj, PSDdataStruct, excludeBadCh, statType, refType)
            % plotAveragePSD - Plots average or median PSD across non-stim recordings
            %
            % Inputs:
            %   PSDdataStruct - optional struct with PSDs (name, date, rerefs, etc.)
            %   excludeBadCh  - optional logical (default = false)
            %   statType      - 'mean' (default) or 'median'
            %   refType       - 'all' (default) or one of rerefTypes
            % Outputs (optional):
            %   plotData      - struct with averaged PSDs
            %   plotDataStd   - struct with standard deviation across recordings

            arguments
                obj
                PSDdataStruct = []
                excludeBadCh logical = false
                statType (1,:) char {mustBeMember(statType, {'mean','median'})} = 'mean'
                refType (1,:) = 'all'
            end

            % Load PSD struct if empty
            if isempty(PSDdataStruct)
                PSDmatFile = fullfile(obj.recordingsMATFolder, 'PSDs.mat');
                if isfile(PSDmatFile)
                    choice = questdlg('Load existing PSDs.mat or compute now?', ...
                        'Load PSDs', 'Load', 'Compute', 'Load');
                    if strcmp(choice, 'Load')
                        data = load(PSDmatFile, 'psdStructs');
                        PSDdataStruct = data.psdStructs;
                        disp('Loaded PSDs.mat');
                    else
                        PSDdataStruct = obj.collectPSDs();
                    end
                else
                    PSDdataStruct = obj.collectPSDs();
                end
            end

            if isempty(PSDdataStruct)
                warning('No PSD data to process.');
                return;
            end

            rerefTypes = obj.rerefMethods;
            if ~strcmp(refType, 'all') && ~ismember(refType, rerefTypes)
                error('Invalid refType "%s". Must be one of: %s', refType, strjoin(rerefTypes, ', '));
            end

            % Filter non-stimulation recordings
            stimMask = [PSDdataStruct.Stim] == 1;
            nonStimPSDs = PSDdataStruct(~stimMask);
            if isempty(nonStimPSDs)
                warning('No non-stimulation PSDs found.');
                return;
            end

            % Initialize plot data
            plotData = struct();
            plotDataStd = struct();
            freqs = nonStimPSDs(1).Frequencies;

            % Loop over reref types
            targetRefs = rerefTypes;
            if ~strcmp(refType, 'all')
                targetRefs = {refType};
            end

            for r = 1:numel(targetRefs)
                ref = targetRefs{r};
                PSDs = [];

                for i = 1:numel(nonStimPSDs)
                    psdMatrix = nonStimPSDs(i).(ref);  % [freq x ch]
                    if isempty(psdMatrix)
                        continue;
                    end

                    % Optionally exclude bad channels
                    if excludeBadCh  && ~strcmp(ref, 'BP')
                        % Match by filename in obj.recordings
                        idx = find(strcmp({obj.recordings.name}, nonStimPSDs(i).name), 1);
                        if isempty(idx)
                            warning('Recording "%s" not found in obj.recordings. Skipping.', nonStimPSDs(i).name);
                            continue;
                        end
                        badCh = obj.recordings(idx).BadChannels.BadChannels;
                        goodIdx = setdiff(1:size(psdMatrix, 2), badCh);
                        psdMatrix = psdMatrix(:, goodIdx);
                    end

                    % Average over channels
                    PSDs(:, end+1) = mean(psdMatrix, 2); %#ok<AGROW>
                end

                % Compute overall stat
                switch statType
                    case 'mean'
                        plotData.(ref) = mean(PSDs, 2);
                    case 'median'
                        plotData.(ref) = median(PSDs, 2);
                end
                plotDataStd.(ref) = std(PSDs, 0, 2);  % std over recordings
            end

            % Plot
            figure('Name', 'Average PSD', 'Position',[249 127 1259 816]);
            hold on;
            colors = lines(numel(targetRefs));
            for r = 1:numel(targetRefs)
                ref = targetRefs{r};
                if isfield(plotData, ref)
                    plot(freqs, pow2db(plotData.(ref)), 'LineWidth', 2, 'Color', colors(r,:));
                end
            end
            hold off;
            grid on;
            xlabel('Frequency (Hz)');
            ylabel('Power [dB]');
            title(sprintf('Average (%s) PSD Across All Non-Stim Recordings', statType));
            legend(targetRefs, 'Location', 'best');
            fontsize(16,'points')
        end

        function plotData = plotRecordingCapabilities(obj, impTable, PSDdataStruct, refType)
            % plotRecordingCapabilities - Combines impedance, PSD, and RMS across days
            %
            % Inputs:
            %   impTable       - optional, uses summarizeImpedance if not provided
            %   PSDdataStruct  - optional, loads or computes using collectPSDs
            %   refType        - 'Native' (default), 'CAR', 'leadCAR', 'leadSpecific'

            arguments
                obj
                impTable = []
                PSDdataStruct = []
                refType (1,:) char {mustBeMember(refType, {'Native','CAR','leadCAR','leadSpecific'})} = 'Native'
            end

            % --- Load impedance summary if not provided
            if isempty(impTable)
                impTable = obj.summarizeImpedance();
            end

            if isempty(impTable)
                warning('No impedance data available.');
                plotData = struct([]);
                return;
            end

            % --- Load PSD struct if not provided
            if isempty(PSDdataStruct)
                PSDmatFile = fullfile(obj.recordingsMATFolder, 'PSDs.mat');
                if isfile(PSDmatFile)
                    choice = questdlg('Load existing PSDs.mat or compute now?', ...
                        'Load PSDs', 'Load', 'Compute', 'Load');
                    if strcmp(choice, 'Load')
                        data = load(PSDmatFile, 'psdStructs');
                        PSDdataStruct = data.psdStructs;
                        disp('Loaded PSDs.mat');
                    else
                        PSDdataStruct = obj.collectPSDs();
                    end
                else
                    PSDdataStruct = obj.collectPSDs();
                end
            end

            % --- Filter non-stim PSDs ---
            PSDdataStruct = PSDdataStruct([PSDdataStruct.Stim] == 0);

            % --- Extract unique days from impedance table ---
            uniqueDates = impTable.Date;

            % --- Extract frequency axis from first valid recording ---
            freqs = [];
            for i = 1:numel(PSDdataStruct)
                if isfield(PSDdataStruct, 'Frequencies') 
                    freqs = PSDdataStruct(i).Frequencies;
                    break;
                end
            end
            if isempty(freqs)
                error('Frequencies not found in SpectralStats.Frequencies for ref %s.', refType);
            end

            plotData = struct();
            nonvalid_days = [];

            for i = 1:numel(uniqueDates)
                thisDate = uniqueDates(i);
                plotData(i).Date = thisDate;

                % Step 1: Impedance
                impRow = impTable(impTable.Date == thisDate, :);
                if isempty(impRow)
                    plotData(i).Impedance = nan(1, 32);
                else
                    chImps = zeros(1, 32);
                    for ch = 1:32
                        colName = sprintf('Impedance_Ch%d', ch);
                        if ismember(colName, impTable.Properties.VariableNames)
                            chImps(ch) = impRow.(colName);
                        else
                            chImps(ch) = NaN;
                        end
                    end
                    plotData(i).Impedance = chImps;
                end

                % Step 2: Extract matching PSDs (same date)
                dateStrings = string({PSDdataStruct.date});
                psdDates = dateshift(datetime(dateStrings, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss'), 'start', 'day');
                matchMask = psdDates == thisDate;
                matchedPSDs = PSDdataStruct(matchMask);

                if isempty(matchedPSDs)
                    plotData(i).PSD = [];
                    plotData(i).RMS = [];
                    nonvalid_days = [nonvalid_days, i];
                    continue;
                end

                % Collect PSDs: [freq x ch x recordings]
                for f = 1:numel(matchedPSDs)
                    thisPSD = matchedPSDs(f).(refType);
                    if f == 1
                        [nFreq, nCh] = size(thisPSD);
                        PSDcube = nan(nFreq, nCh, numel(matchedPSDs));
                    end
                    PSDcube(:, :, f) = thisPSD;
                end
                plotData(i).PSD = mean(PSDcube, 3);  % avg over recordings

                % Step 3: Extract RMS from obj.recordings
                matchedNames = {matchedPSDs.name};
                matchedRMS = nan(numel(matchedNames), 32);
                for r = 1:numel(matchedNames)
                    idx = find(strcmp({obj.recordings.name}, matchedNames{r}), 1);
                    if ~isempty(idx) && isfield(obj.recordings(idx), refType) && ...
                            isfield(obj.recordings(idx).(refType), 'SpectralStats')
                        rms = obj.recordings(idx).(refType).SpectralStats.TotalRMS;
                        matchedRMS(r, :) = rms;
                    end
                end
                plotData(i).RMS = mean(matchedRMS, 1, 'omitnan');
            end

            % Drop days without PSDs (for example when only stim recordings were done)
            plotData(nonvalid_days) = [];

            % Plot recording capabilities development for each channel
            channel_ReCap_browser(plotData,freqs)

            % Plot summarized recording capabilities for each lead
            plotReCap(plotData)
        end

        function plotStats(obj, refType, statType)
            % plotSpectralStats - Visualizes the development of a stat across channels and time
            %
            % Inputs:
            %   refType  - rereferencing method (default = 'Native')
            %   statType - one of:
            %       'Average', 'Median', 'STD', 'Variance', 'RMS',
            %       'Skewness', 'Kurtosis', 'ZeroCrossRate',
            %       'TotalRMS', 'DeltaRMS', 'ThetaRMS', 'AlphaRMS', 'BetaRMS', 'GammaRMS'

            arguments
                obj
                refType (1,:) char {mustBeMember(refType, ...
                    {'Native','CAR','leadCAR','leadSpecific'})} = 'Native'
                statType (1,:) char {mustBeMember(statType, ...
                    {'Average','Median','STD','Variance','RMS','Skewness','Kurtosis','ZeroCrossRate', ...
                    'TotalRMS','DeltaRMS','ThetaRMS','AlphaRMS','BetaRMS','GammaRMS'})} = 'TotalRMS'
            end

            % Decide source field: Stats or SpectralStats
            statsFields = {'Average','Median','STD','Variance','RMS','Skewness','Kurtosis','ZeroCrossRate'};
            isStatsField = ismember(statType, statsFields);

            if isempty(obj.recordings)
                warning('No recordings available.');
                return;
            end

            % --- Initialize ---
            recs = obj.recordings;
            recs = recs(~arrayfun(@(r) isfield(r, 'StimInfo') && ~isempty(r.StimInfo), recs)); % Exclude stim recordings
            nCh = 32;

            % Convert and clean dates
            recDates = datetime({recs.date}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
            recDates = dateshift(recDates, 'start', 'day');
            if ~isempty(obj.excludedDates)
                mask = ~ismember(recDates, obj.excludedDates);
                recs = recs(mask);
                recDates = recDates(mask);
            end

            uniqueDates = unique(recDates);
            nDays = numel(uniqueDates);
            plotData = struct('date', [], 'stat', []);

            % --- Loop through days ---
            for d = 1:nDays             
                day = uniqueDates(d);
                plotData(d).date = day;

                dayMask = recDates == day;
                recsOnDay = recs(dayMask);
                % stack = [];

                stack = nan(numel(recsOnDay), 32);
                for r = 1:numel(recsOnDay)
                    try
                        if isStatsField
                            stack(r, :) = recsOnDay(r).(refType).Stats.(statType);
                        else
                            stack(r, :) = recsOnDay(r).(refType).SpectralStats.(statType);
                        end
                    catch
                        continue;
                    end
                end
                if ~isempty(stack)
                    plotData(d).stat = mean(stack, 1, 'omitnan');
                else
                    plotData(d).stat = nan(1, nCh);
                end
            end

            % Preallocate matrices
            dataMat = nan(nCh, nDays);       % [channel x day]

            for i = 1:nDays
                dataMat(:, i) = plotData(i).stat(:);
            end

            % --- Plot ---
            figure('Name', ['Spectral Stat - ', statType, ' (', refType, ')'],'Position', [100 100 1500 800]);
            imagesc(uniqueDates, 1:nCh, dataMat);
            % set(gca, 'YDir', 'normal');
            colormap turbo;
            c = colorbar; c.Label.String = statType;
            c.Location='northoutside';%'eastoutside'; %'northoutside';
            xlabel('Recording Date');
            ylabel('Channel'); yticks(1:32); 
            title(['Timecourse of ', statType, ' - ', refType]);
            % datetick('x', 'mmm-yyyy', 'keeplimits');
            fontsize(14, 'points');
            axis tight;
        end

        function plotData = plotLongtermStim(obj, refType, plotGroup)
            % plotLongtermStim Visualize long-term stimulation feature summaries.
            %
            % Inputs
            %   refType   - rereferencing type, default = 'Native'
            %   plotGroup - optional cell array of group names to plot
            arguments
                obj
                refType (1,:) char {mustBeMember(refType, ...
                    {'Native','CAR','leadCAR','leadSpecific','BP'})} = 'Native'
                plotGroup cell = {}
            end

            % Initialize output variable
            plotData = struct();

            % Ensure collected features exist
            if ~isprop(obj, 'longtermStimFeatures') || isempty(obj.longtermStimFeatures) || ...
                    ~isfield(obj.longtermStimFeatures, 'groups') || isempty(obj.longtermStimFeatures.groups)
                [obj, stimFeatures] = obj.collectLongtermStimFeatures();
            else
                stimFeatures = obj.longtermStimFeatures;
            end
            grps = stimFeatures.groups;

            % Validate selected groups
            if ~isempty(plotGroup)
                plotGroup = matlab.lang.makeValidName(plotGroup);
                invalidGroups = setdiff(plotGroup, grps);
                if ~isempty(invalidGroups)
                    error('Invalid plotGroup(s): %s', strjoin(invalidGroups, ', '));
                end
                grps = plotGroup;
            end

            % ----- Identify baseline/rest groups -----
            baselineKeywords = {'baseline','bsln','rest'};
            isBaselineGroup = false(size(grps));
            for g = 1:numel(grps)
                grpLower = lower(grps{g});
                isBaselineGroup(g) = any(cellfun(@(x) contains(grpLower, x), baselineKeywords));
            end
            baselineGroups = grps(isBaselineGroup);
            stimGroups = grps(~isBaselineGroup);

            % ----- Extract baseline features -----
            baselinePSDList = {};
            baselineSpikeCounts = [];
            baselineBandPower = struct();
            baselineFreqs = [];

            for g = 1:numel(baselineGroups)
                grp = baselineGroups{g};
                if ~isfield(stimFeatures, grp) || isempty(stimFeatures.(grp))
                    continue
                end

                grpStruct = stimFeatures.(grp);

                for i = 1:numel(grpStruct)
                    if ~isfield(grpStruct(i), refType)
                        continue
                    end
                    rr = grpStruct(i).(refType);
                    if isfield(rr, 'PSD') && isfield(rr.PSD, 'meanAcrossTrials') && ~isempty(rr.PSD.meanAcrossTrials)
                        baselinePSDList{end+1} = rr.PSD.meanAcrossTrials(:); %#ok<AGROW>
                        if isempty(baselineFreqs) && isfield(rr.PSD, 'Frequencies')
                            baselineFreqs = rr.PSD.Frequencies(:);
                        end
                    end

                    if isfield(rr, 'normSpikeCount') && ~isempty(rr.normSpikeCount)
                        baselineSpikeCounts(end+1,1) = rr.normSpikeCount; %#ok<AGROW>
                    end

                    if isfield(rr, 'BandPowerSummary')
                        bandNames = fieldnames(rr.BandPowerSummary);
                        for b = 1:numel(bandNames)
                            bandName = bandNames{b};
                            if ~isfield(baselineBandPower, bandName)
                                baselineBandPower.(bandName) = [];
                            end
                            baselineBandPower.(bandName)(end+1,1) = ...
                                rr.BandPowerSummary.(bandName);
                        end
                    end
                end
            end

            if ~isempty(baselinePSDList)
                avgBaselinePSD = mean(cat(2, baselinePSDList{:}), 2, 'omitnan');
            else
                avgBaselinePSD = [];
            end

            % Prepare summary collection
            summaryData = struct();
            summaryData.SpikeCount.values = [];
            summaryData.SpikeCount.group = {};

            bandSummary = struct();

            % Iterate over stimulation groups
            for g = 1:numel(stimGroups)
                grp = stimGroups{g};
                if ~isfield(stimFeatures, grp) || isempty(stimFeatures.(grp))
                    warning('Group "%s" is empty or missing.', grp);
                    continue
                end
                grpStruct = stimFeatures.(grp);

                restList = {};
                stimList = {};
                recDates = datetime.empty(0,1);
                recDurations = duration.empty(0,1);
                freqs = [];

                % Extract recording-level values
                for i = 1:numel(grpStruct)
                    recDates(i,1) = grpStruct(i).date;
                    recDurations(i,1) = grpStruct(i).duration;

                    if ~isfield(grpStruct(i), refType)
                        continue
                    end

                    rr = grpStruct(i).(refType);

                    % Rest PSD
                    if isfield(rr, 'PSD') && isfield(rr.PSD, 'meanAcrossTrials') && ~isempty(rr.PSD.meanAcrossTrials)
                        restList{end+1} = rr.PSD.meanAcrossTrials(:); %#ok<AGROW>
                        if isempty(freqs) && isfield(rr.PSD, 'Frequencies')
                            freqs = rr.PSD.Frequencies(:);
                        end
                    end

                    % Stim PSD
                    if isfield(rr, 'StimPSD') && isfield(rr.StimPSD, 'meanAcrossStimTrials') && ~isempty(rr.StimPSD.meanAcrossStimTrials)
                        stimList{end+1} = rr.StimPSD.meanAcrossStimTrials(:); %#ok<AGROW>
                        if isempty(freqs) && isfield(rr.StimPSD, 'Frequencies')
                            freqs = rr.StimPSD.Frequencies(:);
                        end
                    end

                    % Normalized spike count
                    if isfield(rr, 'normSpikeCount') && ~isempty(rr.normSpikeCount)
                        summaryData.SpikeCount.values(end+1,1) = rr.normSpikeCount; %#ok<AGROW>
                        summaryData.SpikeCount.group{end+1,1} = grp; %#ok<AGROW>
                    end

                    % Band power summary
                    if isfield(rr, 'BandPowerSummary')
                        bandNames = fieldnames(rr.BandPowerSummary);
                        for b = 1:numel(bandNames)
                            bandName = bandNames{b};
                            if ~isfield(bandSummary, bandName)
                                bandSummary.(bandName).values = [];
                                bandSummary.(bandName).group = {};
                            end
                            bandSummary.(bandName).values(end+1,1) = rr.BandPowerSummary.(bandName);
                            bandSummary.(bandName).group{end+1,1} = grp;
                        end
                    end
                end

                % Daily duration table
                if isfield(stimFeatures, 'RecDays') && isfield(stimFeatures.RecDays, grp)
                    recDays = stimFeatures.RecDays.(grp);
                else
                    recDays = table(recDates, recDurations, 'VariableNames', {'Date','Duration'});
                end

                % Average PSDs
                meanRestPSD = [];
                meanStimPSD = [];
                if ~isempty(restList)
                    meanRestPSD = mean(cat(2, restList{:}), 2, 'omitnan');
                end
                if ~isempty(stimList)
                    meanStimPSD = mean(cat(2, stimList{:}), 2, 'omitnan');
                end

                if isempty(freqs)
                   freqs = 0:(obj.samplingRate/obj.spectralParameters.NFFT):(obj.samplingRate/2);
                   freqs = freqs';
                end

                % ----- Figure 1: overview for current group -------
                figure('Name', sprintf('Long-Term Stimulation - %s', grp), ...
                    'Position', [200 100 1400 850]);

                tl = tiledlayout(2,2,  'TileSpacing','compact', 'Padding','compact');

                sgtitle(sprintf('%s | %s', grp, refType), 'Interpreter','none', ...
                    'FontSize', 18, 'FontWeight','bold');

                % Plot 1: daily recording durations
                ax1 = nexttile(tl,1);
                if ~isempty(recDays) && height(recDays) > 0
                    % bar(ax1, recDays.Date, minutes(recDays.Duration));
                    bar(ax1, 1:numel(recDays.Date), minutes(recDays.Duration));
                end
                xlabel(ax1, 'Recording day [-]');
                ylabel(ax1, 'Duration [min]');
                title(ax1, 'Recording duration');
                grid(ax1, 'on');

                % Plot 2: stim PSD evolution
                ax2 = nexttile(tl,3);
                if ~isempty(stimList) && ~isempty(freqs)
                    stimMat = cat(2, stimList{:});
                    imagesc(ax2, 1:size(stimMat,2), freqs, pow2db(stimMat + eps));
                    set(ax2, 'YDir', 'normal');
                    colormap(ax2, turbo);
                    cb = colorbar(ax2);
                    cb.Label.String = 'Power [dB]';
                    xlabel(ax2, 'Recording index');
                    ylabel(ax2, 'Frequency [Hz]');
                    title(ax2, 'Stim-period PSD evolution');
                else
                    title(ax2, 'No stim-period PSD available');
                end

                % Plot 3: baseline vs rest PSD
                ax3 = nexttile(tl,2);
                hold(ax3, 'on');

                % Plot rest PSDs and their average
                if ~isempty(avgBaselinePSD) && ~isempty(freqs)
                    % plot(ax3, freqs, pow2db(cat(2, baselinePSDList{:}) + eps), ...
                    %     'Color', [0 0.4470 0.7410 0.2], 'HandleVisibility', 'off')
                    % hold on
                    plot(ax3, freqs, pow2db(avgBaselinePSD + eps), 'Color', [0 0.4470 0.7410], ...
                        'LineWidth', 2, 'DisplayName', 'Baseline/rest group');
                end

                if ~isempty(meanRestPSD) && ~isempty(freqs)
                    % plot(ax3, freqs, pow2db(cat(2, restList{:}) + eps), ...
                    %     'Color', [0.8500    0.3250    0.0980 0.2], 'HandleVisibility', 'off')
                    % hold on;
                    plot(ax3, freqs, pow2db(meanRestPSD + eps), 'Color', [0.8500    0.3250    0.098], ...
                        'LineWidth', 2,'DisplayName', 'Rest periods in stim recordings');
                end

                grid(ax3, 'on');
                xlabel(ax3, 'Frequency [Hz]');
                ylabel(ax3, 'Power [dB]');
                title(ax3, 'Rest PSD comparison');
                legend(ax3, 'Location', 'best');
                xlim(ax3, 'tight');

                % Plot 4: stim vs rest PSD
                ax4 = nexttile(tl,4);
                hold(ax4, 'on');

                if ~isempty(meanStimPSD) && ~isempty(freqs)
                    plot(ax4, freqs, pow2db(meanStimPSD + eps), ...
                        'LineWidth', 2, ...
                        'DisplayName', 'Stim periods');
                end

                if ~isempty(meanRestPSD) && ~isempty(freqs)
                    plot(ax4, freqs, pow2db(meanRestPSD + eps), ...
                        'LineWidth', 2, ...
                        'DisplayName', 'Rest periods');
                end

                grid(ax4, 'on');
                xlabel(ax4, 'Frequency [Hz]');
                ylabel(ax4, 'Power [dB]');
                title(ax4, 'Stim vs rest PSD');
                legend(ax4, 'Location', 'best');
                xlim(ax4, 'tight');
                fontsize(gcf, 14, 'points');
            end

            %% Figure 2: summary boxplots
            
            % ---- Spike counts
            figure('Name', 'Long-Term Stimulation Summary', 'Position', [250 150 1000 700]);
            groupOrder = [baselineGroups(:); stimGroups(:)];
            categs = categorical(...
                [repelem(baselineGroups,numel(baselineSpikeCounts))'; ...
                summaryData.SpikeCount.group], groupOrder, 'Ordinal',true);
            cats = categories(categs);
            vals = [baselineSpikeCounts ; summaryData.SpikeCount.values];
            
            boxchart(categs, vals); hold on
            for k = 1:numel(cats)
                meanSpikeCount(k) = mean(vals(categs == cats{k}), 'omitnan');
            end
            % plot(1:numel(cats), meanSpikeCount, '-o', 'LineWidth', 1.5);
            ylabel('Spikes / (channel × min)');
            title(sprintf('Normalized spike count | %s', refType), ...
                'Interpreter','none', 'FontSize', 18, 'FontWeight','bold');
            xaxisproperties= get(gca, 'XAxis');
            xaxisproperties.TickLabelInterpreter = 'none';
            grid('on');

            % ----- Band power summaries Figure ------
            figure('Name', 'Long-Term Stimulation Summary', 'Position', [250 150 1400 700]);
            bandNames = fieldnames(bandSummary);
            nFeatures = numel(bandNames);
            tl = tiledlayout(1, nFeatures, 'TileSpacing','compact',  'Padding','compact');

            sgtitle(sprintf('Long-term stimulation summary | %s', refType), ...
                'Interpreter','none', 'FontSize', 18, 'FontWeight','bold');
            
            for b = 1:numel(bandNames)
                bandName = bandNames{b};
                ax = nexttile(tl,b);

                vals = [baselineBandPower.(bandName); bandSummary.(bandName).values];
                cats = categorical( ...
                      [repelem(baselineGroups,numel(baselineBandPower.(bandName)))'; ...
                      bandSummary.(bandName).group], groupOrder, 'Ordinal',true);

                if ~isempty(vals)
                    boxchart(ax, cats, vals);
                    ylabel(ax, 'Band power [\muV^2]');
                    title(ax, bandName, 'Interpreter','none');
                    ax.TickLabelInterpreter = 'none';
                    grid(ax, 'on');
                    ylim(prctile(vals,[2.5 97.5]))
                else
                    title(ax, sprintf('%s unavailable', bandName), 'Interpreter','none');
                end
            end

            fontsize(gcf, 14, 'points');

            % ------ Populate output variable ------
            plotData.BandPower.Baseline = baselineBandPower;
            plotData.BandPower.Stim = bandSummary;
            plotData.SpikeCount.cats = categorical(...
                [repelem(baselineGroups,numel(baselineSpikeCounts))'; ...
                summaryData.SpikeCount.group], groupOrder, 'Ordinal',true); 
            plotData.SpikeCount.vals = [baselineSpikeCounts ; summaryData.SpikeCount.values];;
            plotData.RecDays = recDays;
        end

        function plotRecordingTimes(obj, dataStruct)
            % plotRecordingTimes - Plots a histogram of recording start times over 24 hours
            %
            % Inputs:
            %   dataStruct - either obj.recordings or obj.longtermStimData
            %
            % The method extracts the 'date' field and plots the time-of-day distribution

            if nargin < 2 || isempty(dataStruct)
                error('You must provide obj.recordings or obj.longtermStimData as input.');
            end

            % Validate presence of 'date' field
            if ~all(isfield(dataStruct, 'date'))
                error('Input structure does not contain a ''date'' field.');
            end

            % Extract datetimes
            try
                datetimes = datetime({dataStruct.date}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
            catch
                error('Failed to parse date strings. Ensure date format is ''dd-MMM-yyyy HH:mm:ss''.');
            end

            % Extract time of day
            timeOfDay = timeofday(datetimes);  % duration format

            % Convert to hours as decimal (e.g. 13.5 = 1:30 PM)
            hoursDecimal = hours(timeOfDay);

            % Plot histogram
            figure;
            histogram(hoursDecimal, 24, ...
                'BinEdges', 0:1:24, ...
                'FaceColor', [0.2 0.5 0.8], ...
                'EdgeColor', 'k');

            title('Recording Times Over 24-Hour Day');
            xlabel('Hour of Day');
            ylabel('Number of Recordings');
            xticks(0:1:24);
            xlim([0 24]);
            grid on;
            fontsize(16,"points")
        end
    
        function plotPSDoverTime(obj, PSDdataStruct, N, excludeBadCh, statType, refType, cmap)
            % plotPSDoverTime - Visualizes PSD evolution across time
            %                  by plotting N averaged PSDs, color-coded by time.
            %
            % Inputs:
            %   PSDdataStruct - optional struct with PSDs (if empty, will load)
            %   N             - number of time bins (default = 10)
            %   excludeBadCh  - logical, whether to exclude bad channels (default = false)
            %   statType      - 'mean' (default) or 'median'
            %   refType       - 'all' (default) or specific rereferencing method
            %   cmap          - colormap (default = winter(N))

            % TO DO: group the PSD lines by leads and mark know bad
            % channels by red
            arguments
                obj
                PSDdataStruct = []
                N (1,1) double = 10
                excludeBadCh logical = false
                statType (1,:) char {mustBeMember(statType, {'mean','median'})} = 'mean'
                refType (1,:) char = 'Native' % {mustBeMember(refType, {'Native','CAR','BP','leadCAR','leadSpecific'})}
                cmap = []  % default to winter(N) below if empty
            end

            % Load PSD struct if empty
            if isempty(PSDdataStruct)
                PSDmatFile = fullfile(obj.recordingsMATFolder, 'PSDs.mat');
                if isfile(PSDmatFile)
                    choice = questdlg('Load existing PSDs.mat or compute now?', ...
                        'Load PSDs', 'Load', 'Compute', 'Load');
                    if strcmp(choice, 'Load')
                        data = load(PSDmatFile, 'psdStructs');
                        PSDdataStruct = data.psdStructs;
                        disp('Loaded PSDs.mat');
                    else
                        PSDdataStruct = obj.collectPSDs();
                    end
                else
                    PSDdataStruct = obj.collectPSDs();
                end
            end

            if isempty(PSDdataStruct)
                warning('No PSD data to process.');
                return;
            end

            rerefTypes = obj.rerefMethods;
            if ~strcmp(refType, 'all') && ~ismember(refType, rerefTypes)
                error('Invalid refType "%s". Must be one of: %s', refType, strjoin(rerefTypes, ', '));
            elseif strcmp(refType,'BP') % Throw an error when selecting reftype BP
                error('Unable to perform assignment for refType "BP" due to the various indexing of channels')
            end


            % Filter non-stimulation recordings
            stimMask = [PSDdataStruct.Stim] == 1;
            nonStimPSDs = PSDdataStruct(~stimMask);
            if isempty(nonStimPSDs)
                warning('No non-stimulation PSDs found.');
                return;
            end

            % Retrieve leads
            leadMap = obj.leadAssignment;
            leads = keys(leadMap);
            nLeads = numel(leads);

            % Convert dates to datetime and align to day
            recDates = dateshift( ...
                datetime({nonStimPSDs.date}, 'InputFormat', 'd-MMM-y HH:mm:ss'), ...
                'start', 'day');
            uniqueDates = unique(recDates, 'sorted')';
            if ~isempty(obj.excludedDates)
                [~,excludeMask] = intersect(uniqueDates,obj.excludedDates);
                uniqueDates(excludeMask) = [];
            end

            % Split into N bins
            if N == 1
                dateIdxs = numel(uniqueDates); % Plot the latest PSDs
            else
                dateIdxs = round(linspace(1, numel(uniqueDates), N));
            end

            % Initialize plot data
            plotData = struct();
            freqs = nonStimPSDs(1).Frequencies;

            for i = 1:numel(dateIdxs)
                selectedDays = find(recDates == uniqueDates(dateIdxs(i)));
                plotData(i).date = uniqueDates(dateIdxs(i));

                for l = 1:nLeads
                    leadName = leads{l};
                    leadNameFld = matlab.lang.makeValidName(leads{l});
                    chIdxsOriginal = leadMap(leadName).indices;
                    avgDayPSD = nan(numel(freqs), numel(selectedDays));

                    for j = 1:numel(selectedDays)
                        psdIdx = selectedDays(j);
                        psdMatrix = nonStimPSDs(psdIdx).(refType);  % [freq x ch]
                        if isempty(psdMatrix)
                            continue;
                        end
                        chIdxs = chIdxsOriginal;

                        if excludeBadCh && ~strcmp(refType, 'BP')
                            idx = find(strcmp({obj.recordings.name}, nonStimPSDs(psdIdx).name), 1);
                            if isempty(idx)
                                warning('Recording "%s" not found in obj.recordings. Skipping.', ...
                                    nonStimPSDs(psdIdx).name);
                                continue;
                            end
                            badCh = obj.recordings(idx).BadChannels.BadChannels;
                            % Remove bad channels from this lead before indexing
                            chIdxs = setdiff(chIdxs, badCh);
                        end

                        % Prevent indexing outside the PSD matrix
                        chIdxs = chIdxs(chIdxs >= 1 & chIdxs <= size(psdMatrix, 2));
                        if isempty(chIdxs)
                            plotData(i).(leadNameFld).PSD = nan(numel(freqs));
                            continue;
                        end

                        leadPSD = psdMatrix(:, chIdxs);
                        % Average over channels within lead
                        if strcmp(statType, 'median')
                            avgDayPSD(:,j) = median(leadPSD, 2, 'omitnan');
                        else
                            avgDayPSD(:,j) = mean(leadPSD, 2, 'omitnan');
                        end

                    end
                    % Average across recordings from the selected day
                    if strcmp(statType, 'median')
                        plotData(i).(leadNameFld).PSD = median(avgDayPSD, 2, 'omitnan');
                    else
                        plotData(i).(leadNameFld).PSD = mean(avgDayPSD, 2, 'omitnan');
                    end

                    plotData(i).(leadNameFld).channels = chIdxsOriginal;
                end
            end

            % Default colormap
            if isempty(cmap)
                colors = flipud(winter(N));
            else
                colors = cmap;
            end

            % Plot
            f = figure('Name', 'PSD Evolution', 'Position',[100 100 1500 800]);
            tl = tiledlayout(f,'flow', 'TileSpacing', 'compact', 'Padding', 'compact');

            for l = 1:nLeads
                leadNameFld = matlab.lang.makeValidName(leads{l});
                ax = nexttile(tl);
                hold(ax, 'on');
                for d = 1:numel(plotData)
                    if d==length(plotData) % Plot the last line more thick
                        plot(ax,freqs, ...
                            smoothdata(pow2db(plotData(d).(leadNameFld).PSD),'movmean',10), ...
                            'LineWidth', 3, ...
                            'Color', colors(d,:), ...
                            'DisplayName', string(plotData(d).date));
                    else
                        plot(ax,freqs, ...
                            smoothdata(pow2db(plotData(d).(leadNameFld).PSD),'movmean',10), ...
                            'LineWidth', 1.2, ...
                            'LineStyle','-',...
                            'Color', colors(d,:), ...
                            'DisplayName', string(plotData(d).date));
                    end
                end
                title(ax, sprintf('%s lead', leads{l}), 'Interpreter', 'none');
                ylabel(ax, 'Power (dB/Hz)');
                grid(ax, 'on');
                % set(ax, 'XScale', 'log');
                xlim('tight')
                axis square
            end
            hold off;
            if l == nLeads
                xlabel(ax, 'Frequency (Hz)');
            else
                ax.XTickLabel = [];
            end

            colormap(f, flipud(colors));
            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.Label.String = 'Recording date';
            nTicks = 11;
            dateIdx = round(linspace(1, N, nTicks)); % Corresponding dates
            dateLbls = [plotData.date];
            cb.TickLabels = fliplr(cellstr(dateLbls(dateIdx)));
            title(tl, sprintf('PSD evolution over time (%s, %s)', refType, statType), ...
                'Interpreter', 'none');

            sgtitle(sprintf('%s PSD Evolution \n %s', statType, refType));
            % legend('Location', 'eastoutside');
            fontsize(16,'points')

        end

        function plotSpikeDistribution(obj, refType)
            % plotSpikeDistribution - Bar chart of spike counts per channel (filtered/sorted)
            % Inputs:
            %   refType - string ('Native', 'CAR', etc.)

            if ~isfield(obj.subjectSummary, 'SpikeCountPerChannel')
                error('subjectSummary does not contain SpikeCountPerChannel field.');
            end

            allRefs = fieldnames(obj.subjectSummary.SpikeCountPerChannel);
            if ~ismember(refType, allRefs)
                error('Ref type "%s" not found in subjectSummary.SpikeCountPerChannel.', refType);
            end

            % Extract and validate
            spikeData = obj.subjectSummary.SpikeCountPerChannel.(refType);

            % --- Normalize input ---
            if isnumeric(spikeData)
                chLabels = arrayfun(@(ch) sprintf('Ch %d', ch), 1:numel(spikeData), 'UniformOutput', false);
                spikeCounts = spikeData(:);
            elseif istable(spikeData)
                chLabels = spikeData.Channel;
                spikeCounts = spikeData.SpikeCount;
            else
                error('Unsupported format for spike count data.');
            end

            % --- Exclude BAD-Ch entries ---
            isBad = startsWith(chLabels, 'BAD-Ch', 'IgnoreCase', true);
            chLabels = chLabels(~isBad);
            spikeCounts = spikeCounts(~isBad);

            % --- Parse channel numbers ---
            chNumbers = zeros(size(chLabels));
            for i = 1:numel(chLabels)
                label = chLabels{i};
                if startsWith(label, 'Ch ')
                    chNumbers(i) = sscanf(label, 'Ch %d');
                elseif contains(label, '-')  % e.g., '5-12'
                    tokens = regexp(label, '^(\d+)-\d+', 'tokens');
                    if ~isempty(tokens)
                        chNumbers(i) = str2double(tokens{1}{1});
                    else
                        chNumbers(i) = NaN;
                    end
                else
                    chNumbers(i) = NaN;  % fallback
                end
            end

            % --- Sort logic ---
            if numel(chLabels) > 32
                [spikeCounts, sortIdx] = sort(spikeCounts, 'descend');
                chLabels = chLabels(sortIdx);
            else
                [chNumbers, sortIdx] = sort(chNumbers);
                chLabels = chLabels(sortIdx);
                spikeCounts = spikeCounts(sortIdx);
            end

            % --- Plot ---
            figure('Name', 'Spike Channel Distribution', 'Position',[249 127 1259 816]);
            bar(spikeCounts, 'FaceColor', [0.2 0.5 0.8]);
            set(gca, 'XTick', 1:numel(chLabels), ...
                'XTickLabel', chLabels, ...
                'XTickLabelRotation', 45);
            xlabel('Channel');
            ylabel('Spike Count');
            title(['Spike Distribution - ', refType]);
            if numel(chLabels) > 32
                xlim([0 35.6])
            end
            grid on;
            fontsize(16,'points')
        end
        
        function plotNormalizedSpikeCounts(obj, refType)
            % plotNormalizedSpikeCounts - Per-channel boxplots of normalized spike rates (spikes/min)
            %
            % Inputs:
            %   refType - rereferencing type (default = 'Native')

            arguments
                obj
                refType (1,:) char {mustBeMember(refType, ...
                    {'Native','CAR','leadCAR','leadSpecific','BP'})} = 'Native'
            end

            restRecs = obj.recordings(~arrayfun(@(r) isfield(r, 'StimInfo') && ~isempty(r.StimInfo), obj.recordings));
            restRates = struct();

            % --- 1. baseline spike counts ---
            for i = 1:numel(restRecs)
                try
                    if isfield(restRecs(i), refType) && ...
                            isfield(restRecs(i).(refType), 'SpikeDetect') && ...
                            isfield(restRecs(i).(refType).SpikeDetect, 'NumberOfDenoisedSpikes')

                        spikeCounts = restRecs(i).(refType).SpikeDetect.NumberOfDenoisedSpikes;
                        durMin = minutes(duration(string(restRecs(i).Duration)));
                        restRates(i).date = dateshift( datetime({restRecs(i).date}, ...
                            'InputFormat', 'd-MMM-y HH:mm:ss'), 'start', 'day');
                        restRates(i).rateTotal = spikeCounts / (durMin * numel(restRecs(i).(refType).ChannelNames)); % Normalized rate per minute per channel
                        restRates(i).chNames = restRecs(i).(refType).ChannelNames;
                    end
                catch
                    continue;
                end
            end

            % ######### NEEDS TO BE CHANGED ############
            % --- 2. Collect channel-wise spike counts from long-term stim ---
            % stimRecs = obj.longtermStimData;
            % stimRates = struct();
            % 
            % for i = 1:numel(stimRecs)
            %     try
            %         if isfield(stimRecs(i), refType) && ...
            %                 isfield(stimRecs(i).(refType), 'SpikeDetect') && ...
            %                 isfield(stimRecs(i).(refType).SpikeDetect, 'stat') && ...
            %                 isfield(stimRecs(i), 'StimInfo')
            % 
            %             spikeCounts = stimRecs(i).(refType).SpikeDetect.stat.num_denoised;
            % 
            %             % Stim-free (rest) duration in minutes
            %             sc = stimRecs(i).StimInfo.StimSequence.Sequence;
            %             restDur = seconds(sum(stimRecs(i).StimInfo.StimSequence.SequenceDuration(sc == 0)) / obj.samplingRate);
            %             stimRates(i).date = dateshift(datetime({stimRecs(i).date}, ...
            %                 'InputFormat', 'd-MMM-y HH:mm:ss'), 'start', 'day');
            %             stimRates(i).rateTotal = sum(spikeCounts) / minutes(restDur);
            %             stimRates(i).ratePerCh = spikeCounts / minutes(restDur);  
            %             stimRates(i).chNames = stimRecs(i).(refType).ChannelNames;
            %         end
            %     catch
            %         continue;
            %     end
            % end
            % % % --- Calculate statistics ----
            % [h,p] = ttest2([restRates.rateTotal], [stimRates.rateTotal]);
            % --- 3.1 Plot Total spike cound ditstriubution ---
            % figure('Name', 'Normalized Spike Counts');
            % % ---- Box plot -------
            % % boxplot([restRates(:); stimRates(:)], ...
            % %     [repmat({'Baseline'}, numel(restRates), 1); ...
            % %     repmat({'Stim'}, numel(stimRates), 1)],'Notch','on');
            % % ---- Box Chart -----
            % boxchart([repmat(1, numel(restRates), 1); ...
            %     repmat(2, numel(stimRates), 1)],...
            %     [[restRates.rateTotal]'; [stimRates.rateTotal]']);
            % xticks([1,2])
            % xticklabels({'Baseline', 'Long-term Stimulation'})
            % % --------
            % ylabel('Spikes per Minute');
            % title(['Normalized Spike Count Comparison - ', refType]);
            % subtitle(['p=',num2str(p)])
            % grid on;
            % fontsize(16,'points')

            % -----3.2 Plot spike distribution time development ----
            restTable = timetable([restRates.date]',[restRates.rateTotal]');
            uniqueRowsTT = unique(sortrows(restTable));
            uniqueTimes = unique(uniqueRowsTT.Time);
            spikesOverTime = retime(uniqueRowsTT,uniqueTimes,'mean');
            
            % Smooth the data
            % dataTrend = smoothdata(spikesOverTime.Var1,'movmean',3);
            % dataTrend = smoothdata(spikesOverTime.Var1,'sgolay',5);
            [~,peakIdx] = findpeaks(spikesOverTime.Var1,'MinPeakProminence',4);  % locate prominent peak
            weights = 0.1 * ones(size(spikesOverTime.Var1));             % Assign weights: high at the peak, low elsewhere
            weights(peakIdx) = 10;  % prioritize peak preservation
            % Fit a smoothing spline with weights
            x = (1:length(uniqueTimes))';  % or just: (1:length(...))(:)
            y = spikesOverTime.Var1(:);
            dataTrend = fit(x, y, 'smoothingspline', 'Weights', weights, 'SmoothingParam', 0.5);
            
            figure('Name', 'Normalized Spike Time Development','Position',[249 127 1259 816]);
            hold on;
            % scatter(spikesOverTime.Time,spikesOverTime.Var1,'o','black'); hold on % Plot original time points
            plot(spikesOverTime.Time, spikesOverTime.Var1, 'x','color', 'black', 'LineWidth', 2)
            plot(spikesOverTime.Time, dataTrend(x),'Color', [0.9294    0.6941    0.1255],'LineWidth',2)
            title('Normalized Spike Count Time Development')
            ylabel('Spikes/(minute x channel)'); xlabel('Date-time')
            legend({'Spike count','Data Trend'})
            xlim('tight')
            grid on;
            fontsize(16,'points')
        end

        function plotSpikeExamples(obj, refType)
            % plotSpikeExamples Plot random spike waveform examples.
            %
            % Syntax
            %   obj.plotSpikeExamples()
            %   obj.plotSpikeExamples(refType)
            %
            % Description
            %   Loads spike examples from obj.spikeEventsFile and plots up to 64
            %   randomly selected spike waveforms for the requested rereferencing type.

            arguments
                obj
                refType (1,:) char {mustBeMember(refType, ...
                    {'Native','CAR','leadCAR','leadSpecific','BP'})} = 'BP'
            end

            % Load spike example events
            if ~isfile(obj.spikeEventsFile)
                error('Spike events file not found: %s', obj.spikeEventsFile);
            end

            tmp = load(obj.spikeEventsFile);
            spikeExamples = tmp.spikeExamples;

            if isempty(spikeExamples)
                error('No spike examples found in %s.', obj.spikeEventsFile);
            end

            if ~isstruct(spikeExamples)
                error('Loaded spikeExamples variable must be a struct.');
            end

            % Validate recordings with requested reref type
            hasRef = arrayfun(@(x) isfield(x, refType), spikeExamples);

            if ~any(hasRef)
                error('No spike examples found for refType "%s".', refType);
            end

            validRecIdx = find(hasRef);
            nRecs = numel(validRecIdx);

            if nRecs < 10
                error('Not enough recordings with refType "%s" to plot examples. Found only %d recordings.', ...
                    refType, nRecs);
            end

            % Create figure
            figure('Name', sprintf('Spike examples - %s', refType), ...
                'Position', [100 100 1400 1000]);

            tl = tiledlayout(8, 8, ...
                'TileSpacing', 'compact', ...
                'Padding', 'compact');

            % Random plotting loop
            nToPlot = 64;
            pltCount = 1;
            iter = 1;
            maxIter = 3 * nToPlot;

            while pltCount <= nToPlot && iter <= maxIter
                iter = iter + 1;
                % Random recording
                rIdx = validRecIdx(randi(nRecs));
                spikeRec = spikeExamples(rIdx).(refType);

                % Validate spikeRec fields
                requiredFields = {'spikes', 'channelLabels', 'timeStamps'};

                if ~all(isfield(spikeRec, requiredFields))
                    continue
                end

                if isempty(spikeRec.spikes) || isempty(spikeRec.channelLabels) || isempty(spikeRec.timeStamps)
                    continue
                end

                if size(spikeRec.spikes, 1) < 1
                    continue
                end

                % Random spike index
                sIdx = randi(size(spikeRec.spikes, 1));
                spike = spikeRec.spikes(sIdx, :);

                % Time axis centered at zero, in samples
                nSamp = numel(spike);
                t_ax = -nSamp/4 :0.5: (nSamp/4 -0.5);

                % Recording name
                if isfield(spikeExamples, 'name') && ~isempty(spikeExamples(rIdx).name)
                    recName = erase(spikeExamples(rIdx).name, '.dat');
                else
                    recName = sprintf('Recording %d', rIdx);
                end

                % Channel label
                if numel(spikeRec.channelLabels) >= sIdx
                    chLabel = spikeRec.channelLabels{sIdx};
                else
                    chLabel = 'Unknown channel';
                end

                % Plot pike
                ax = nexttile(tl, pltCount);
                plot(ax, t_ax, spike, ...
                    'Color', [0.1 0.1 0.1], ...
                    'LineWidth', 1); 
                hold on
                % Plot gray line marking 100 µV
                scaleIdx = round(0.1 * numel(spike));   % position near the left side
                scaleIdx = max(1, min(scaleIdx, numel(spike)));
                xScale = t_ax(scaleIdx);
                yCenter = spike(scaleIdx);
                plot(ax, [xScale xScale], ...
                    yCenter + [-50 50], ...
                    'LineWidth', 2, 'Color', [0.651 0.651 0.651]);

                xlim(ax, [t_ax(1), t_ax(end)]);
                title(ax, sprintf('%s | %s', recName, chLabel), ...
                    'Interpreter', 'none', ...
                    'FontSize', 7);
                pltCount = pltCount + 1;
            end
            % Fallback warning
            if pltCount <= nToPlot
                warning('Only %d/%d spike examples could be plotted for refType "%s".', ...
                    pltCount-1, nToPlot, refType);
            end
            fontsize(gcf, 10, 'points');

            sgtitle(sprintf('Random spike examples (%s)', refType), ...
                'Interpreter', 'none', ...
                'FontSize', 16, ...
                'FontWeight', 'bold');
        end

        function plotImpedanceDevelopment(obj, impTable)
            % If implant table was not inserted then call  summarize
            % impedance
            if nargin <2
                impTable = obj.summarizeImpedance();
            end
            chnls = unique([values(obj.leadAssignment).indices]);
            leads = keys(obj.leadAssignment);
            % Color coded Impedance Development
            f1 = figure('Name','Color Coded Impedance development', 'Position', [100 100 1600 900]);
            clrlim = [50 10e3];
            isc = imagesc(impTable.Date,chnls,impTable{:,[chnls]+2}',clrlim); xlabel('Recording date [-]');
            xticks(impTable.Date(1) : calmonths(2) : impTable.Date(end));
            xtickformat('MMM yyyy')
            colormap turbo; c = colorbar; c.Label.String = 'Impedance [Ω]';
            % c.Location='northoutside';
            ylabel('Channel [-]'); yticks(1:32);
            title('Impedance development')
            fontsize(16,'points');
            axis tight

            % Average Lead imp
            f2 = figure('Name','Impedance across the time', 'Position', [150 150 1000 800]);
            clrs = lines(numel(leads));
            hold on;
            for l = 1:numel(leads)
                idxs = obj.leadAssignment(leads(l)).indices;
                % leadAvg = mean(impTable{:,idxs+2}./1000,2); % /.1000 convert to kΩ and Skip first two columns
                % leadSTD = std(impTable{:,idxs+2}./1000,0,2); 
                % lowBound = max(leadAvg-leadSTD,0.001);  % Clip negative values to zero
                % upperBound = flipud(leadAvg+leadSTD);

                leadAvg = median(impTable{:,idxs+2}./1000,2);
                leadVar =  prctile(impTable{:,idxs+2}./1000,[20 80],2);
                lowBound = leadVar(:,1);
                upperBound =flipud(leadVar(:,2));

                % Patch between mean ± std
                fill([impTable.Date; flipud(impTable.Date)], [lowBound; upperBound], ...
                    clrs(l,:), ...
                    'FaceAlpha', 0.2, ...
                    'EdgeColor', 'none', ...
                    'HandleVisibility', 'off');  % Don't show shaded area in legend

                plot(impTable.Date, leadAvg, ...
                    'LineWidth',2.5, ...
                    'Color',clrs(l,:), ...
                    'DisplayName',leads(l));
            end
            grid on;
            xlim('tight')
            xlabel('Date');
            ylabel('Impedance (k\Omega)');
            legend('Location', 'best');
            title('Average Lead Impedance Over Time');
            fontsize(16,'points');
        end

        function [plotData] = plotBadChannelOccurence(obj, impTable)
            % If implant table was not inserted then call  summarize
            % impedance
            if nargin <2
                impTable = obj.summarizeImpedance();
            end
            
            % ~~~~~ Plot overall summary of bad channel occurence ~~~~~
            % Extract numeric part of channel names
            overallBadChOccurenece = obj.subjectSummary.BadChannelOccurrence;
            channelNums = sscanf(sprintf('%s\n', overallBadChOccurenece.Channel{:}), 'Ch %d\n');
            % Sort by the numeric values
            [sortedNums, sortIdx] = sort(channelNums);
            % Apply the sort to the table
            T_sorted = overallBadChOccurenece(sortIdx, :);

            % Plot as bar chart
            figure('Name','Bad Channel Occurenece Count', 'Position', [100 100 1500 800]);
            bar(T_sorted.Count)
            xticks(1:height(T_sorted))
            xticklabels(T_sorted.Channel)
            xtickangle(45)
            ylabel('Count')
            title('Bad Channel Detection Count')
            fontsize(16,"points")

            % ~~~~~ Plot month by month bad channel occurence ~~~~~
            monthlyBadChannels = obj.subjectSummary.MonthlyBadChannelOccurrence;
            monthFields = fieldnames(monthlyBadChannels);
            plotData = table('Size', [numel(monthFields) 3], ...
                 'VariableTypes', {'datetime','double','cellstr'}, ...
                 'VariableNames', {'Month','TotalGoodCh','BadChLabels'});
            for m = 1:numel(monthFields)
                month = datetime(monthFields{m},'InputFormat','MMM_yyyy','Format','MMM-yyyy');
                if isempty(monthlyBadChannels.(monthFields{m}).BadChannelOccurence)
                    totalGoodCh = 32;
                    badChLabels = {''};
                else
                    % Channel is classified as bad if was identified bad at
                    % least in 10% of recordings;
                    badChIdx = monthlyBadChannels.(monthFields{m}).BadChannelOccurence.Percentage >= 10;
                    totalGoodCh = 32 - sum(badChIdx);
                    badChnls = monthlyBadChannels.(monthFields{m}).BadChannelOccurence.Channel(badChIdx);
                    badChLabels = {strjoin("Ch" + badChnls, '; ')};
                end
                plotData.Month(m) = month;
                plotData.TotalGoodCh(m) = totalGoodCh;
                plotData.BadChLabels(m) = badChLabels;
            end
            figure('Name','Channel Count', 'Position', [100 100 1500 800]);
            stairs(plotData.Month,plotData.TotalGoodCh,'LineWidth',2,'Marker','*')
            yline(32,'LineStyle','--')
            ylim([0 33])
            xticks(plotData.Month(1) : calmonths(2) : plotData.Month(end));
            xlim('tight')
            % xtickangle(45)
            grid on
            ylabel('Count')
            title('Channel Count')
            fontsize(16,"points")
        end

%% ~~~~~~~~~~~~~~~~~~ Helper Functions ~~~~~~~~~~~~~~~~~~~~~~~~
        function montage = getMontage(obj)
            % Returns channel montage for each  lead based on the subject construct
            switch obj.rerefSpecs.construct
                case "surface"
                    g1 = [(1:6)', (7:12)'];
                    g2 = [(13:17)', (18:22)'];
                    g3 = [(23:27)', (28:32)'];
                    montage = {g1; g2; g3};

                case "depth"
                    montage = {1:8;
                        9:16;
                        17:24;
                        25:32};

                case "combined"
                    montage = {1:8;
                        9:16;
                        17:24;
                        25:32};
            end
        end

%% ~~~~~~~~~~~~~~~~~~~~ Dependencies ~~~~~~~~~~~~~~~~~~~~~~~~~~
        function [obj, report] = checkDependencies(obj, group, mode)
            % checkDependencies Check whether required dependencies are available.
            % Usage:
            %   [obj, report] = obj.checkDependencies()
            %   [obj, report] = obj.checkDependencies('core')
            %   [obj, report] = obj.checkDependencies('Spectral', 'warning')
            %
            % Inputs:
            %   group - 'core' (default), 'Spectral', 'SpikeDetection',
            %           'LongtermStim', or 'all'
            %   mode  - 'warning' (default), 'error', or 'silent'
            %
            % Outputs:
            %   obj    - updated object with DependencyReport.(group)
            %   report - struct containing dependency check results

            if nargin < 2 || isempty(group)
                group = 'core';
            end

            if nargin < 3 || isempty(mode)
                mode = 'warning';
            end

            validModes = {'warning','error','silent'};

            if ~ismember(lower(mode), validModes)
                error('Invalid mode "%s". Must be one of: %s', mode, strjoin(validModes, ', '));
            end

            deps = obj.getDependencyGroup(group);

            items = repmat(struct( ...
                'name', "", ...
                'type', "", ...
                'required', false, ...
                'notes', "", ...
                'found', false, ...
                'resolvedPath', ""), numel(deps), 1);

            missingRequired = strings(0,1);
            missingOptional = strings(0,1);

            for k = 1:numel(deps)
                dep = deps(k);
                [found, resolvedPath] = obj.resolveDependency(dep);
                items(k).name = string(dep.name);
                items(k).type = string(dep.type);
                items(k).required = dep.required;
                items(k).notes = string(dep.notes);
                items(k).found = found;
                items(k).resolvedPath = string(resolvedPath);

                if ~found
                    if dep.required
                        missingRequired(end+1,1) = string(dep.name); %#ok<AGROW>
                    else
                        missingOptional(end+1,1) = string(dep.name); %#ok<AGROW>
                    end
                end
            end

            report = struct();
            report.group = string(group);
            report.checkedAt = datetime('now');
            report.items = items;
            report.ok = isempty(missingRequired);
            report.missingRequired = missingRequired;
            report.missingOptional = missingOptional;
            groupField = matlab.lang.makeValidName(char(string(group)));
            obj.DependencyReport.(groupField) = report;

            % Message handling
            if ~report.ok
                msg = sprintf('Missing required dependencies for group "%s": %s', ...
                    group, strjoin(cellstr(missingRequired), ', '));

                switch lower(mode)
                    case 'warning'
                        warning(msg);
                    case 'error'
                        error(msg);
                    case 'silent'
                        % do nothing
                end

            elseif ~isempty(missingOptional) && strcmpi(mode, 'warning')
                warning('Missing optional dependencies for group "%s": %s', ...
                    group, strjoin(cellstr(missingOptional), ', '));
            end
        end

        function printDependencyReport(obj, group)
            % printDependencyReport Print dependency check report to command window.
            % Usage:
            %   obj.printDependencyReport()
            %   obj.printDependencyReport('core')

            if nargin < 2 || isempty(group)
                group = 'core';
            end

            groupField = matlab.lang.makeValidName(char(string(group)));
            if ~isfield(obj.DependencyReport, groupField)
                fprintf('No dependency report stored for group "%s".\n', group);
                return;
            end

            rpt = obj.DependencyReport.(groupField);
            fprintf('\nDependency report for group: %s\n', rpt.group);
            fprintf('Checked at: %s\n', char(rpt.checkedAt));
            fprintf('Overall status: %s\n', string(rpt.ok));

            for k = 1:numel(rpt.items)
                status = "MISSING";
                if rpt.items(k).found
                    status = "FOUND";
                end
                fprintf('  [%s] %-28s (%s)', status, rpt.items(k).name, rpt.items(k).type);

                if strlength(rpt.items(k).resolvedPath) > 0
                    fprintf(' -> %s', rpt.items(k).resolvedPath);
                end
                fprintf('\n');

            end
            fprintf('\n');

            %%
        end
    end

    methods (Static, Access = private)
        function deps = getDependencyGroup(group)
            % getDependencyGroup Return dependency registry for a given group.
            core = BIC_data.getCoreDependencies();
            optional = BIC_data.getOptionalDependencies();

            switch lower(char(string(group)))
                case 'core'
                    deps = core;
                case 'spectral'
                    deps = optional.Spectral;
                case 'spikedetection'
                    deps = optional.SpikeDetection;
                case 'longtermstim'
                    deps = optional.LongtermStim;
                case 'all'
                    deps = [core; optional.Spectral; optional.SpikeDetection; optional.LongtermStim];
                otherwise
                    error('Unknown dependency group "%s".', group);
            end
        end

        function [found, resolvedPath] = resolveDependency(dep)
            % resolveDependency Resolve dependency availability.
            found = false;
            resolvedPath = "";
            switch lower(dep.type)
                case {'function','class','mex'}
                    resolvedPath = which(dep.name);
                    found = ~isempty(resolvedPath);
                case 'folder'
                    found = isfolder(dep.name);
                    if found
                        resolvedPath = dep.name;
                    end
                case 'toolbox'
                    v = ver;
                    found = any(strcmpi({v.Name}, dep.name));
                    if found
                        resolvedPath = dep.name;
                    end
                otherwise
                    error('Unknown dependency type "%s" for "%s".', dep.type, dep.name);
            end
        end

        function deps = getCoreDependencies()
            %% getCoreDependencies Core dependencies needed by the class.
            mk = @(name,type,required,notes) struct( ...
                'name', name, ...
                'type', type, ...
                'required', required, ...
                'notes', notes);

            deps = [ ...
                mk('load_bcidat',              'function', true,  'BCI2000 .dat loading'); ...
                mk('extractImpedanceValues',    'function', true,  'Extracts Impedace values from .txt or .log files'); ...
                mk('initializeDataStruct',     'function', true,  'Initializes data structure for recordings'); ...
                mk('custom_reref',             'function', true,  'Re-referencing helper'); ...
                mk('interpolate_lost_samples', 'function', true,  'Lost sample interpolation'); ...
                mk('calculatePacketLoss',      'function', true,  'Calculate average packet-loss'); ...
                mk('parseBCI2000_parameters_values', 'function', true, 'BCI2000 parameter extraction helper'); ...
                mk('runlength',                'function', true,  'Parses the vector x into runs of constant value. Available from https://www.mathworks.com/matlabcentral/fileexchange/241-runlength-m');...
                mk('identifyBadChannels',      'function', true,  'Identifies bad channels in the recording'); ...
            ];

        end

        function optional = getOptionalDependencies()
            %% getOptionalDependencies Optional dependencies grouped by feature.
            mk = @(name,type,required,notes) struct( ...
                'name', name, ...
                'type', type, ...
                'required', required, ...
                'notes', notes);

            optional = struct();

            optional.FileProcessing = [ ...
                mk('progressbar',              'function', false, 'Progressbar in Command Window; Available from https://www.mathworks.com/matlabcentral/fileexchange/167996-matlab-progress-bar'); ...
                mk('struct2single',            'function',  true,  'Converts elements of struct to single');...
            ];

            optional.Spectral = [ ...
                mk('extractSpectralFeatures', 'function', true,  'Spectral feature extraction'); ...
                mk('channelCrossPSD',         'function', true,  'Cross spectral density'); ...
                mk('cpsdToCoherence',         'function', true,  'Coherence from CPSD') ...
            ];

            optional.SpikeDetection = [ ...
                mk('Spike_detection_function', 'function', true, 'Spike detection'); ...
                mk('pSpike_elimination_RFOMP', 'function', true, 'Spike denoising / RFOMP') ...
            ];

            optional.LongtermStim = [ ...
                mk('parsave', 'function', true, 'Parallel save helper') ...
            ];
        end

    end
end
