function expandedStruct = initializeDataStruct(existingStruct, rerefTypes)
% initializeDataStruct - Expands an existing dir() structure with new fields

% Validate input
requiredFields = {'name', 'folder', 'date', 'bytes', 'datenum'};
if ~all(ismember(requiredFields, fieldnames(existingStruct)))
    error('Input structure is missing required dir() fields.');
end

% Preallocate output
nFiles = numel(existingStruct);
expandedStruct = existingStruct;  % start with the existing fields

% Define empty initial values for new fields
for i = 1:nFiles
    expandedStruct(i).SamplingRate = [];
    expandedStruct(i).NativeRef = [];
    expandedStruct(i).Duration = "";
    expandedStruct(i).PacketLoss = [];
    expandedStruct(i).AmplFactor = '';
    expandedStruct(i).Ground = [];
    expandedStruct(i).Impedance = {};
    expandedStruct(i).StimInfo = [];
    expandedStruct(i).BadChannels = [];
    expandedStruct(i).ProcessingErrorMsg = {};

    % Initialize rereferencing fields
    for ridx = 1:length(rerefTypes)
        ref = rerefTypes{ridx};
        expandedStruct(i).(ref) = struct();
    end
end
end