function dataFiles = findDataFiles(path, existing_files)
    % ############################################
    % Function to search and return all the data files in the selected
    % directory and subdirectories
    % Inputs: path - path to directory containing the data
    %         existing_files - structure with already analyzed/discarded
    %         data
    % Outputs: dataFiles - struct containing the name, path and date of the
    %                      (new) files
    % ############################################
    if nargin < 2
        disp(['Searching for Data files in ' path]) % Print sentence into the console
        % Get a list of all files matching the pattern
        dataFiles = dir(fullfile(path, '**','*.dat'));
    else
        disp(['Searching for new Data in ' path]) % Print sentence into the console
        % Get a list of all files matching the pattern
        dataFiles = dir(fullfile(path, '**','*.dat'));
        % Compare and Remove existing files
        [new, new_idx] = setdiff({dataFiles.name},{existing_files.name});
        dataFiles = dataFiles(new_idx);
    end
    
    emptyIdx = [dataFiles.bytes] == 0;   % Search for empty files
    dataFiles(emptyIdx) = [];                           % Remove the empty files
    dataFiles = rmfield(dataFiles,'isdir');             % Remove isdir field
    
end