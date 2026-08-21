function ImpFiles = findImpFiles(path, existing_files)
    % ############################################
    % Function to search and return all the impedance files in the selected
    % directory and subdirectories
    % Inputs: path - path to directory containing the Imp data
    %         existing_files - structure with already analyzed/discarded
    %         Imedance data
    % Outputs: ImpFiles - struct containing the name, path and date of the
    %                     files
    % ############################################
    disp(['Searching for Impedance files in ' path]) % Print sentence into the console

    % Specify all the possible patterns for impedance files
    fls_pattern = {
        '*_impedances*.txt';
        '*_info*.txt';
        'Grd_*.log'
        };

    % Get a list of all files matching the first pattern
    ImpFiles = dir(fullfile(path, '**',fls_pattern{1}));
    for pat = 2:length(fls_pattern)
        ImpFiles = [ImpFiles; dir(fullfile(path, '**',fls_pattern{pat}))];
    end

    if nargin > 1
        [new, new_idx] = setdiff({ImpFiles.name},{existing_files.name});
        ImpFiles = ImpFiles(new_idx);
    end


    emptyIdx = [ImpFiles.bytes] == 0;   % Search for empty files
    ImpFiles(emptyIdx) = [];                           % Remove the empty files
    ImpFiles = rmfield(ImpFiles,'isdir');              % Remove isdir field
end