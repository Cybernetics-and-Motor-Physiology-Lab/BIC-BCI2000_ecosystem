function uniqueStruct = removeDuplicateRecordings(recordingsStruct)
    % removeDuplicateRecordings - Removes duplicate entries from a recordings structure
    %
    % Duplicates are identified by matching both 'name' and 'folder' fields

    if isempty(recordingsStruct)
        uniqueStruct = recordingsStruct;
        return;
    end

    % Create composite keys from name + folder
    keys = strcat({recordingsStruct.folder}, '|', {recordingsStruct.name});
    
    % Use unique to find first occurrence
    [~, uniqueIdx] = unique(keys, 'stable');
    
    % Keep only unique entries
    uniqueStruct = recordingsStruct(uniqueIdx);

    % Report
    nRemoved = numel(recordingsStruct) - numel(uniqueStruct);
    if nRemoved > 0
        fprintf('Removed %d duplicate recordings.\n', nRemoved);
    else
        disp('No duplicates found.');
    end
end