function impedanceValues = extractImpedanceValues(filename)

    format = strsplit(filename, '.'); 
    format = format{end};
    
    % Initialize a cell to store impedance values
    impedanceValues = {};
    
    % Open the file for reading
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Failed to open the file: %s', filename);
    end
    
    switch format
        case 'txt'
            % Read the file line by line
            tline = fgetl(fileID);
            while ischar(tline)
                % Match the lines that contain channel and impedance values
                tokens = regexp(tline, '^Ch\s*(\d+)\s+([\d\.inf]+)$', 'tokens');
                if ~isempty(tokens)
                    channel = str2double(tokens{1}{1});
                    value = str2double(tokens{1}{2});
                    % Store the impedance value in the structure
                    impedanceValues = [impedanceValues; {channel, value}];
                end
                % Read the next line
                tline = fgetl(fileID);
            end

        case 'log'
            % Read the file line by line
            tline = fgetl(fileID);
            while ischar(tline)
                % Match the lines that contain channel and impedance values
                tokens = regexp(tline, 'ch\s*(\d+)\s*=\s*([\d+\.inf]+)', 'tokens');
                if ~isempty(tokens)
                    channel = str2double(tokens{1}{1})+1; % adjust the indexing from 0 to 1
                    value = str2double(tokens{1}{2});
                    % Store the impedance value in the structure
                    impedanceValues = [impedanceValues; {channel, value}];
                end
                % Read the next line
                tline = fgetl(fileID);
            end
    end
  
    % Close the file
    fclose(fileID);
end