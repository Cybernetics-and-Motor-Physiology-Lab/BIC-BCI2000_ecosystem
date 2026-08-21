function mustBeFile(file)
% Helper function for BIC_data 
    if ~isfile(file)
        error('File "%s" does not exist.', file);
    end
end