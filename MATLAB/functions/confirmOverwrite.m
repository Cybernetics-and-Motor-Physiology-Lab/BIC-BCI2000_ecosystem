function proceed = confirmOverwrite(filepath)
    % ----- Helper for overwriting confirmation -----
    if isfile(filepath)
        choice = questdlg(sprintf('File "%s" already exists. Overwrite?', ...
            filepath), 'Confirm Overwrite', 'Yes', 'No', 'No');
        proceed = strcmp(choice, 'Yes');
    else
        proceed = true;
    end
end