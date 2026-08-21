function file_path = get_paths(dir_path,extension, method)
% --------- Function to return path to the data based on conditions ------
% -------------------------------------------------------------------------
% --Inputs:
%   Recquired: -> dir_path: -char, path to the directory which to list
%              -> extension: -char, regular expression to filter the list
%              of the data based on the extension
%   Optional: -> method: -char ('last', 'first', 'biggest')
%                        -returns path to file fullfiling the condidtion
%
% --Output:
%   -> file_path: - cell/char
%       - if no method is specified returns all the datapaths matching the stated expression in cell
%       - if method is specified, returns the dtapath to file matching the condition
% Frederik Lampert 2024

fl_list = struct2table(dir([dir_path filesep '*' extension]));
switch nargin
    case 2
        pths = fl_list{:,'folder'}; names = fl_list{:,'name'};
        file_path = strcat(pths, names);
    case 3
        switch method
            case 'last'
                fl_list = sortrows(fl_list, 'datenum', 'desc');
                file_path = char(strcat(fl_list{1,'folder'}, filesep, fl_list{1,'name'}));
            case 'first'
                fl_list = sortrows(fl_list, 'datenum', 'asc');
                file_path = char(strcat(fl_list{1,'folder'}, filesep, fl_list{1,'name'}));
            case 'biggest'
                fl_list = sortrows(fl_list, 'bytes', 'desc');
                file_path = char(strcat(fl_list{1,'folder'}, filesep, fl_list{1,'name'}));
            otherwise
                error(['Invalid input for method "%s", select ' ...
                    'from one of these: -"first" -"last" -"biggest"'], method)
        end
end