function [reref_signal, chnls_names] = custom_reref(signal,construct,bad_chnls,method,specification, silent)
%--------------------------------------------------------------------------
% Function: custom_reref
%--------------------------------------------------------------------------
% Description:
%   This function performs re-referencing of multi-channel EEG/ECoG signals
%   based on different referencing strategies. It supports different
%   constructs (surface, depth, combined) and methods (CAR, BP, leadCAR, leadSpecific).
%
% Usage:
%   [reref_signal, chnls_names] = custom_reref(signal, construct, 
%                                              bad_chnls, method, specification)
%
% Input Arguments:
%   - signal (double matrix)       : Input signal matrix (samples x channels)
%   - construct (string)           : Type of construct ('surface', 'depth', 'combined')
%   - bad_chnls (vector, optional) : Indices of bad channels to exclude (default: [])
%   - method (string, optional)    : Re-referencing method ('CAR', 'BP', 
%                                    'leadCAR', 'leadSpecific') (default: 'CAR')
%   - specification (optional)     : Additional parameter depending on method
%                                    (string 'NN/horizontal'/'vertical' for surface 'BP';
%                                     string specifying type ('D'-directional, 'S'-standard) of DBS 
%                                     electrodes separated by dash for depth 'BP' (eg.'D-S-D-S'); 
%                                     or a vector of reference channels (eg. [1,13,23] for 'leadSpecific')
%   - silent (optional)            : True/False statement whether the
%                                   function should print the warnings
%                                   (preset to false)
%
% Output Arguments:
%   - reref_signal (double matrix) : Re-referenced signal
%   - chnls_names (cell array)     : Updated channel names after re-referencing
%
% Methods:
%   - CAR           : Common Average Reference
%   - BP            : Differential Pair (Bi-polar derivation) within leads
%   - leadCAR       : Common Average Reference within leads
%   - leadSpecific  : Referenced to the specified channel within the lead
%
% Dependencies:
%   - "car" function (common average reference)
%   - "surface_BP", "depth_BP" for BP re-referencing
%
% ~~~~~~~~~~~~~~~~~~ Chanel placement on surface grids ~~~~~~~~~~~~~~~~
% ###### Grid1 ######## | ####### Grid2 ####### | ###### Grid3 ########
% ###### 1 # 7 ######## | ##### 13 # 18 ####### | ##### 23 # 28 #######
% ##### 2 ### 8 ####### | #### 14 ### 19 ###### | #### 24 ### 29 ######
% #### 3 ##### 9 ###### | ### 15 ##### 20 ##### | ### 25 ##### 30 #####
% ### 4 ####### 10 #### | ## 16 ####### 21 #### | ## 26 ####### 31 ####
% ## 5 ########## 11 ## | # 17 ########## 22 ## | # 27 ########## 32 ##
% #6 ############# 12 # |
% ~~~~~~~~~~~~~~~~~~ Chanel placement on DBS electrods ~~~~~~~~~~~~~~~~
% ###### Conventional ###### | ##### Directional ####
% ######   Tail     ######## | #######   Tail  #######
% ######    (8)     ######## | #######   (8)   #######
% ######    (7)     ######## | ##### (5)|(6)|(7) ######
% ######    (6)     ######## | ##### (2)|(3)|(4) ######
% ######    (5)     ######## | #######   (1)   #######
% ######    (4)     ######## | #######  Tip   #######
% ######    (3)     ######## | Differential pairs are between
% ######    (2)     ######## | the Tip alectrodes and all of the 
% ######    (1)     ######## | first ring electrodes
% ######    Tip     ######## | 
% ############################## 
% ~~~~~~~~~~~~~~~~~~~~~ SURFACE BP ORIENTATION ~~~~~~~~~~~~~~~~~~~~~~~~~
%  -> NN - searches for nearest neighbor channels, if the
%       neighbor is bad, skips the differential pair
%  -> horizontal - horizontal piars, if neighbor is bad, searches for
%       another channels
%  -> vertical - horizontal piars, if neighbor is bad, searches for
%       another channels
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author:
%     Frederik Lampert, Mayo Clinic , 2025.
%--------------------------------------------------------------------------

    arguments
        signal {mustBeVectorOrMatrixDouble}                      % Must be a non-empty double vector or matrix
        construct (1,:) string {mustBeMember(construct, ["surface","depth","combined", "surface_wide"])}
        bad_chnls (1,:) double {mustBeInteger} = []             % Default empty
        method (1,:) string {mustBeMember(method, ["CAR","BP","leadCAR","leadSpecific"])} = "leadCAR"
        specification {mustBeStringOrVectorOfInts} = []         % Default empty
        silent (1,1) logical = false
    end

    %----------------------------------------------------------
    % Step 1: Define Grid Layouts for Surface, Depth, and Combined Constructs
    %----------------------------------------------------------
    switch construct
        case "surface"
            g1 = [(1:6)', (7:12)'];
            g2 = [(13:17)', (18:22)'];
            g3 = [(23:27)', (28:32)'];
            leads = {g1; g2; g3};

        case "surface_wide"
            g1 = [(1:8)', (9:16)'];
            g2 = [(17:24)', (25:32)'];
            leads = {g1; g2};

        case "depth"
            leads = {1:8;
                9:16;
                17:24;
                25:32};

        case "combined"
            leads = {1:8;
                9:16;
                17:24;
                25:32};

    end


    %----------------------------------------------------------
    % Step 2: Mark Bad Channels in Each Grid
    %----------------------------------------------------------
    for l_idx = 1:length(leads)
        bad_idx = ismember(leads{l_idx},bad_chnls);  % find indices of bad channels
        leads{l_idx}(bad_idx) = NaN;
    end

    %----------------------------------------------------------
    % Step 3: Assign Channel Names
    %----------------------------------------------------------
    % Predefine channel names  (might be changed based on the selected method)
    chnls_names =  arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(signal,2), 'UniformOutput', false);
    chnls_names(bad_chnls) = arrayfun(@(n) sprintf('BAD-Ch %d', n), bad_chnls, 'UniformOutput', false);

    % disp([sprintf('Re-referncing signal using %s, ommiting channels: ', method), ...
    %     sprintf('%d ', bad_chnls)]);

    %----------------------------------------------------------
    % Step 4: Apply Re-referencing Method
    %----------------------------------------------------------
    switch method
        case "CAR"
            reref_signal = zeros(size(signal));             % Initialize output matrix
            reref_signal(:,bad_chnls) = signal(:,bad_chnls);% Put back the signal from bad channels to perserve signal dimensions
            g_idx = sort(rmmissing(reshape(cell2mat(leads),1,[])));       % extractct indices of good channels only in correct order
            reref_signal(:,g_idx) = car(signal(:,g_idx));      % Common average rerefernece with exclusion of the bad channels

        case "BP"
            switch construct
                case {"surface", "surface_wide"}
                    % Check if specification was enetered correctly
                    if ~any(strcmp(specification, {'NN','horizontal', 'vertical'}))
                        orientation = questdlg(['Invalid orientation for re-referencing within grids. ' ...
                            'Please select one of the following options'],'Specify orientation', ...
                            'NN','horizontal','vertical','NN');
                    else
                        orientation = specification;
                    end
                    % Perfor re-referncing
                    [reref_signal, chnls_names] = surface_BP(signal,leads,orientation,silent);

                case "depth"
                    % Check if specification was enetered correctly
                    if isempty(regexp(specification, '^[DS]-[DS]-[DS]-[DS]$', 'once'))
                        error(['Invalid DBS electrode type specification!', ...
                            'Please specify the DBS electrode types as a string separated by dashes, where:', newline, ...
                            '"S" - represents a STANDARD electrode', newline, ...
                            '"D" - represents a DIRECTIONAL electrode', newline, ...
                            'For example: "D-D-S-S" specifies Channels 1-16 as channels from directional leads,',...
                            'and Channels 17-32 as from standard leads.']);
                    else
                        leadTypes = specification;
                    end
                    [reref_signal, chnls_names] = depth_BP(signal,leads,leadTypes, silent);

                case 'combined'
                    % Will be added once we will have data from combined construct
                    % Check if specification was enetered correctly
                    if isempty(regexp(specification, '^[GDS]-[GDS]-[GDS]-[GDS]$', 'once'))
                        error(['Invalid DBS electrode type specification!', ...
                            'Please specify the electrode lead types as a string separated by dashes, where:', newline, ...
                            '"G" - represents a surface ECoG GRID', newline, ...
                            '"S" - represents a STANDARD DBS electrode', newline, ...
                            '"D" - represents a DIRECTIONAL DBS electrode', newline, ...
                            'For example: "G-G-S-D" specifies Channels 1-16 as channels from ECoG grids,',...
                            'and Channels 17-32 as from standard DBS leads.']);
                    else
                        leadTypes = strsplit(upper(specification),'-'); % Extract the lead types
                    end

                    reref_signal = []; % Iniialize empty matrix for ouput
                    chnls_names = {};   % Initialize empty cell for channel labels
                    for i = 1:length(leadTypes)
                        temp_sig = [];
                        temp_names = {};
                        switch leadTypes{i}
                            case 'G'
                                N = floor(length(leads{i})/2);
                                g = [leads{i}(1:N)', leads{i}(N+1:end)'];
                                [temp_sig, temp_names] = surface_BP(signal,{g},'NN',silent);
                            case {'S', 'D'}
                                [temp_sig, temp_names] = depth_BP(signal,leads(i),leadTypes{i},silent);                    
                        end
                        reref_signal = [reref_signal temp_sig];
                        chnls_names = [chnls_names temp_names];
                    end
            end

        case "leadCAR"
            reref_signal = zeros(size(signal));             % Initiate output matrix
            reref_signal(:,bad_chnls) = signal(:,bad_chnls);    % Put back the signal from bad channels to perserve signal dimensions
            for l_idx = 1:length(leads)
                lead_chnls = sort(rmmissing(reshape(leads{l_idx},1,[]))); % Get indices of good channels within a lead
                if numel(lead_chnls) < 2
                    if ~silent
                        warning('Skipping lead %d for leadCAR: fewer than 2 good channels.', l_idx);
                    end
                    continue
                end
                reref_signal(:,lead_chnls) = car(signal(:,lead_chnls)); % Common average rerefernece within the leadwith exclusion of the bad channels
            end


        case "leadSpecific"
            % Check if specification was enetered correctly
            if ~isnumeric(specification)
                reference_chnls = inputdlg('Enter the reference channels for each grid', ...
                    'Eneter grid reference', ...
                    [1 45], ...
                    {'1 9 17 25'});
                reference_chnls = str2num(reference_chnls{1});
            else
                reference_chnls = specification;
            end

            if numel(reference_chnls) ~= length(leads)
                if ~silent
                warning('Reference channel was not specified for each lead')
                end
            end

            reref_signal = zeros(size(signal));             % Initiate output matrix
            reref_signal(:,bad_chnls) = signal(:,bad_chnls);    % Put back the signal from bad channels to perserve signal dimensions

            for l_idx = 1:length(leads)
                lead_chnls = sort(rmmissing(reshape(leads{l_idx},1,[])));   % Get indices of good channels within a lead

                if numel(lead_chnls) < 2
                    if ~silent
                        warning('Skipping lead %d for leadSpecific: fewer than 2 good channels.', l_idx);
                    end
                    continue
                end

                if length(intersect(reference_chnls,lead_chnls)) > 1
                    error('More than 1 reference channel was specified for the lead #%d', l_idx)
                elseif isempty(intersect(reference_chnls,lead_chnls))
                    error('No reference channel was specified for the lead #%d', l_idx)
                else
                    ref_ch = intersect(reference_chnls,lead_chnls);
                    lead_chnls(lead_chnls==ref_ch) = []; % Ommit itself
                end
                reref_signal(:, lead_chnls) = signal(:,lead_chnls) - signal(:,ref_ch);
                chnls_names{ref_ch} = sprintf('REF-Ch %d', ref_ch);
            end
    end
    
    % Final safety check: only true fatal failure
    if isempty(reref_signal) || size(reref_signal,2) == 0
        error('custom_reref:NoOutput', ...
        'Unable to create any re-referenced output signal. Bad channels prevent all possible re-referencing.');
    end
end

%% Helper Functions
function mustBeVectorOrMatrixDouble(x)
    if ~(isnumeric(x) && ismatrix(x) && isa(x,'double') && ~isempty(x))
        error("'signal' must be a non-empty double vector or matrix.");
    end
end

function mustBeStringOrVectorOfInts(x)
    % Accepts either a string or a vector of integers

    % If it's char, convert to string for consistent checks
    if ischar(x)
        x = string(x);
    end

    if isstring(x)
        % valid
        return;
    elseif isnumeric(x)
        % must be integer-valued
        if any(mod(x,1) ~= 0)
            error("'specification' must be a string or a vector of integers.");
        end
    else
        error("specification must be either:\n%s\n%s", ...
            "-> string specifying direction (""NN"" ""horizontal"" / ""vertical"") for method leadBP", ...
            "-> vector of integers specifying reference channels within each grid");
    end
end

function [data]=car(data)
    %this function calculates and returns the common avg reference of 2-d matrix "data".
    % Created by Kai J. Miller 12/07
    
    data=double(data);
    
    transflag=0;
    if size(data,1)<size(data,2)
        data=data';
        transflag=1;
    end
    
    num_chans=size(data,2);
    
    % create a CAR spatial filter
    spatfiltmatrix=[];
    spatfiltmatrix=-ones(num_chans);
    for i=1:num_chans
        spatfiltmatrix(i, i)=num_chans-1;
    end
    spatfiltmatrix=spatfiltmatrix/num_chans;
    
    % perform spatial filtering
    if (isempty(spatfiltmatrix) ~= 1)
        % fprintf(1, 'Spatial filtering\n');
        data=data*spatfiltmatrix;
        if (size(data, 2) ~= size(spatfiltmatrix, 1))
            fprintf(1, 'The first dimension in the spatial filter matrix has to equal the second dimension in the data');
        end
    end
    
    if transflag==1, data=data'; end

end

function [grid_BP_ref, BP_labels] = surface_BP(inp_signal,grid_spec,dp_orientation, silent)
    % cortec_surface_dp - Re-references signal data by differentiating channels
    % based on the surface grid design
    % Syntax:
    %   reref_signal = cortec_surface_dp(signal)
    %   reref_signal = cortec_surface_dp(signal, orientation)
    % Orientation:
    %   NN - searches for nearest neighbor channels, if the neighbor
    %       is bad, skips the BP
    %   horizontal - horizontal piars, if neighbor is bad, searches for
    %       another channels
    %   vertical - horizontal piars, if neighbor is bad, searches for
    %       another channels
    % Author:
    %   Frederik Lampert, Mayo Clinic, 2025.
    
    
    grid_BP_ref = [] ;     % Initiate empty variable (in future prealocate matrix instead of changing it size iteratively)
    BP_labels = {};
    
    switch dp_orientation
        case 'NN'
            for i = 1:length(grid_spec)              % iterate over grids
                if sum(isnan(grid_spec{i}),'all') >= numel(grid_spec{i})-1         % If all channels on the grid are bad - skip the grid (or if only one is good)
                    continue
                end
                % Vertical
                for j = 1:size(grid_spec{i},2)       % iterate over channel columns
                    col_chnls = grid_spec{i}(:,j);   % Retrieve column channels from the grid
                    ch_i1 = col_chnls(~isnan([diff(col_chnls); NaN])); % Extract the positions of channels, which are not followed by NaN
                    ch_i2 = ch_i1+1;                 %  Retrieve the second position
                    if length(ch_i1) >= 1  % If there's at least one good pair
                        grid_BP_ref = [grid_BP_ref (inp_signal(:,ch_i1) - inp_signal(:,ch_i2))];    % Substract the channels ( Bipolar re-reference)
                        pairs = num2cell([ch_i1, ch_i2],2);  % Mark down differnital pairs
                        BP_labels = [BP_labels pairs'];         % Save differential pairs as channel names for future reference
                    else
                    end
                end
                % horizontal
                for j = 1:size(grid_spec{i},1)       % iterate over grid rows
                    selected_chnls = grid_spec{i}(j,:); % Retrieve channels that will be used for bipolar re-referncing
                    if sum(isnan(selected_chnls)) >= 1          % If at least one of the channels is bad, skip the channels
                        continue
                    else            % If both channels are good -> bipolar re-refernce
                        grid_BP_ref = [grid_BP_ref diff(inp_signal(:,selected_chnls),[],2)]; % Bipolar re-referencing
                        BP_labels = [BP_labels {selected_chnls}];         % Save differential pairs as channel names for future reference
                    end
                end               
            end

        case 'vertical'
            % Specify channel to be differntiated
            for i = 1:length(grid_spec)              % iterate over grids
                if sum(isnan(grid_spec{i}),'all') >= numel(grid_spec{i})-1         % If all channels on the grid are bad - skip the grid (or if only one is good)
                    continue
                end
                for j = 1:size(grid_spec{i},2)       % iterate over channel columns
                    selected_chnls = grid_spec{i}(:,j); % Retrieve channels that will be used for bipolar re-referncing
                    selected_chnls(isnan(selected_chnls)) = [];     % Drop the bad channels
                    if length(selected_chnls) >= 2
                        grid_BP_ref = [grid_BP_ref diff(inp_signal(:,selected_chnls),[],2)]; % Bipolar re-reference the data
                        pairs = num2cell([selected_chnls(1:end-1) selected_chnls(2:end)],2);  % Mark down differnital pairs
                        BP_labels = [BP_labels pairs'];         % Save differential pairs as channel names for future reference
                    else
                        if ~silent
                            warning(['Skipping Grid %d for %s bipolar re-referencing: ', ...
                                'too many bad channels.'], i, dp_orientation);
                        end
                        continue
                    end
                end
            end
    
        case 'horizontal'
            % Specify channel to be differntiated
            for i = 1:length(grid_spec)              % iterate over grids
                if sum(isnan(grid_spec{i}),'all') >= numel(grid_spec{i})-1         % If all channels on the grid are bad - skip the grid (or if only one is good)
                    continue
                end
                for j = 1:size(grid_spec{i},1)       % iterate over channel columns
                    selected_chnls = grid_spec{i}(j,:); % Retrieve channels that will be used for bipolar re-referncing
                    if sum(isnan(selected_chnls)) == 2          % If both channels are bad throw error
                        if ~silent
                            warning(['Skipping row %d in Grid %d for %s bipolar re-referencing: ', ...
                                'both channels are bad.'], j, i, dp_orientation);
                        end
                        continue
                    elseif sum(isnan(selected_chnls)) == 1     % If one channel is bad look for a subsitute
                        if isnan(selected_chnls(1))            % If the NaN is at 1st columne
                            temp = fillmissing(grid_spec{i}(:,1),'nearest'); % Replace with the closest non Nan Value the same side
                            selected_chnls(1) = temp(j);
                            grid_BP_ref = [grid_BP_ref diff(inp_signal(:,selected_chnls),[],2)]; % Bipolar re-referencing
                            BP_labels = [BP_labels {selected_chnls}];         % Save differential pairs as channel names for future reference
                        else
                            temp = fillmissing(grid_spec{i}(:,2),'nearest'); % Replace with the closest non Nan Value the same side
                            selected_chnls(2) = temp(j);
                            grid_BP_ref = [grid_BP_ref diff(inp_signal(:,selected_chnls),[],2)]; % Bipolar re-referencing
                            BP_labels = [BP_labels {selected_chnls}];         % Save differential pairs as channel names for future reference
                        end
                    else            % If both channels are good -> bipolar re-refernce
                        grid_BP_ref = [grid_BP_ref diff(inp_signal(:,selected_chnls),[],2)]; % Bipolar re-referencing
                        BP_labels = [BP_labels {selected_chnls}];         % Save differential pairs as channel names for future reference
                    end
                end
            end
        otherwise
            error('Invalid surface BP orientation: %s', dp_orientation);
    end
    BP_labels = cellfun(@(v) sprintf('Ch%d-%d', v(1), v(2)), BP_labels, 'UniformOutput', false);
end

function [lead_BP_ref, BP_labels] = depth_BP(inp_signal,ch_numbering,lead_types, silent)
    nch = 8;    % Deafult number of channels in lead, in futur maybe potentialychanged to 16
    lead_types = strsplit(upper(lead_types),'-'); % Extract the lead types
    if numel(lead_types) ~= numel(ch_numbering)
        error('Invalid DBS electrode type specification! Number of lead types in specification must equal to the number of leads')
    end
    lead_BP_ref = []; % Iniialize empty matrix for ouput
    BP_labels = {};   % Initialize empty cell for channel labels
    for i = 1:length(ch_numbering)
         if sum(isnan(ch_numbering{i}),'all') >= numel(ch_numbering{i})-1         
             continue       % If all channels on the grid are bad - skip the grid (or if only one is good)
         end
         selected_chnls = sort(rmmissing(ch_numbering{i})); 
         switch lead_types{i}
             case 'S'
                 % Standard lead bipolar pairs.
                 % Allow pairs that are adjacent or have only one missing contact between them.
                 % Examples allowed: Ch5-Ch6  -> gap = 1 or Ch5-Ch7  -> gap = 2, one bad channel in between
                 % Example omitted: Ch5-Ch8  -> gap = 3, more than one missing channel
                 selected_chnls = sort(rmmissing(ch_numbering{i}));
                 if numel(selected_chnls) < 2
                     if ~silent
                         warning('Unable to perform bipolar re-referencing: fewer than 2 good contacts on standard lead.');
                     end
                     continue
                 end

                 chDiff = diff(selected_chnls);
                 validPairs = chDiff <= 2;

                 if ~any(validPairs)
                     if ~silent
                         warning('Unable to create bipolar pairs on standard lead: all gaps are larger than one missing contact.');
                     end
                     continue
                 end

                 ch_pairs = [selected_chnls(find(validPairs)).', ...
                     selected_chnls(find(validPairs) + 1).'];

                 lead_BP_ref = [lead_BP_ref, ...
                     inp_signal(:, ch_pairs(:,1)) - inp_signal(:, ch_pairs(:,2))];

                 BP_labels = [BP_labels, ...
                     arrayfun(@(a,b) sprintf('Ch%d-%d', a, b), ...
                     ch_pairs(:,1), ch_pairs(:,2), 'UniformOutput', false)'];

             case 'D'
                 pos = mod(selected_chnls, nch);

                 bottomContact = selected_chnls(pos == 1);      % bottom contact
                 lowerRing = selected_chnls(ismember(pos, 2:4));
                 upperRing = selected_chnls(ismember(pos, 5:7));
                 topContact = selected_chnls(pos == 0);      % top contact
                 ch_pairs = [];

                 % Bottom contact to all available lower-ring contacts
                 if ~isempty(bottomContact) && ~isempty(lowerRing)
                     ch_pairs = [ch_pairs; ...
                         repmat(bottomContact(:), numel(lowerRing), 1), lowerRing(:)];
                 end

                 % Matched vertical segmented pairs:
                 for seg = 1:3
                     lowerContact = selected_chnls(pos == seg + 1);   % 2,3,4
                     upperContact = selected_chnls(pos == seg + 4);   % 5,6,7
                     if ~isempty(lowerContact) && ~isempty(upperContact)
                         ch_pairs = [ch_pairs; lowerContact(:), upperContact(:)];
                     end
                 end

                 % Top contact to all available upper-ring contacts
                 if ~isempty(topContact) && ~isempty(upperRing)
                     ch_pairs = [ch_pairs; ...
                         repmat(topContact(:), numel(upperRing), 1), upperRing(:)];
                 end

                 if isempty(ch_pairs)
                     if ~silent
                         warning('Unable to perform bipolar re-referencing due to too many bad channels.');
                     end
                     continue
                 end

                 lead_BP_ref = [lead_BP_ref, ...
                     inp_signal(:, ch_pairs(:,1)) - inp_signal(:, ch_pairs(:,2))];
                 BP_labels = [BP_labels, arrayfun(@(a,b) sprintf('Ch%d-%d', a, b), ...
                     ch_pairs(:,1), ch_pairs(:,2), 'UniformOutput', false)'];
             otherwise
                 error('Invalid depth lead type "%s". Must be "S" or "D".', lead_types{i});
         end
    end
end
