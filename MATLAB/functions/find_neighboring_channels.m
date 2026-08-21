function neighbors = find_neighboring_channels(selectedChannel, construct, specification, method)
% Function that will take as an input specified construct and returns
% nearest neighbors for various electrode types and leads
% Dependencies:
%   - "surface_BP", "depth_BP" for BP re-referencing
%
% ~~~~~~~~~~~~~~~~~~ Chanel placement on 12 channel surface grids ~~~~~~~~~~~~~~~~
% ###### Grid1 ######## | ####### Grid2 ####### | ###### Grid3 ########
% ###### 1 # 7 ######## | ##### 13 # 18 ####### | ##### 23 # 28 #######
% ##### 2 ### 8 ####### | #### 14 ### 19 ###### | #### 24 ### 29 ######
% #### 3 ##### 9 ###### | ### 15 ##### 20 ##### | ### 25 ##### 30 #####
% ### 4 ####### 10 #### | ## 16 ####### 21 #### | ## 26 ####### 31 ####
% ## 5 ########## 11 ## | # 17 ########## 22 ## | # 27 ########## 32 ##
% #6 ############# 12 # |

% ~~~~~~~~~~~~~~~~~~ Chanel placement on 8 channel surface grids ~~~~~~~~~~~~~~~~
% #### Grid1 ### | #### Grid2 #######
% ### 1 ## 5 ### | ### 9  ## 13 #####
% ### 2 ## 6 ### | ### 10 ## 14 #####
% ### 3 ## 7 ### | ### 11 ## 15 #####
% ### 4 ## 8 ### | ### 12 ## 16 #####


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
%     Frederik Lampert, Mayo Clinic , 2026.
%--------------------------------------------------------------------------
    arguments
        selectedChannel (1,1) double {mustBeInteger, mustBePositive}
        construct (1,:) string {mustBeMember(construct, ["surface","depth","combined","surface_wide"])}
        specification {mustBeStringOrVectorOfInts} = []
        method (1,:) string {mustBeMember(method, ["NN","horizontal","vertical"])} = "NN"
    end

    switch construct
        case "surface"
            g1 = [(1:6)',  (7:12)'];
            g2 = [(13:17)', (18:22)'];
            g3 = [(23:27)', (28:32)'];
            leads = {g1; g2; g3};

        case "surface_wide"
            g1 = [(1:8)',  (9:16)'];
            g2 = [(17:24)', (25:32)'];
            leads = {g1; g2};

        case {"depth","combined"}
            leads = {1:8; 9:16; 17:24; 25:32};
    end

    switch construct
        case {"surface","surface_wide"}
            neighbors = surfaceNeighbors(selectedChannel, leads, method);

        case "depth"
            if isempty(regexp(char(specification), '^[DS]-[DS]-[DS]-[DS]$', 'once'))
                error(['Invalid DBS electrode type specification. ', ...
                    'Use a string such as "D-D-S-S".']);
            end

            leadTypes = strsplit(upper(char(specification)), '-');
            selLeadIdx = findLeadIndex(selectedChannel, leads);
            neighbors = depthNeighbors(selectedChannel, leads, leadTypes{selLeadIdx});

        case "combined"
            if isempty(regexp(char(specification), '^[GDS]-[GDS]-[GDS]-[GDS]$', 'once'))
                error(['Invalid combined electrode type specification. ', ...
                    'Use a string such as "G-G-S-D".']);
            end

            leadTypes = strsplit(upper(char(specification)), '-');
            selLeadIdx = findLeadIndex(selectedChannel, leads);
            leadType = leadTypes{selLeadIdx};

            switch leadType
                case 'G'
                    % For combined construct, treat this lead as a 4x2 surface grid.
                    leadVals = leads{selLeadIdx};
                    g = [leadVals(1:4)', leadVals(5:8)'];
                    neighbors = surfaceNeighbors(selectedChannel, {g}, method);

                case {'S','D'}
                    neighbors = depthNeighbors(selectedChannel, leads(selLeadIdx), leadType);

                otherwise
                    error('Invalid lead type "%s".', leadType);
            end
    end
    neighbors = unique(neighbors(:).');
end


function neighborChnls = surfaceNeighbors(selCh, leads, orientation)
    leadIdx = findLeadIndex(selCh, leads);
    grid = leads{leadIdx};
    [r, c] = find(grid == selCh, 1);

    if isempty(r)
        error('Selected channel %d was not found in the surface layout.', selCh);
    end

    neighborChnls = [];
    switch char(orientation)
        case 'NN'
            offsets = [-1 0; 1 0; 0 -1; 0 1];
        case 'vertical'
            offsets = [-1 0; 1 0];
        case 'horizontal'
            offsets = [0 -1; 0 1];
        otherwise
            error('Invalid surface orientation: %s', orientation);
    end

    for k = 1:size(offsets,1)
        rr = r + offsets(k,1);
        cc = c + offsets(k,2);
        if rr >= 1 && rr <= size(grid,1) && cc >= 1 && cc <= size(grid,2)
            neighborChnls(end+1) = grid(rr,cc); %#ok<AGROW>
        end
    end
end


function neighborChnls = depthNeighbors(selCh, ch_numbering, lead_type)
    leadIdx = findLeadIndex(selCh, ch_numbering);
    lead = ch_numbering{leadIdx};
    pos = find(lead == selCh, 1);   % local contact position 1-8

    if isempty(pos)
        error('Selected channel %d was not found in the depth layout.', selCh);
    end

    switch upper(char(lead_type))
        case 'S'
            localNeighbors = [pos-1, pos+1];
            localNeighbors = localNeighbors(localNeighbors >= 1 & localNeighbors <= numel(lead));
            neighborChnls = lead(localNeighbors);
        case 'D'
            switch pos
                case 1
                    localNeighbors = [2 3 4];
                case 2
                    localNeighbors = [1 3 4 5];
                case 3
                    localNeighbors = [1 2 4 6];
                case 4
                    localNeighbors = [1 2 3 7];
                case 5
                    localNeighbors = [2 6 7 8];
                case 6
                    localNeighbors = [3 5 7 8];
                case 7
                    localNeighbors = [4 5 6 8];
                case 8
                    localNeighbors = [5 6 7];
                otherwise
                    localNeighbors = [];
            end
            neighborChnls = lead(localNeighbors);
        otherwise
            error('Invalid depth lead type "%s". Must be "S" or "D".', lead_type);
    end
end


function leadIdx = findLeadIndex(ch, leads)
    leadIdx = [];
    for i = 1:numel(leads)
        if any(leads{i}(:) == ch)
            leadIdx = i;
            return
        end
    end
    error('Selected channel %d was not found in any lead.', ch);
end


function mustBeStringOrVectorOfInts(x)
    if isempty(x)
        return
    end
    if ischar(x) || isstring(x)
        return
    elseif isnumeric(x)
        if any(mod(x,1) ~= 0)
            error('specification must be a string or vector of integers.');
        end
    else
        error('specification must be a string or vector of integers.');
    end
end