function S = struct2single(S)
% Recursively convert all numeric fields in a struct to single

    if isstruct(S)
        fn = fieldnames(S);
        for i = 1:numel(fn)
            S.(fn{i}) = struct2single(S.(fn{i}));
        end

    elseif isnumeric(S)
        S = single(S);

    elseif iscell(S)
        for i = 1:numel(S)
            S{i} = struct2single(S{i});
        end
    end
end