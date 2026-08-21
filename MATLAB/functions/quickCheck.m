function quickCheck(chnls, file_path)

switch nargin
    case 2
        % --------------- load data from path -------------
        [signal, states, parameters ] ...
            = load_bcidat(file_path, '-calibrated');         % load using MEX package
        % https://www.bci2000.org/mediawiki/index.php/User_Reference:Matlab_MEX_Files
        [tmp1,fl_id,tmp2] = fileparts(file_path);
        clear tmp1 tmp 2
    case 1 
        % ---------- Select file in GUI ---------------
        [fl, fl_pth] = uigetfile([['/Users' filesep] '*.dat'], 'Select a File');     % select file
        [signal, states, parameters] = load_bcidat([fl_pth, fl], '-calibrated');         % load using MEX package
        % https://www.bci2000.org/mediawiki/index.php/User_Reference:Matlab_MEX_Files
        % get the name of the file
        fl_id = strsplit(fl, '.'); fl_id = fl_id(1);       % get file name and extension

    otherwise
        error('Select channels to plot')
end


%% Plot data

FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
t_ax = linspace(0,size(signal,1)/FS, size(signal,1))';   % Creating time axis for plotting signal

% conversion to table
varnames = cell(1,size(signal,2));
for i=1:size(signal,2)
    varnames(i) = {sprintf('Chanel %d', i)};
end

% Plot data
figure
    if parameters.EnableStimulation.NumericValue == 1          % check if there was stimulatio
        sig_tab = array2table([signal, t_ax, double(states.ImplantStimulation)], 'VariableNames', [varnames {'Time'} {'Stimulation'}]);
        s = stackedplot(sig_tab(:,[chnls size(sig_tab,2)-1:size(sig_tab,2)]),"XVariable",'Time');
    else            % if no then don't dispay stimulation states
        sig_tab = array2table([signal, t_ax], 'VariableNames', [varnames {'Time'}]);
        s = stackedplot(sig_tab(:,[chnls size(sig_tab,2)]),"XVariable",'Time');
    end
xlabel('Time [s]'); 
sgtitle(fl_id, 'Interpreter', 'none')

end
