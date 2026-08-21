function plotReCap(plotDataStruct)
    % plotRecordingWaterfall - 3D plot of RMS per channel across time,
    % color-coded by impedance.

    if isempty(plotDataStruct)
        error('Input plotDataStruct is empty.');
    end

    dates = [plotDataStruct.Date];
    nDays = numel(dates);
    nCh = 32;

    % Preallocate matrices
    RMS = nan(nCh, nDays);       % [channel x day]
    Imp = nan(nCh, nDays);       % [channel x day]

    for i = 1:nDays
        RMS(:, i) = plotDataStruct(i).RMS(:);
        Imp(:, i) = plotDataStruct(i).Impedance(:);
    end
    
    [Y, X] = meshgrid(1:nCh, datenum(dates));  % meshgrid in [rows x cols]
    Z = RMS';  % [days x channels]
    C = Imp';  % same size as Z

    figure('Name', 'Recording Waterfall: RMS and Impedance','Position', [249 127 1259 816]);
    plt = mesh(X, Y, Z, C, 'LineWidth',3.5, 'FaceAlpha',0.4); %'FaceAlpha',0.6; 'FaceColor','none'
    plt.MeshStyle = 'column';

    colormap(turbo);
    clim([100 5e3]);
    cb = colorbar;
    cb.Label.String = 'Impedance [Ω]';

    datetick('x', 'dd-mmm', 'keeplimits');
    xlim('tight')
    xlabel('Recording Date');
    ylabel('Channel');
    zlabel('RMS Power');
    title('Recording Capabilities: RMS Power & Impedance');
    grid on;
    view(15, 30);
    fontsize(16,'points')

    %% ---------------------------
    % [X, Y] = meshgrid(datenum(dates), 1:32);
    % Z = RMS;
    % C = Imp;  % same size as Z
    % 
    % figure;
    % plt = waterfall(X, Y, Z, C);
    % plt.LineWidth = 1.5;
    % plt.FaceAlpha = 0.7;
    % datetick('x', 'dd-mmm', 'keeplimits');
    % colorbar
    % colormap(turbo);
    % xlabel('Date');
    % ylabel('Channel');
    % zlabel('RMS Power');
    % title('Recording Capabilities: Waterfall View');
    %  view(10, 60);

end