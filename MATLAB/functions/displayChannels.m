function displayChannels(signal, sampling_frequency, reference_channel, channels_to_plot, amp_multiplier, fig_title)

    if nargin < 4
        channels_to_plot = 1:size(signal,2);
    elseif isempty(channels_to_plot)
        channels_to_plot = 1:size(signal,2);
    else
        channels_to_plot = sort(channels_to_plot);
    end
    
    if nargin<5
        amp_multiplier = 1;
    end
    
    
    nChan = length(channels_to_plot);
    N = size(signal,1);
    t_ax = linspace(0,N/sampling_frequency,N);
    
    y_spacing = 10000;
    if mod(nChan,2)==1
        shift_s=linspace(0-floor(nChan/2),0+floor(nChan/2),nChan)*y_spacing;
    else
        shift_s=linspace(0-floor(nChan/2),0+floor(nChan/2)-1,nChan)*y_spacing;
    end
    signal = signal(:,channels_to_plot)*amp_multiplier;
    shifted_signal = signal-repmat(shift_s,N,1);
    
    f = figure('Position', [100 100 1500 800]);
    main_ax = axes('Parent', f);
    hold(main_ax, 'on')

    for ch = 1:length(channels_to_plot)
        if ismember(channels_to_plot(ch), reference_channel)
            plt = plot(main_ax,t_ax,shifted_signal(:,ch),'Color', '#f29220', ...
                'DisplayName', sprintf('Ch %d',ch)); hold on
            % text(t_ax(end-2), shifted_signal(end,ch),'Reference ch.','Color', 'r','FontSize',13);
            % text(xLabelPos, shifted_signal(end, ch), ...
            %     sprintf('Ch %d', channels_to_plot(ch)), ...
            %     'FontSize', 12, 'HorizontalAlignment', 'right', ...
            %     'VerticalAlignment', 'middle', 'Color', 'r');
            plt.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Channel ',repelem(channels_to_plot(ch),N));
        else
            plt = plot(main_ax, t_ax,shifted_signal(:,ch), ...
                'DisplayName', sprintf('Ch %d',ch), 'Color','#2b2b2b'); hold on %'Color', [0 0.4470 0.7410]
            % text(xLabelPos, shifted_signal(end, ch), ...
            %     sprintf('Ch %d', channels_to_plot(ch)), ...
            %     'FontSize', 12, 'HorizontalAlignment', 'right', ...
            %     'VerticalAlignment', 'middle', 'Color', 'k');
            plt.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Channel ',repelem(channels_to_plot(ch),N));
        end
    end
    
    
    % Compute major and minor tick locations
    yticks_major = flip(shift_s) * (-1);
    % yticks_minor = sort([yticks_major - 50*amp_multiplier, yticks_major + 50*amp_multiplier]);
    ylimits = [min(shifted_signal(:))-y_spacing/2, max(shifted_signal(:))+y_spacing/2];

    % LEFT Y-axis: Channel names
    main_ax.YLim = ylimits;
    main_ax.YTick = yticks_major;
    main_ax.YTickLabel = compose('Ch %d', flip(channels_to_plot));
    main_ax.YLabel.String = 'Channel';

    % X-axis
    main_ax.XLim = [min(t_ax), max(t_ax)];
    main_ax.XLabel.String = 'Time [s]';    
    
    if exist('fig_title', 'var')
        title(fig_title,'Interpreter','none');
    else
        title('Raw signal');
    end
    
    grid on;
    fontsize(16, "points")

    % -------- Draw scale bar (bottom right corner) ------
    % Define scale sizes
    scale_amp = 100 * amp_multiplier;
    scale_time = 1;  % seconds

    % Position: bottom-right, just below lowest trace
    x0 = main_ax.XLim(2) - scale_time;
    y0 = min(shifted_signal(:)) - 0.3 * y_spacing;

    % Vertical line (amplitude scale)
    line(main_ax, [x0, x0], [y0, y0 + scale_amp], 'Color', 'k', 'LineWidth', 2);
    text(main_ax, x0 + 0.05 * scale_time, y0 + scale_amp/2, '100 \muV', ...
        'FontSize', 10, 'VerticalAlignment', 'middle');

    % Horizontal line (time scale)
    line(main_ax, [x0, x0 + scale_time], [y0, y0], 'Color', 'k', 'LineWidth', 2);
    text(main_ax, x0 + scale_time/2, y0 - 0.05 * y_spacing, ...
        sprintf('%.1f s', scale_time), ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
end



%% Old code
    % f = figure('Position', [100 100 1500 800]); 
    % tlo = tiledlayout(nChan,1, TileSpacing="none");
    % h = gobjects(1,nChan);
    % for ch = 1:nChan
    %     ax = nexttile(tlo);
    %     if ch == reference_channel
    %         h(ch) = plot(t_ax,signal(:,ch),'Color', 'r','DisplayName', sprintf('Ch %d',ch));
    %         title(['Reference Channel' num2str(reference_channel)]);
    %     else
    %         h(ch) = plot(t_ax,signal(:,ch),'DisplayName', sprintf('Ch %d',ch));
    %     end
    %     xlim([0 length(t_ax)/sampling_frequency]);
    %     yticks([0]);
    % 
    %     set(ax,'box', 'off');
    % 
    % end
    % xlabel('Time [s]'); ylabel('U [uV]');
    % title(tlo,'Channels');

    % ~~~~

        % xLabelPos = t_ax(end) + 0.04 * range(t_ax);  % Slightly right of the plot;  % Left of plot area
        % plot([N/sampling_frequency, N/sampling_frequency], ...
    %     [signal(end,end)-(50*amp_multiplier), ...
    %     signal(end,end)+(50*amp_multiplier)], ...
    %     '-_','Color','k','LineWidth',1.5);
    % text(N/sampling_frequency+0.1, signal(end,end),'100\muV', 'FontSize',13)
