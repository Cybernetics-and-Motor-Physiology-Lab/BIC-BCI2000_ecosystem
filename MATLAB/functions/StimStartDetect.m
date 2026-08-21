function stim_idx = StimStartDetect(signal, states, parameters, window, method, plot_detection)
    
    stim_start = find(diff(double(states.ImplantStimulation)) == 1);       % states.ImplantStimulation in uint8 (allowing only positive int)
    stim_end = find(diff(double(states.ImplantStimulation)) == -1);
    stim_start = stim_start+floor((stim_end-stim_start)/2);

    window = round(window/2,3);
    t1_wnd = stim_start - window;   t1_wnd(t1_wnd<1)=1;                % Calculate the window beggining index (and take care of border values)
    t2_wnd = stim_start + window;   t2_wnd(t2_wnd>length(signal))=length(signal);  % Calculate the window end index (and take care of border values)

    stim_idx = nan(size(stim_start));

    for wnd = 1:length(stim_start)      % iterate through the stimulation windows
        stim_wnd = signal(t1_wnd(wnd):t2_wnd(wnd)); % take the appropriate signal window

        switch method
            case 'center'
                stim_idx = stim_start;
                
            case 'tkeo'
                % ~~~~~ Teager Keiser Energy Operator ~~~~~~~~
                % https://www.mathworks.com/matlabcentral/fileexchange/45406-teager-keiser-energy-operator-vectorized)
                y=[0;diff(stim_wnd)];
                squ=y(2:length(y)-1).^2;
                oddi=y(1:length(y)-2);
                eveni=y(3:length(y));
                ey=squ - (oddi.*eveni);
                % [x[n]] = x^2[n] - x[n - 1]x[n + 1]
                % operator ex
                squ1=stim_wnd(2:length(stim_wnd)-1).^2;
                oddi1=stim_wnd(1:length(stim_wnd)-2);
                eveni1=stim_wnd(3:length(stim_wnd));
                ex=squ1 - (oddi1.*eveni1);
                ex = [ex(1); ex; ex(length(stim_wnd)-2)]; %make it the same length

                % ~~~~ Detect Peaks ~~~~~~~~~
                [pks,loc] = findpeaks(ex,'NPeaks',1, 'SortStr','descend');
                stim_idx(wnd) = t1_wnd(wnd) + loc-1;
                % stim_idx(wnd) = t1_wnd(wnd) + find(ex==max(ex))-1;
                % figure
                % plot(stim_wnd); hold on
                % stem(loc,stim_wnd(loc))
            
            case 'min'
                stim_idx(wnd) = t1_wnd(wnd) + find(stim_wnd==min(stim_wnd))-1;

            case 'max'
                stim_idx(wnd) = t1_wnd(wnd) + find(stim_wnd==max(stim_wnd))-1;

        end
    end

    if plot_detection == true
        figure('Position', [100 100 1500 800]);
        plot(linspace(0,length(signal), length(signal)), signal); hold on
        stem(stim_start,signal(stim_start), '--o','filled', 'Color','#77AC30' ); hold on 
        stem(stim_idx, signal(stim_idx),'--*','Color','r') %LineWidth=1
        for i=1:length(stim_idx)
            text(stim_idx(i)+50,signal(stim_idx(i))+70, num2str(i))
        end
        plot((double(states.ImplantStimulation)*0.25*max(signal)+min(signal)))
        lgd = legend({'Channel for detection', 'Center of ImplantStimulation State', 'Detected Stimulation', 'ImplantStimulation'});
        fontsize(lgd,14,'points');
        xlim([0,length(signal)])
        title('Stim detection')
    end
end




    % % estimate trend using polyfit
    % poly_est = polyfit(t,ch_resp_corrected(resp_idx,:),4);
    % trend_est = polyval(poly_est,t);
    % ch_resp_corrected(resp_idx,:) = ch_resp_corrected(resp_idx,:)-trend_est;
    % plot(t,ch_resp(resp_idx,:),t,trend_est,t,ch_resp_corrected(resp_idx,:));
    % legend('data','estimated trend','trend removed','Location','NorthWest');