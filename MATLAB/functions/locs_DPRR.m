function [dp_channels, dp_locs] = locs_DPRR(locs, dist, plotfig)

% This function performs bipolar re-refencing of the channels. 
%
% KJM, CaMP lab, 2022

%% Calculate contact distance

sep_locs = sum(diff(locs).^2, 2) .^ .5;
dp_channels = find(sep_locs < dist);

%% Differential pairs are located between contacts

dp_locs = (locs(1:(end-1), :) + locs(2:end, :)) / 2;
dp_locs = dp_locs(dp_channels, :);
    
%% Plot figures showing retained and excluded differential pairs

if plotfig == 'y'
    figure, plot(sep_locs, 'ro'), hold on, plot(dp_channels, sep_locs(dp_channels, :),'go'),
    xlabel('electrode difference pairs'), ylabel('distance between pair')
    %
    figure, plot3(locs(:, 1),locs(:, 2),locs(:, 3), 'b.')
    hold on, plot3(dp_locs(:, 1), dp_locs(:, 2), dp_locs(:, 3),'r.')
end

end




