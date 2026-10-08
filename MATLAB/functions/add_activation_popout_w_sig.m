function add_activation_popout_w_sig(locs,bad_channels,r2,pvals,theta,phi,colorMap)

% Modified "rb_scaled_circle_add_popout_w_sig" function from KJM
% note - it will make a small white circle if wts is 0

msize=10;

good_chnls = 1:length(locs);        % Create vector of all channel numbers
good_chnls(bad_channels) = [];      % Specify only good channels

r2(isnan(r2)) = 0;

%% scale to maximum across slices
wm=max(abs(r2));

%% offset locs in direction of view
a_offset=.1*max(abs(locs(:,1)))*[cosd(theta-90)*cosd(phi) sind(theta-90)*cosd(phi) sind(phi)];
locs = locs + repmat(a_offset,size(locs,1),1);
g_locs = locs(good_chnls,:);        % Get locations of good channels

%% plot weights
hold on
%     border_grayscale=.98;
border_grayscale=0;


for q=1:size(g_locs,1)% add activity colorscale
    if pvals(q)==0 % doesn't meet p-threshold
        if abs(r2(q))==0 % not significant
            % circle w border
            plot3(g_locs(q,1),g_locs(q,2),g_locs(q,3),'o',...
                'MarkerSize',msize,...
                'LineWidth',1,...                'MarkerEdgeColor',border_grayscale*[1 1 1],...
                'MarkerEdgeColor','none',...
                'MarkerFaceColor',.35*[1 1 1])
            %
        elseif r2(q)>0
            % white circle
            plot3(g_locs(q,1),g_locs(q,2),g_locs(q,3),'o',...
                'MarkerSize',msize*abs(r2(q))/wm+msize,...
                'LineWidth',1,...
                'MarkerEdgeColor','none',... 'MarkerEdgeColor',border_grayscale*[1 1 1],...
                'MarkerFaceColor',.99*[1 1-r2(q)/wm 1-r2(q)/wm])
            %
        elseif r2(q)<0
            % white circle
            plot3(g_locs(q,1),g_locs(q,2),g_locs(q,3),'o',...
                'MarkerSize',msize*abs(r2(q))/wm+msize,...
                'LineWidth',1,...
                'MarkerEdgeColor','none',... 'MarkerEdgeColor',border_grayscale*[1 1 1],...
                'MarkerFaceColor',.99*[1+r2(q)/wm 1+r2(q)/wm 1])
        end

    elseif pvals(q)==1 % it is statistically significant
        if r2(q)>0
            % white circle
            plot3(g_locs(q,1),g_locs(q,2),g_locs(q,3),'o',...
                'MarkerSize',msize*abs(r2(q))/wm+msize,...
                'LineWidth',2,...
                'MarkerEdgeColor',border_grayscale*[1 1 1],... 'MarkerEdgeColor',[1 1 0],...
                'MarkerFaceColor',.99*[1 1-r2(q)/wm 1-r2(q)/wm])
            %
        elseif r2(q)<0
            % white circle
            plot3(g_locs(q,1),g_locs(q,2),g_locs(q,3),'o',...
                'MarkerSize',msize*abs(r2(q))/wm+msize,...
                'LineWidth',2,...
                'MarkerEdgeColor',border_grayscale*[1 1 1],... 'MarkerEdgeColor',[1 1 0],...
                'MarkerFaceColor',.99*[1+r2(q)/wm 1+r2(q)/wm 1])
        end
    end
end

% Plot bad channels with "X" markers
plt = plot3(locs(bad_channels,1),locs(bad_channels,2),locs(bad_channels,3),'x',...
    'MarkerSize',msize,...
    'LineWidth', 2,...
    'MarkerEdgeColor','r');
% colormap(plt,colorMap)
% colorbar(plt);

hold off
      
%% rotate view
loc_view(theta,phi)

