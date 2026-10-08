function rb_scaled_circle_add_popout(locs, wts, th, phi)
    

% note - it will make a small white circle if wts is 0

    msize=6; 

    wts(isnan(wts)) = 0;
    
%% scale to maximum across slices
    wm=max(abs(wts));

%% offset locs in direction of view
a_offset=.1*max(abs(locs(:,1)))*[cosd(th-90)*cosd(phi) sind(th-90)*cosd(phi) sind(phi)];

locs=locs+...
    repmat(a_offset,size(locs,1),1);

locs(isnan(locs))=0;

%% plot weights
    hold on
%     border_grayscale=.98;
    border_grayscale=0;    
    for q=1:size(locs,1)% add activity colorscale
        %
        if abs(wts(q))==0 % not significant
            % circle w border
            plot3(locs(q,1),locs(q,2),locs(q,3),'o',...
            'MarkerSize',msize,...
            'LineWidth',.5,...
            'MarkerEdgeColor',border_grayscale*[1 1 1],... 
            'MarkerFaceColor',.35*[1 1 1])  
            %
        elseif wts(q)>0
            % white circle
            plot3(locs(q,1),locs(q,2),locs(q,3),'o',...
            'MarkerSize',msize*abs(wts(q))/wm+msize,...
            'LineWidth',.5,...
            'MarkerEdgeColor',border_grayscale*[1 1 1],... 
            'MarkerFaceColor',.99*[1 1-wts(q)/wm 1-wts(q)/wm])
            %
        elseif wts(q)<0
            % white circle
            plot3(locs(q,1),locs(q,2),locs(q,3),'o',...
            'MarkerSize',msize*abs(wts(q))/wm+msize,...
            'LineWidth',.5,...
            'MarkerEdgeColor',border_grayscale*[1 1 1],... 
            'MarkerFaceColor',.99*[1+wts(q)/wm 1+wts(q)/wm 1])
        end
    end
    hold off

        
%% rotate view
loc_view(th,phi)

