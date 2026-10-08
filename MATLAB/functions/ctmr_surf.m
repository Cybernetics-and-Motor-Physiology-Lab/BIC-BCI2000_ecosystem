function bsurf=ctmr_surf(cortex, alp, clr_map)
% function [electrodes]=ctmr_gauss_plot(cortex,electrodes,weights)
% projects electrode locations onto their cortical spots in the 
% left hemisphere and plots about them using a gaussian kernel
% for only cortex use: 
% ctmr_gauss_plot(cortex,[0 0 0],0)
% rel_dir=which('loc_plot');
% rel_dir((length(rel_dir)-10):length(rel_dir))=[];
% addpath(rel_dir)
%   Created by:
%   K.J. Miller & D. Hermes, 
%   Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%   Version 1.1.0, released 26-11-2009

% Set default values for plotting



if nargin < 2
    alp = 0.7;  %transparency
    clr_map = 'gray';   % color map
end
if nargin < 3
    clr_map = 'gray';   % color map
end

c=zeros(length(cortex.vert(:,1)),1).'; % set all to the same color

%%
% c=(c/max(c));
bsurf=tripatch(cortex, 'nofigure', c');
shading interp;
a=get(gca);
%%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
d=a.CLim;
set(gca,'CLim',[-4 4])
l=light;
colormap(clr_map)
lighting gouraud; %play with lighting...
% material dull;
material([.3 .8 .1 10 1]);
axis off
set(gcf,'Renderer', 'zbuffer')

% view(45,45);
set(l,'Position',[1 1 1])
%
set(bsurf,'FaceAlpha',alp)
set(gcf,'Color','w')

% %exportfig
% exportfig(gcf, strcat(cd,'\figout.png'), 'format', 'png', 'Renderer', 'painters', 'Color', 'cmyk', 'Resolution', 600, 'Width', 4, 'Height', 3);
% disp('figure saved as "figout"');