function wiggledisplay(Data,x,t,style,showmax,wscale,axisorder)
%
% WIGGLEDISPLAY : Wiggle/variable area/image or combined plot of GPR (or
%                 seismic) data.Also has the capability of projecting a
%                 colour-scale image of the data beneath the wiggles.   
%
%         Usage : wiggledisplay(Data,dmax,t,x,style,showmax);  
%
%       Inputs  :
%        Data   : [ns,ntr ] matrix of the GPR or seismic section 
%           x   : [1, ntr]  vector of the hozirontal (spatial) trace
%                           coordinates 
%           t   : [1, ns]   vector of the time coordinates - MUST BE ROW
%                           VECTOR for wiggledisplay to work properly.
%        style  : 'vararea' creates variable area plots 
%                 'wiggle'  creates wiggle-trace plots
%                 'wig&img' creates wiggle / image plots (shows data
%                           image beneath wiggles)
%                 'var&img' creates variable area / image plots (shows
%                           data image beneath wiggles)
%      showmax  : scalar, max number of traces to display. Default = 100.
%       wscale  : scaling factor, can be set empty!. Optimize by trial and
%                 error 
%
%      Cretits  : Some programming tips taken from the M-code of an
%                 unknown author, presumably Thomas Mejer Hansen of the
%                 Niels Bohr Institute for Astronomy, Physics and
%                 Geophysics, University of Copenhagen, Denmark.
%
%        Author : Andreas Tzanis, 
%                 Department of Geophysics, 
%                 University of Athens
%                 atzanis@geol.uoa.gr
%
% Copyright (C) 2005, Andreas Tzanis. All rights reserved.
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%
global ENVAR
plImage = 0;              % Default is: No image under wiggles!
if nargin == 1,
    x=[1:1:size(Data,2)];
    t=[1:1:size(Data,1)];
    style  ='wiggle';
    showmax=100;
    wscale = 1;
    axisorder=1;
end
% max data value to be used in constructing wiggle scale factor
dmax=max(nanmax(abs(Data)));
if nargin == 2,
    t=[1:1:size(Data,1)];
    style='wiggle';
    showmax=100;
    wscale = 1;
    axisorder=1;
end
if nargin == 3, 
    style='wiggle'; 
    showmax=100;
    wscale = 1;
    axisorder=1;
end
if nargin == 4, 
    showmax=100;
    wscale = 1;
	axisorder=1;
end
if nargin == 5,
    wscale = 1;
    axisorder=1;
end
if nargin == 6,
    axisorder=1;
end
if strcmp(style,'wig&img'),
    style = 'wiggle';
    plImage = 1;          % Draw image under wiggles
end
if strcmp(style,'var&img'),
    style = 'vararea';    % Draw image under var. area plot
    plImage = 1;
end
if plImage == 1,
    imagesc(x,t,Data);
    colormap(ENVAR.colormap);
    hold on
end
% force t to be a row vector
[nr,~] = size(t);
if nr > 1,
    t=t';
end
if (showmax > 0),
    dx = x(2) - x(1);
    ntraces=length(x);
    d=ntraces/showmax;
    if d <= 1; 
        d = 1; 
    end
    d=round(d);
    dmax=dmax/(wscale*d*abs(dx));
    for i=1:d:ntraces
        if i > ntraces,
            break
        end
        xt=Data(:,i)'./dmax;
        if strcmp(style,'vararea')==1,
            xt1=xt;
            xt1(xt1 < 0) = 0;
            fill(x(i) + [xt,fliplr(xt1)],[t,fliplr(t)],[0 0 0]);
        else
            plot(x(i) + xt, t ,'k-','linewidth',.1);
        end   
        if i==1, 
            hold on;
        end
    end
    if axisorder == 1
        axis([min(x(1:ntraces))-(x(2)-x(1)) max(x(1:ntraces))+(x(2)-x(1)) min(t) max(t)])
    else
        axis([min(x) max(x) 0 max(t)])
    end
end
hold off;  
set(gca,'Ydir','revers')
set(gca,'XDir','reverse');
view(-90,90)
% End function wiggledisplay
