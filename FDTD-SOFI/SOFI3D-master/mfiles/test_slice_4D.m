close all; clear all;    
[x,y,z] = meshgrid(0:.5:10,0:.5:10,0:.5:10);
    c = x.^2+y.^2+z.^2;
    xs = 0:0.5:10;
    ys = xs;
    zs = xs;
    c(7:15,7:15,13:21)=NaN;
    h = slice(x,y,z,c,xs,ys,zs);
    set(h,'FaceColor','interp',...
        'EdgeColor','none')
    camproj perspective
    box on
    view(-70,70)
    colormap hsv
    colorbar