%creates a binary formation model for 3D seismic modeling
close all;
clear all;
clc;

%general parameters
writefiles=1;%1=yes else=no
modname='../par/model/2layer_TQ';
plotmodels=1;%1=yes else=no
grid_spacing=0.1;

%dimensions of base model
xmax_base=200;%gridpoints
ymax_base=100; %(vertical)
zmax_base=50;

x=grid_spacing:grid_spacing:xmax_base*grid_spacing;
y=grid_spacing:grid_spacing:ymax_base*grid_spacing;
z=grid_spacing:grid_spacing:zmax_base*grid_spacing;

%change in material properties, i.e. fault zone
x_start_fault=1; % begin of new material, at fault zone in grid points
y_start_fault=50; % 
z_start_fault=1; % 

%velocity, density base model
vp_base=400; %m/s
vs_base=200;
rho_base=1800; %kg/m3

%velocity, density fault zone
vp_fault=800.0; %m/s
vs_fault=400.0;
rho_fault=2000.0; %g/cm3


%----------------VP model
%creating model
mod=(zeros(xmax_base,ymax_base,zmax_base));
%filling in base parmaters
mod(:,:)=vp_base;
%filling in fault parameters (i.e. change in material properties)
mod(x_start_fault:end,y_start_fault:end,z_start_fault:end)=vp_fault;
clim = [min(min(min(mod))) max(max(max(mod)))];
%plot models

top_margin = 0.1; % top margin
btm_margin = 0.1; % bottom margin
left_margin = 0.1;% left margin
right_margin = 0.15;% right margin
 
fig_margin = 0.08; % margin beween figures(sub) 
 
row = 2; % rows
col = 2; % cols
 
% Calculate figure height and width according to rows and cols 
fig_h = (1- top_margin - btm_margin - (row-1) * fig_margin) / row;
fig_w = (1 - left_margin - right_margin - (col-1) * fig_margin) / col;

if plotmodels==1
    h1=figure(1);
    for i = 1 : row
        for j = 1 : col
            % figure position: you can refer to 'help axes' to review the
            % parameter meaning, note that original point is lower-left
            position = [left_margin + (j-1)*(fig_margin+fig_w), ...
                1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
                fig_w, fig_h];
            axes('position', position)
            % draw colorful pictures...
            if i==1 && j== 1
                imagesc(z,y,squeeze(mod(1,:,:)),clim);
                shading flat;
                title('vp-velocity Model (y-z-plain)');
                xlabel('z in m (horizontal)');
                ylabel('y in m (vertical)');
                set(gca,'YDir','reverse');
                axis equal tight;
            elseif  i==1 && j== 2
                imagesc(x,z,squeeze(mod(:,1,:))',clim);
                shading flat;
                title('vp-velocity Model (x-z-plain)');
                xlabel('x in m');
                ylabel('z in m (horizontal)');
                axis equal tight;
            elseif  i==2 && j== 1
                imagesc(x,y,squeeze(mod(:,:,1))',clim);
                shading flat;
                title('vp-velocity Model (x-y-plain)');
                xlabel('x in m');
                ylabel('y in m (vertical)');
                set(gca,'YDir','reverse');
                axis equal tight;
            else
                [z_grid,x_grid,y_grid] = meshgrid(z,x,y);
                mod_T=permute(mod,[1 3 2]); % transpose the 3D matrix
                h = slice(z_grid,x_grid,y_grid,mod_T,z,x,y);
                set(h,'FaceColor','interp',...
                    'EdgeColor','none')
                camproj perspective
                box on
                view(50,30)
                title('vp-velocity Model (3D)');
                xlabel('z in m');
                ylabel('x in m');
                zlabel('y in m (vertical)');
                set(gca,'ZDir','reverse');
                axis equal tight;
                
            end
        end
        
        % draw colorbar
        axes('position', [1-right_margin-fig_margin, btm_margin, 0.2,...
            (1-(top_margin+btm_margin))*2/5]);
        axis off; 
        colorbar();caxis(clim);
    end
end

%write models to disk
if writefiles==1
    mod=permute(mod,[2 1 3]);
    fid_mod=fopen(strcat(modname,'.vp'),'w');
    fwrite(fid_mod,mod,'float');
    fclose(fid_mod);
end

clear i j h mod_T mod x_grid y_grid z_grid clim;

%----------------VS model
%creating model
mod=(zeros(xmax_base,ymax_base,zmax_base));
%filling in base parmaters
mod(:,:)=vs_base;
%filling in fault parameters (i.e. change in material properties)
mod(x_start_fault:end,y_start_fault:end,z_start_fault:end)=vs_fault;
clim = [min(min(min(mod))) max(max(max(mod)))];

%plot models
if plotmodels==1
    h2=figure(2);
    for i = 1 : row
        for j = 1 : col
            % figure position: you can refer to 'help axes' to review the
            % parameter meaning, note that original point is lower-left
            position = [left_margin + (j-1)*(fig_margin+fig_w), ...
                1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
                fig_w, fig_h];
            axes('position', position)
            % draw colorful pictures...
            if i==1 && j== 1
                imagesc(z,y,squeeze(mod(1,:,:)),clim);
                shading flat;
                title('vs-velocity Model (y-z-plain)');
                xlabel('z in m (horizontal)');
                ylabel('y in m (vertical)');
                set(gca,'YDir','reverse');
                axis equal tight;
            elseif  i==1 && j== 2
                imagesc(x,z,squeeze(mod(:,1,:))',clim);
                shading flat;
                title('vs-velocity Model (x-z-plain)');
                xlabel('x in m');
                ylabel('z in m (horizontal)');
                axis equal tight;
            elseif  i==2 && j== 1
                imagesc(x,y,squeeze(mod(:,:,1))',clim);
                shading flat;
                title('vs-velocity Model (x-y-plain)');
                xlabel('x in m');
                ylabel('y in m (vertical)');
                set(gca,'YDir','reverse');
                axis equal tight;
            else
                [z_grid,x_grid,y_grid] = meshgrid(z,x,y);
                mod_T=permute(mod,[1 3 2]); % transpose the 3D matrix
                h = slice(z_grid,x_grid,y_grid,mod_T,z,x,y);
                set(h,'FaceColor','interp',...
                    'EdgeColor','none')
                camproj perspective
                box on
                view(50,30)
                title('vs-velocity Model (3D)');
                xlabel('z in m');
                ylabel('x in m');
                zlabel('y in m (vertical)');
                set(gca,'ZDir','reverse');
                axis equal tight;
                
            end
        end
        
        % draw colorbar
        axes('position', [1-right_margin-fig_margin, btm_margin, 0.2,...
            (1-(top_margin+btm_margin))*2/5]);
        axis off; 
        colorbar();caxis(clim);
    end
end

%write models to disk
if writefiles==1
    mod=permute(mod,[2 1 3]);
    fid_mod=fopen(strcat(modname,'.vs'),'w');
    fwrite(fid_mod,mod,'float');
    fclose(fid_mod);
end

clear i j h mod_T mod x_grid y_grid z_grid clim;

%----------------Rho model
%creating model
mod=(zeros(xmax_base,ymax_base,zmax_base));
%filling in base parmaters
mod(:,:)=rho_base;
%filling in fault parameters (i.e. change in material properties)
mod(x_start_fault:end,y_start_fault:end,z_start_fault:end)=rho_fault;
clim = [min(min(min(mod))) max(max(max(mod)))];

%plot models
if plotmodels==1
    h3=figure(3);
    for i = 1 : row
        for j = 1 : col
            % figure position: you can refer to 'help axes' to review the
            % parameter meaning, note that original point is lower-left
            position = [left_margin + (j-1)*(fig_margin+fig_w), ...
                1- (top_margin + i * fig_h + (i-1) * fig_margin), ...
                fig_w, fig_h];
            axes('position', position)
            % draw colorful pictures...
            if i==1 && j== 1
                imagesc(z,y,squeeze(mod(1,:,:)),clim);
                shading flat;
                title('Density Model (y-z-plain)');
                xlabel('z in m (horizontal)');
                ylabel('y in m (vertical)');
                set(gca,'YDir','reverse');
                axis equal tight;
            elseif  i==1 && j== 2
                imagesc(x,z,squeeze(mod(:,1,:))',clim);
                shading flat;
                title('Density Model (x-z-plain)');
                xlabel('x in m');
                ylabel('z in m (horizontal)');
                axis equal tight;
            elseif  i==2 && j== 1
                imagesc(x,y,squeeze(mod(:,:,1))',clim);
                shading flat;
                title('Density Model (x-y-plain)');
                xlabel('x in m');
                ylabel('y in m (vertical)');
                set(gca,'YDir','reverse');
                axis equal tight;
            else
                [z_grid,x_grid,y_grid] = meshgrid(z,x,y);
                mod_T=permute(mod,[1 3 2]); % transpose the 3D matrix
                h = slice(z_grid,x_grid,y_grid,mod_T,z,x,y);
                set(h,'FaceColor','interp',...
                    'EdgeColor','none')
                camproj perspective
                box on
                view(50,30)
                title('Density Model (3D)');
                xlabel('z in m');
                ylabel('x in m');
                zlabel('y in m (vertical)');
                set(gca,'ZDir','reverse');
                axis equal tight;
                
            end
        end
        
        % draw colorbar
        axes('position', [1-right_margin-fig_margin, btm_margin, 0.2,...
            (1-(top_margin+btm_margin))*2/5]);
        axis off; 
        colorbar();caxis(clim);
    end
end

%write models to disk
if writefiles==1
    mod=permute(mod,[2 1 3]);
    fid_mod=fopen(strcat(modname,'.rho'),'w');
    fwrite(fid_mod,mod,'float');
    fclose(fid_mod);
end

clear i j h x_grid y_grid z_grid;