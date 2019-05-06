%creates a binary formation model for 3D seismic modeling
close all;
clear all;
clc;

%general parameters
writefiles=1;%1=yes else=no
modname='../par/model/test';
plotmodels=1;%1=yes else=no
grid_spacing=0.1;

%dimensions of base model
xmax_base=400;%gridpoints
ymax_base=100; %(vertical)

%change in material properties, i.e. fault zone
x_start_fault=1; % begin of new material, at fault zone in grid points
y_start_fault=50; % 

%velocity, density base model
vp_base=400; %m/s
vs_base=200;
rho_base=1700; %kg/m3

%velocity, density fault zone
vp_fault=800.0; %m/s
vs_fault=400.0;
rho_fault=1900.0; %kg/m3

%----------------VP model
%creating model
mod_vp=(zeros(xmax_base,ymax_base));
mod_vs=(zeros(xmax_base,ymax_base));
mod_rho=(zeros(xmax_base,ymax_base));
%filling in base parmaters
mod_vp(:,:)=vp_base;
mod_vs(:,:)=vs_base;
mod_rho(:,:)=rho_base;
%filling in fault parameters (i.e. change in material properties)
mod_vp(x_start_fault:end,y_start_fault:end)=vp_fault;
mod_vs(x_start_fault:end,y_start_fault:end)=vs_fault;
mod_rho(x_start_fault:end,y_start_fault:end)=rho_fault;

%plot models
if plotmodels==1
    
   h1=figure(1);
  
   subplot(2,2,1);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod_vp(:,:))');
   shading flat;
   colorbar
   title('vp-velocity Model (m/s)');
   xlabel('x in m');
   ylabel('y in m ');set(gca, 'YDir', 'reverse')
   axis equal tight;   
   
   subplot(2,2,2);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod_vs(:,:))');
   shading flat;
   colorbar
   title('vs-velocity Model (m/s)');
   xlabel('x in m');
   ylabel('y in m ');set(gca, 'YDir', 'reverse')
   axis equal tight;  
   
   subplot(2,2,3);
   pcolor(grid_spacing:grid_spacing:xmax_base*grid_spacing,grid_spacing:grid_spacing:ymax_base*grid_spacing,squeeze(mod_rho(:,:))');
   shading flat;
   colorbar
   title('Density Model (kg/m^3)');
   xlabel('x in m');
   ylabel('y in m ');set(gca, 'YDir', 'reverse')
   axis equal tight;  

end

%write models to disk
if writefiles==1
    mod_vp=permute(mod_vp,[2 1 3]);
    fid_mod_vp=fopen(strcat(modname,'.vp'),'w');
    fwrite(fid_mod_vp,mod_vp,'float');
    fclose(fid_mod_vp);
    
    mod_vs=permute(mod_vs,[2 1 3]);
    fid_mod_vs=fopen(strcat(modname,'.vs'),'w');
    fwrite(fid_mod_vs,mod_vs,'float');
    fclose(fid_mod_vs);
    
    mod_rho=permute(mod_rho,[2 1 3]);
    fid_mod_rho=fopen(strcat(modname,'.rho'),'w');
    fwrite(fid_mod_rho,mod_rho,'float');
    fclose(fid_mod_rho);
end

% clear mod;

