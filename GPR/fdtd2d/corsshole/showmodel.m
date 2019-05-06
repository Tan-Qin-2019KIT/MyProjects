clear
load model1
load data1
% model
figure
axes('position',[0.15,.2,.6646,.72])
imagesc(x,z,ep')
xlabel('x (m)')
z=ylabel('z (m)');
set(z,'position',[-1.7 5.5220941883767525 ])
colormap(fliplr(jet));
axis image
hold on
plot(srcx,srcz,'ow','markersize',3);
plot(recx,recz,'xw','markersize',3);
% colorbar
axes('position',[0.775,.25,.1,.5])
axis off
colormap(fliplr(jet))
h=colorbar;
set(gca,'Clim',[5 9]);
set(h,'ytick',[5:2:9]);
set(h,'FontSize',10,'FontName','Times New Roman')
set(get(h,'Title'),'string','\it\epsilon_r\rm [-]','FontSize',10,'FontName','Times New Roman');
set(gcf, 'position', [0 0 325 300].*0.65);
% data
figure
axes('position',[0.28,.2,.7,.75])
data=reshape(gather(:,:,15),1601,21)';
data=satu2(data,0.1);
imagesc(t.*1e9,1:21,data)
xlabel('{\itt} (ns)')
ylabel('Trace No.')
% set(gca,'ytick',1:2:21)
hold on
plot(d_obs(15:21:end),1:21,'r')
colormap(gray)
set(gcf, 'position', [0 0 225 300].*0.65);