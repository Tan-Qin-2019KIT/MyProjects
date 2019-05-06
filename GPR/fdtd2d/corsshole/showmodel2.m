clear
load model2
load data2
% % model
% figure
% axes('position',[0.08,.15,.75,.8])
% imagesc(x,z,ep')
% set(gca,'xtick',[0:1:5]);
% xlabel('x (m)')
% z=ylabel('z (m)');
% set(z,'position',[-0.8 6])
% colormap(fliplr(jet));
% axis image
% axis([0 5 0 12])
% hold on
% plot(srcx+0.025,srcz,'ow','markersize',1);
% plot(recx+0.025,recz,'xw','markersize',3);
% % colorbar
% axes('position',[0.68,.25,.1,.5])
% axis off
% colormap(jet)
% h=colorbar;
% set(gca,'Clim',[22 28]);
% set(h,'ytick',[22:3:28]);
% set(h,'FontSize',10,'FontName','Times New Roman')
% set(get(h,'Title'),'string','\it\epsilon_r\rm [-]','FontSize',10,'FontName','Times New Roman');
% set(gcf, 'position', [0 0 300 300.*1.6].*0.65./1.2188);

% % data
figure
axes('position',[0.23,.15,.75,.8])
SS=22;
data=reshape(gather(:,:,SS),1601,45)';
data=satu2(data,0.5);
imagesc(t(1:800).*1e9,1:45,data(:,1:800))
xlabel('{\itt} (ns)')
z=ylabel('Trace No.');
set(z,'position',[-25 23.132743])
% set(gca,'ytick',1:2:21)
hold on
plot(d_obs((SS-1)*45+1:SS*45),1:45,'r')
colormap(gray)
set(gcf, 'position', [0 0 260 300.*1.6].*0.65./1.2188);