clear all;close all;
fileout = 'model3.out';
% fileout = 'pec.sca';
[Header,Fields]=gprmax(fileout);
% N=1:Header.NSteps;               %�ƶ�����
% Position=Header.dx*Header.tx+(N-1)*(Header.dx*Header.TxStepX); %����ÿ������λ��

Data(:,:)=Fields.ez(1:2:end,20,:);   %ת��������?
dt=(Fields.t(2)-Fields.t(1))*1e9*2;
[m,n]=size(Data);
% imagesc(1:n,(0:m-1).*dt,Data);%��ͼ
% colorbar
% xlabel('Trace No.');
% ylabel('t (ns)');
% hold on
mcmc=pickfisrt(Data,51);
fa=mcmc.*dt;
plot(1:n,fa,'r')
figure
wiggledisplay(-Data,1:n,(0:m-1).*dt,'vararea',n,2);
ylabel('\itt\rm (ns)');
set(gca,'YTick',0:40:(m-1)*dt);
xlabel('Trace No.');
hold on
plot(1:n,fa,'r')
title('Model 3')
set(gcf, 'position', [0 0 600 800]);
set(gcf, 'position',[810.0000  209.0000  283.0000  458.0000]);
% set(gca,'Position',[.09 .06 .9 .92]);
set(gca,'Position',[0.1600    0.1000    0.8000    0.8400]);