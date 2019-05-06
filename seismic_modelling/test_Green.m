%this program is for calculating the Green function under spherical tissue
%model.

clear;clc;

ua=0.02;us=15;g=0.94; %��������
a=10; %�뾶
r1=5; %��Դ
r2=a^2/r1; %���Դ
D=1/(3*(ua+us*(1-g))); %����ϵ��
k=sqrt(-ua*D); %k^2=-Q
r = linspace(0,10,50); %x=r
th = linspace(0,2*pi,120);
[th,r] = meshgrid(th,r);

c = -a/r1*exp(i*k*(1-a/r1)*sqrt(a^2+r1^2-2*a*r1.*cos(th)));
G0 = (1/(4*pi) * exp(i*k*sqrt(abs(r.^2 + r1^2-2*r*r1.*cos(th)))))./sqrt(abs((r.^2+r1^2)-(2*r.*r1).*cos(th)));
G1 = ((1/(4*pi) * exp(i*k*sqrt(abs(r.^2 + r2^2-2*r.*r2.*cos(th))))).*c)./sqrt(abs(r.^2+r2^2-(2*r.*r2).*cos(th)));

G = (G0 + G1)/D; %��ͨ��

%�������»�����ͼ
[X,Y] = pol2cart(th,r);
surfc(X,Y,G);

axis([-10,10,-10,10]);
axis equal;
colormap(gray);
shading interp;
%view([0,90]); 