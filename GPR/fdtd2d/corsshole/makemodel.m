clear;
x=-0.5:0.025:10.5;z=0:0.025:11;
c=0.299792458;
background=9;
anomaly=5;
ep=background.*ones(length(z),length(x));
x01=2.75;x02=7.25;
y01=5.25;y02=5.75;
a=1.25;b=1;
for i=1:length(z)
    for j=1:length(x)
        if ((((x(j)-x01)/a)^2+((z(i)-y01)/b)^2)<=1 || (((x(j)-x02)/a)^2+((z(i)-y02)/b)^2)<=1)
            ep(i,j)=anomaly;
        end
    end
end
ep=ep';
mu=ones(size(ep));
sig=1e-3.*ones(size(ep));
imagesc(x,z,ep')
axis image
colormap(jet)