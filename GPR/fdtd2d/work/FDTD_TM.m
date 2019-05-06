% TM wave simulation based on FDTD compiling by Tan Qin in March 2019
% the differential order is O(2,2)
clear all; close all;
% basic parameters
epsilon_0=8.8541878176e-12; mu_0=4*pi*1e-7;

% model space and grid size
dx=0.2; dy=dx;
nx=81; ny=97;
x=(1:nx)*dx; y=(1:ny)*dy;
% parameters distribution
epsilon_r=1;mu_r=1;
epsilon=epsilon_r*epsilon_0; 
mu=mu_r*mu_0; 

% stability criterion
Z=sqrt(mu./epsilon);
c=sqrt(1./(mu.*epsilon));
dtmax = finddt(1,1,dx,dy);
dt=floor(dtmax*1e9*1e2)./(1e9*1e2);
dtau=1./2*dx;
% dt=dtau./c;
nt=2.^7; t=(0:nt-1)*dt;
alpha=8;

% source wavelet
t0=15e-9; fc=60e6;
% w=zeros(nt,1);
% theta=x(50)-10*alpha+c*t;
% for it=1:nt      
%     if theta(it)>=0 && theta(it)<=8*alpha
%         w(it,1)=sin(theta(it)*pi./(8*alpha));
%     end
% end
% [ t,w ] = func_ricker( fc,nt,t0,dt );
w = blackharrispulse(fc,t);

Ez=zeros(nx,ny); Hx=Ez; Hy=Ez;

x0=50*dx; y0=30*dy;
ix0=x0./dx;jy0=y0./dy;
for it=1:nt
%     there is no absorbing boundary condition in this code
    
    ix=2:nx-1;    jy=2:ny-1; % utility of the matrix operation
    Hx(ix,jy)=Hx(ix,jy)-1/Z*dtau/dy*(Ez(ix,jy+1)-Ez(ix,jy));
    Hy(ix,jy)=Hy(ix,jy)+1/Z*dtau/dx*(Ez(ix+1,jy)-Ez(ix,jy));
    Ez(ix,jy)=Ez(ix,jy)+Z*dtau/dx*(Hy(ix,jy)-Hy(ix-1,jy))...
        -Z*dtau/dy*(Hx(ix,jy)-Hx(ix,jy-1));

    Ez(ix0,jy0)=Ez(ix0,jy0)+w(it);
    Hy(ix0,jy0)=1./Z*Ez(ix0,jy0);
%     Hx(ix0,jy0)=Hy(ix0,jy0);
    if mod(it,10)==0
        figure (1)
        snap=imagesc(x,y,Ez'/max(max(abs(Ez))));grid on;hold on;
%         plot(x0,y0,'w','markersize',50);
        set(gca,'YDir','normal'); axis equal;
        title(['TM-snapshot (' num2str(t(it)*1e9) 'ns)']);
        xlabel('x (m)');
        ylabel('y (m)');
        colorbar; colormap(jet);
        set(gca, 'CLim', [-1 1]);
        pause(0.1)
    end
end
