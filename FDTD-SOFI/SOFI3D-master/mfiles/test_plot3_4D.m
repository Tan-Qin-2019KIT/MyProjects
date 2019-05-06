close all; clear all; clc

x=linspace(0,1,100);

y=linspace(0,1,100);

[X,Y]=meshgrid(x,y);

for i=0:1

    Z=linspace(i,i,100);

    plot3(X,Z,Y,'r');hold on; 

    plot3(Y,Z,X,'r');hold on;

    plot3(Z,X,Y,'y');hold on; 

    plot3(Z,Y,X,'y');hold on; 

    plot3(X,Y,Z,'b');hold on; 

    plot3(Y,X,Z,'b');hold on; 

end

xlabel('X','fontsize',20);ylabel('Y','fontsize',20);

zlabel('Z','fontsize',20)