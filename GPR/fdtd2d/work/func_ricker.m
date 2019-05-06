function [ t,a ] = func_ricker( fc,number,t0,dt )
    t=(0:number-1)*dt;
    a=(1-2*(pi*fc*(t-t0)).^2).*exp(-(pi*fc*(t-t0)).^2);
end