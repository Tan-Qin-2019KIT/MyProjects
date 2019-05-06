function [ t,a ] = func_ricker( fc,t,t0,dt )
    
    a=(1-2*(pi*fc*(t-t0)).^2).*exp(-(pi*fc*(t-t0)).^2);
end