%TM-UWB radar autocorrelation receiver detection charateristic
%-------------------------------------CFAR----------------------
clear
clc
%-------------------optimal correlation detection receiver--------------
% Q(x)=erfc(x)
pfa=10^(-2);   beta=erfcinv(pfa);

Tw=60*10^(-9);%range gate=60ns,observed range=60ns*15cm=9m

B=2*10^9; %Bandwidth=2GHz

df=B*Tw;%degreeth of freedom

SNR=-20:1:10; gamma_g=10 .^(SNR/ 10);

N=10; %the number of pulse cummlative

temp=sqrt(2*N*gamma_g );

pd=erfc(beta-temp);

[m,n]=size(pd);

for i=1:n if pd(i) > 1 ,    pd(i)=1; end ; end

figure(1);semilogy(SNR,pd,'ko-');hold on

%-------first receiver:differential self-correlation receiver(DTR)---
% Q(x)=erfc(x)

%pfa=10^(-4);   beta=erfcinv(pfa);

Tw=60e-9;%range gate=60ns,observed range=60ns*15cm=9m

B=2e9;   %Bandwidth=2GHz

df=B*Tw; %degreeth of freedom

SNR=[-20:1:10]; %dB

gamma_g=10 .^(SNR/ 10);

N=10; %the number of pulse cummlative

temp1 = beta - gamma_g .* sqrt(2 * N / df);

temp2 = sqrt( 1 + 2 * gamma_g ./ df );

pd2=erfc(temp1 ./ temp2);

[m,n]=size(pd2);

for i=1:n
    
    if pd2(i) > 1 ,    pd2(i)=1; end ;
    
end

figure(1);plot(SNR,pd2,'k*-');hold on;

%-------second receiver:differential averaged self-correlation receiver(WADTR)---
%--------------------------build target model------------------
%-----------------two bodyes targets---one 1m second 3m range gate 6m(40ns)-----------------------
%------dispersive central 1.55m,1.66m,1.75m,3m,3.08m,3.13m----------------

%---------------------second-order gaussian pulse (fc,tm,tau)---------------
tau=0.5e-9;   tm=1e-9;%1ns

fc=50e9;%50GHz

dt=1/fc; over=floor(tm/dt);   e=mod(over,2);

kbk=floor(over/2);

tmp=linspace(dt,tm/2,kbk);

s=(1-4*pi* (tmp./tau).^2 ) .* exp(-2*pi* (tmp./tau).^2);

if e  %sample=odd number
    for k=1:length(s)
        
        y(kbk+1)=1;        y(kbk+1+k)=s(k);        y(kbk+1-k)=s(k);
        
    end
else %sample=even number
    for k=1:length(s)
        
        y(kbk+k)=s(k);        y(kbk+1-k)=s(k);
        
    end
end

E=sum((y.^2).*dt );%pulse energy
w0=y./ sqrt(E); % energy normalization
Td=60e-9; %period of pulse =TD=100ns
power=0; % average transmitted power (dBm watt)
power1=(10^(power/10))/1000;%% average transmitted power (watt)
Ex=power1*Td; %energy per pulse
w0=w0.*sqrt(Ex);

%w0=y;figure(2); plot(w0),xlabel('ns');ylabel('second order Gaussian pulse')

%-----------------generate one dimensional high resolution range profile -
Range=[1.55,1.66,1.75,3,3.08,3.13];   gain=[1,0.85,0.6,1,0.85,0.6];

%Range=[1.55,1.66,1.75];gain=[1,0.85,0.6,];


C=3*10^8;   Range_delay=2*Range/C;

Range_position_sample=floor(Range_delay/dt);

Tw=60e-9;  %gate duration=60ns

gate_samples=floor(Tw/dt);

gate_seq=zeros(1,gate_samples);

for i=1:length(Range_delay)    

    k=Range_position_sample(i);
    
 if e   %sample=odd number
       gate_seq(k-kbk:k+kbk)=gate_seq(k-kbk+1:k+kbk+1) + w0 * gain(i);
else %sample=even number
       gate_seq(k-kbk+1:k+kbk)=gate_seq(k-kbk+1:k+kbk) + w0 * gain(i);
 end
    
end

%----------------signal ------------------

HRRP= gate_seq;

%-------------introduce additive white Gaussian noise over signal

ebno=[-20:1:10];                     %SNR=E/NO

E_Tw = sum( abs( HRRP ).^2) /gate_samples; % average range gate energy

EbNo=10.^(ebno./10);

N0=E_Tw ./ EbNo; % averge noise power

nstdv=sqrt(N0./2); 

%------------------generate noise signal in order to check correct of noise energy-----------------------------

for j=1:length(EbNo)
    
     noise(j,:)= nstdv(j).*randn(1,gate_samples) ;  N00(j)= sum(abs(noise(j,:)).^2)/gate_samples;   N11(j)=std( noise(j,:)).^2;
   
end
%figure(2),plot(HRRP_add_noise(1,:));

%----------generate reference module signal---------average advanced M group receiver signal of range gate 
M=3;

HRRP_add_noise1=zeros(length(EbNo),gate_samples);

for m=1:M
    
    for j=1:length(EbNo)
        
     HRRP_add_noise1(j,:)=HRRP_add_noise1(j,:) + HRRP + nstdv(j).*randn(1,gate_samples) ;
    
    end
 
end

module_noise=HRRP_add_noise1./M;

%-----------------------differential crosscorrelation----------------
Tm=20e-9;        % 10ns----sub-interval width duration ----------------------

Tm_sample=floor(Tm/dt);

Number_sub_interval=floor(gate_samples/Tm_sample);

move_interval=10e-9;%start point interval = 1ns

move_sample=floor(move_interval/dt);

Number_sub_interval1=floor((gate_samples - Tm_sample) / move_sample);

if ( (gate_samples - Number_sub_interval1 * move_sample )== Tm_sample)
    
    Number_sub_interval=Number_sub_interval1+1;
    
else
    
    Number_sub_interval=Number_sub_interval1;
    
end
Number_sub_interval
Number_sub_interval1

%---------------------accumulated N receive_siganl and siganl+niose_energy ---------------------
  
N=10;     %--number of pulse cummlative------------- 

receive_HRRP_add_noise1=zeros(length(EbNo),gate_samples);

for kk=1:N
%----------receiver signal---------
for j=1:length(EbNo)

    receive_HRRP_add_noise1(j,:) = receive_HRRP_add_noise1(j,:) + HRRP + nstdv(j).*randn(1,gate_samples) ;
    
end

end

%-----------------------------siganl+niose_energy of sub_interval---------
for j=1:length(EbNo)
    
for i=1:Number_sub_interval
    
    %------------------------move windows detection-----------
    ll = module_noise  (j, (i-1) * move_sample + 1 : (i-1) * move_sample +  Tm_sample);
    
    mm = receive_HRRP_add_noise1(j, (i-1) * move_sample + 1 : (i-1) * move_sample +  Tm_sample);
    
    energy_sub_interval(j,i)= (ll * mm') /gate_samples;

end
end
%figure(4),plot(energy_sub_interval(10,:));
 
%-------------------if echo signal have target signal, then the energy integration_sub_interval of signal---------------------
for j=1:length(EbNo)
    
for i=1:Number_sub_interval

     %------------------------move windows detection-----------
    mmm = HRRP((i-1) * move_sample + 1 : (i-1) * move_sample + Tm_sample );
    
    signal_energy_sub_interval(j,i)= N * (mmm * mmm') / gate_samples;
    
end

end

%----------------------calculate detection probability----------------
%pfa=10^(-4);   beta=erfcinv(pfa);

B=2*10^9; %Bandwidth=2GHz

df=B*Tm;%degreeth of freedom

%--------MRC(maximum ratio combination)---------------

weight_MRC1=(energy_sub_interval);   

%weight_MRC2=sqrt(abs(energy_sub_interval));    %weight_MRC1=weight_MRC2;

%--------EGC(equal gain combination)---------------

%weight_EGC=ones(length(EbNo),Number_sub_interval);  weight_MRC1=weight_EGC;

%---------------SD(select diversity)------------
th=0.03;  
weight_SD=zeros(length(EbNo),Number_sub_interval);
for j=1:length(EbNo)
   for i=1:Number_sub_interval
       if energy_sub_interval(j,i)>th
        weight_SD(j,i)=1;
    end
   end
end
%weight_MRC1=weight_SD;

%-----------calculate detection probability----------------------------------------------

for j=1:length(EbNo)
    
temp1 = beta * sqrt( 1/(2*M) * N * df * N0(j)^2 * weight_MRC1(j,:) * weight_MRC1(j,:)' ) - ....
           weight_MRC1(j,:) * signal_energy_sub_interval(j,:)';

temp2 = sqrt( ( 1 + 1/M ) *  N0(j) / 2 * ( (weight_MRC1(j,:).^2 ) * signal_energy_sub_interval(j,:)') ....
        + 1/( 2 * M) * N * df * N0(j)^2 * (weight_MRC1(j,:) * weight_MRC1(j,:)' ) );

pd3(j)= erfc(temp1 /temp2);

if pd3(j) > 1 ,    pd3(j)=1; end

end ;
%%%-----------------plot ---------------------------
figure(1)

PT=plot(ebno,pd3,'r^-');hold on,

set(PT,'LineWidth',[1]);set(gca,'yscale','linear');set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1]);

X=xlabel('SCR (dB)'); set(X,'FontSize',12);

Y=ylabel('Detection Probability (Pd)');  set(Y,'FontSize',12);
grid
