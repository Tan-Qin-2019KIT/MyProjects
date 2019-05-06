                                                                                                                                                                                                                                                                                                                                                        %TM-UWB radar self-correlation receiver detection charateristic using
%autocorrelation receiver :Montel Carlo simulation
%-------------------------------------CFAR----------------------
clear
clc
%-------second receiver:differential averaged auto-correlation receiver(WADTR)---
%--------------------------build target model------------------
%-----------------two bodyes targets---one 1m second 3m range gate 6m(40ns)-----------------------
%------dispersive central 1.55m,1.66m,1.75m,3m,3.08m,3.13m----------------
%---------------------second-order gaussian pulse (fc,tm,tau)---------------

tau=0.5e-9;   tm=1e-9;%1ns

fc=50e9;%50GHz
 
dt=1/fc;

over=floor(tm/dt);

e=mod(over,2);   

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

%----------Analyzing UWB signal bandwidth, and get range resolution-------------

[Ess,f_high,f_low,BW] = cp0702_bandwidth(w0,dt,-3);

%-------------generate one dimensional high resolution range profile of single target-------

%Range=[1.55,1.66,1.75,3,3.08,3.13];   gain=[1,0.85,0.6,0.5,0.3,0.2];

Range=[2.55,2.66,2.75];               gain=[1,0.85,0.6]; 

%----------------------BW=1.5GHz,range resolution =10cm,-------------------

C=3*10^8;      Range_resolution_Cell = C/(2*BW); 
  
Range_resolution_Cell_samples = floor( 2 * Range_resolution_Cell / C / dt) ;
 
Range_delay = 2 * Range/ C;

Range_position_sample = floor( Range_delay / dt);

Range_Gate_Time = 60e-9;         

Range_Gate_samples = floor(Range_Gate_Time / dt); 

Range_Cell_Number=floor(Range_Gate_samples / Range_resolution_Cell_samples);

gate_seq=zeros(1,Range_Gate_samples);

for i=1:length(Range_delay)    

    k=Range_position_sample(i);
    
 if e   %sample=odd number
   
     gate_seq(k-kbk:k+kbk)=gate_seq(k-kbk+1:k+kbk+1) + w0 * gain(i);
     
else %sample=even number
    
    gate_seq(k-kbk+1:k+kbk)=gate_seq(k-kbk+1:k+kbk) + w0 * gain(i);
    
 end
 
end
%----------------signal ------------------

HRRP = gate_seq;

%-----normalization HRRP--------------------
HRRP = gate_seq  ./  max(HRRP);   % figure(2), plot(HRRP)

%-------------introduce additive white Gaussian noise over signal--------
ebno = [-20:2:10];                     %  SNR=E/NO

Target_Length_sample =( Range_position_sample(length(Range_delay))-Range_position_sample(1)) /2;

% -----average range gate energy of receive signal----------------

E_Tw = sum(HRRP.^2) / (Target_Length_sample + over) ; 

%------------------------averge noise power----
EbNo = 10.^(ebno./10);

N0 = E_Tw ./ EbNo;                

nstdv = sqrt(N0./2);

%------Given Pf, Generate decision detection threshold---------------------

Number_Mont_carlo=1000;

Acculating_Pulse_Number = 1;    %acculating the number of echo 

%----------------------------------------------------------------------------

Window_Cell_Number= 3 ;         %Window Number

%Window_Cell_Number = Range_Cell_Number; %ED

Moving_Window_Cell_Number = 1;  % Moving Window Cell

Window_Number = floor ( 1+ ( Range_Cell_Number - Window_Cell_Number ) / Moving_Window_Cell_Number);

%---------------------------------------------------------------------

CFAR_noise_output=zeros(length(EbNo) ,Number_Mont_carlo);

%-------according to SNR, generating the noise-------------------
for j=1:length(EbNo)    
  
 for kk=1:Number_Mont_carlo
     
 %-------------------------------initial Output of  pulses  accemulated
    
     Output_Window_Noise=zeros(1, Window_Number);
     
     %-------------------Generate several Noise signal----------------------
   
   for nn = 1: Acculating_Pulse_Number
    
       Receive_Noise  =  nstdv(j) .*  randn(1 , Range_Gate_samples) ;  
       
       Receive_Module_product_noise = Receive_Noise .*  Receive_Noise / Range_Gate_samples;       
            
  % ----------------Generate Range Profile----------------------
 
 for  mmm = 1 : Range_Cell_Number
 
    Range_Profile(mmm) =sum( Receive_Module_product_noise( (mmm-1) * Range_resolution_Cell_samples...
                                                     + 1 : mmm * Range_resolution_Cell_samples));
 end
 
 %---------------------Window Cell Processing------------------------------------   
 
for i=1 : Window_Number
   
%---------------------------------Sliding windows detection-----------   
 
Output_Window_Noise(i) = Output_Window_Noise(i) + sum( Range_Profile((i-1) * Moving_Window_Cell_Number .....
                                + 1 : (i-1) * Moving_Window_Cell_Number + Window_Cell_Number ));  
 
end
 
end %%% end for Average_Echo_Number
 
%--------------------Find Maximum Value of Window Noise----------------

CFAR_noise_output(j,kk) = max(Output_Window_Noise);

%---------------second: Select sevral larger Value of Output_Window_Signal-----

alfa=0.9; temp_Th = alfa * max(Output_Window_Noise);   temp2=0;

for ii=1:Window_Number
    
    if Output_Window_Noise(ii) >= temp_Th
        
    temp2 = temp2 + Output_Window_Noise(ii);
    
    end
    
end

 %CFAR_noise_output(j,kk)=temp2;
 
 %--------------------Three: -------------------------------------

%CFAR_noise_output(j,kk) = sum(Output_Window_Noise.^2);

%---------------------------------------------------------------

end  %end number_Mont_carlo

end %end for EbNo
%-------------------Pf for each SNR, get threshold------------
th_l1=zeros(1,length(EbNo));

for jj = 1 : length(EbNo)

xx = sort( CFAR_noise_output(jj,:) );

yy = length( CFAR_noise_output(jj,:) );

 %---------------Pf=0.01,----------------------

zz = floor(Number_Mont_carlo/100);     

th_l1(jj) =( xx( yy - zz + 1) + xx( yy - zz )  )/2 ;

end %end for

%--------------------------------------------------------------------------
%------Given Pf, Generate decision detection threshold: th_l, and calculate
%pd-----------------------------------------------------------------------
%-------according to SNR, generating the noise-------------------

CFAR_Signal_output=zeros(length(EbNo), Number_Mont_carlo);

for j = 1 :length(EbNo)    
  
  for kk = 1 : Number_Mont_carlo
     
   %--------------initial Output of  N pulses  accemulated------------------
    
     Output_Window_Signal =  zeros(1, Window_Number); 
       
   %-------------------Generate several signal + Noise signal----------------------
   
   for nn = 1: Acculating_Pulse_Number
    
       Receive_Signal_Noise  = HRRP + nstdv(j) .*  randn(1 , Range_Gate_samples) ;  
       
       Receive_Module_product_Signal = Receive_Signal_Noise .*   Receive_Signal_Noise / Range_Gate_samples;
     
  % ----------------Generate Range Profile----------------------
 
 for  mmm = 1 : Range_Cell_Number
 
    Range_Signal_Profile(mmm) =sum( Receive_Module_product_Signal( (mmm-1) * Range_resolution_Cell_samples......
                                                     + 1 : mmm * Range_resolution_Cell_samples));
 end

 %---------------------accemuating again, Window accemulated Processing of N Pulse------------------------------------   
 
for i=1 : Window_Number
   
 Output_Window_Signal(i) =  Output_Window_Signal(i) + sum( Range_Signal_Profile((i-1) * Moving_Window_Cell_Number .....
                                + 1 : (i-1) * Moving_Window_Cell_Number + Window_Cell_Number ));  
 end
 
end %%% end for Acculating_Pulse_Number
       
%--------------------first: Find Maximum Value of Window Noise----------------

CFAR_Signal_output(j,kk) = max(Output_Window_Signal);

%---------------second: Select sevral larger Value of Output_Window_Signal-----

temp1 = alfa * max(Output_Window_Signal); temp3=0;

for ii=1:Window_Number
    
    if Output_Window_Signal(ii) >= temp1
        
    temp3 = temp3 + Output_Window_Signal(ii);
    
    end
    
end

%CFAR_Signal_output(j,kk) = temp3;

%--------------------Three: -------------------------------------

%CFAR_Signal_output(j,kk) = sum(Output_Window_Signal.^2);

%--------------------------------------------------
end  %end number_Mont_carlo

end %end for EbNo
 

%--------------target decision detection-----

Pd1=zeros(1,length(EbNo));

for jj=1:length(EbNo)
   
    target_Pf=0;
   
for kk=1:Number_Mont_carlo
    
        if  CFAR_Signal_output(jj,kk) > th_l1(jj) ,    target_Pf = target_Pf + 1; end
        
    end %end for decision detection
%-------------------Pd for each SNR----------------------------------------

  Pd1(jj) = target_Pf / Number_Mont_carlo;
end
%%%-----------------plot ---------------------------------------------------
figure(3),

PT=plot(ebno,Pd1,'ro-');hold on,

set(PT,'LineWidth',[1]);set(gca,'yscale','linear');set(gca,'ytick',[0 0.2 0.4 0.6 0.8 1]);

X=xlabel('SCR /dB'); set(X,'FontSize',12);

Y=ylabel('Detection Probability (Pd)');  set(Y,'FontSize',12);

grid