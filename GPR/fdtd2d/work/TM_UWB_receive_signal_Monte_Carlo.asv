                                                                                                                                                                                                                                                                                                                                                        %TM-UWB radar self-correlation receiver detection charateristic using
%Montel Carlo simulation 
%--------receiver signal= anntenna couping siganl + front wall reflecting signal+ target refraction signal----
clear
clc

%--------------------------build target model------------------
%-----------------two bodyes targets---one 1m second 3m range gate 6m(40ns)-----------------------
%------dispersive central 1.55m,1.66m,1.75m,3m,3.08m,3.13m----------------

%---------------------second-order gaussian pulse (fc,tm,tau)---------------
tau=0.5e-9;

tm=1e-9;%---------pulse_width = 1ns

fc=50e9;%-----------sampling rate = 50GHz         

dt=1/fc;%-----------sampling interval=20ps

%--------------------------------------------------------------------------
over=floor(tm/dt); e=mod(over,2); kbk=floor(over/2); tmp=linspace(dt,tm/2,kbk);

s = (1-4*pi* (tmp./tau).^2 ) .* exp(-2*pi* (tmp./tau).^2);
if e  %sample=odd number
    for k=1:length(s)
        y(kbk+1)=1;         y(kbk+1+k)=s(k);        y(kbk+1-k)=s(k);
    end
else %sample=even number
    for k=1:length(s)
        y(kbk+k)=s(k);        y(kbk+1-k)=s(k);
    end
end

w0=y; 
%figure(2); plot(w0),xlabel('ns');ylabel('second order Gaussian pulse')

%-----------------generate one dimensional high resolution range profile-----
%------dispersive central 1.55m,1.66m,1.75m,3m,3.08m,3.13m, and gains individually----------------

%-------------according to SNR, generating the additive white Gaussian noise over signal----
ebno=[-20];            %SNR=E/NO (dB)

% the target velocity of V=0.274m/s, Tupdate=0.296s
 target_v=0.274 ;
 
 Time_update_scans=0.296;
 
 scans=100;
  
 Range_initional=0.5; C=3*10^8;

Tw=60*10^(-9);   %gate duration=60ns, Observe range=9m

gate_samples=floor(Tw/dt);

receive_signal=zeros( scans ,gate_samples);

for iii = 1 : scans  
  
    Range_target = Range_initional +  target_v * Time_update_scans * (iii - 1);
    
    Range_delay = 2 * Range_target/C;

    Range_position_sample = floor(Range_delay/dt);
    
    gate_seq=zeros(1,gate_samples);
    
    for i=1:length(Range_delay)    
    
    k= Range_position_sample(i);
    
      if e   %sample=odd number
     
            gate_seq(k-kbk:k+kbk)=gate_seq(k-kbk+1:k+kbk+1) + w0;
   
          else %sample=even number
    
           gate_seq(k-kbk+1:k+kbk)=gate_seq(k-kbk+1:k+kbk) + w0;
    
         end  %end if
   end  %end for range_delay

HRRP= gate_seq;

E_Tw=sum(abs(HRRP).^2)/gate_samples; % average range gate energy

EbNo=10.^(ebno./10);

N0=E_Tw ./ EbNo;                     % averge noise power

nstdv=sqrt(N0);
    
cumlative_number_of_pulse=100;

%-------according to SNR, generating the noise-------------------

temp = zeros(1,gate_samples);

for j=1:cumlative_number_of_pulse

 temp    =  temp  +  HRRP +  nstdv .* randn(1,gate_samples) ;

  end  %end for j
  
 receive_signal(iii,:) = temp ;  
 
end %end for iii scans
 % figure(2),plot(receive_signal(3,:));

%figure(1),disp_time=[1: gate_samples]*dt * 1e+9;
%plot(disp_time,receive_signal ,'r-'),set(gca,'yscale','linear');  ylabel('Radar receiver signal ');  xlabel('Time (ns)'); title('Through-wall ')
 
  yaxis=[ 1 : gate_samples ] * dt * 1e+9;
  
  receive_signal_db=20.*log10(abs(receive_signal));
  
  drpmax=max(max(receive_signal_db));
  
  figure(1),  
 imagesc( (1:scans) * Time_update_scans ,yaxis,receive_signal_db'-drpmax,[-30 0]);  hold on
 colormap('gray');
  ylabel( 'sampling time (ns)/ each scan  ');
  xlabel('scan (s) , each scan Tupdate=0.296s , v=0.0274m/s');
 title(' range profile ')
  colorbar('vert');  
  
%---------------------Get move target receive signal--------------------
  
for iii =20 : scans  
    
     move_target_siganl(iii,:) = receive_signal(iii,:) - receive_signal (iii-19,:) + 0.0001;
     
 end
 
  yaxis=[ 1 : gate_samples ] * dt * 1e+9;
  
  move_target_siganl_db=20.*log10(abs(move_target_siganl));
  
  drpmax=max(max(move_target_siganl));
  
  figure(2),  
 imagesc( (1:scans-20) * Time_update_scans ,yaxis,move_target_siganl_db'-drpmax,[-30 0]); 
  colormap('gray');
  ylabel( 'sampling time (ns)/ each scan  ');
  xlabel('scan (s) , each scan Tupdate=0.296s , v=0.0274m/s');
 title(' range profile ')
  colorbar('vert');  
 
  
  