% A MATLAB implementation of the 2D FDTD theory for scattering off a PEC cylinder,
% as described in Chapter 3, Section 3.2, 
% D.B.Davidson, "Computational Electromagnetics for RF and Microwave Engineering", CUP 2005.
% Usage: fdtd_2D_demo_v_1_2

% Options are provided at run-time to include the  PEC cylinder or not; and
% also to refine the mesh in a simple fashion. Enter 1 for reasonable
% results, with reasonable execution time. 

% The first plot shows the position of the cylinder, and the
% scattered/total field boundary. Subsequent plots show the pulse
% propagating through the mesh, and scattering off the cylinder. This is
% repeated as a movie.

% Also shown is a plot in the time domain of the scattered field at coordinates (point1_x,point1_y).
% Data thus saved can be used to produce results similar to Fig. 3.13, using the
% "cyl_plot.m" script.

% Commented out sections of code indicate how systematic testing may be
% undertaken, as discussed in Section 3.2.6.

% Author: D.B.Davidson
% Originally written 22 Feb 2003, revised 2 Feb 2005.
% Incorporates corrections 09 Oct 2007. E and H are now dimensioned
% (slightly) differently. 

clear;
format compact;
% Selectively includes PEC cylinder, centered at (N_centre_x,N_centre_y)
cyl_present=input('Include PEC cylinder? (1 or true)')

% Physical constants

c = 2.997925e8;    % [m/s] Speed of light in vacuum
eps_0 = 8.854e-12; % [F/m] epsilon_0
mu_0 = 4*pi*1e-7;  % [H/m] mu_0
eta_0 = sqrt(mu_0/eps_0);
                   % [Ohm] Wave impedance of free space
% Set parameters for simulation

refine = input('Factor to refine mesh? 1 standard')    % A quick method to refine the mesh 
pulse_compress = 1 % A quick way to shorten the pulse

N_x = refine*400     % number of cells in x-direction
N_y = refine*200     % ditto y
M = refine*1024       % Number of time steps
L = round(N_x/2)        % scat/tot field boundary
delta_s = 0.005*2/refine % [m] spatial step. This generates a physically larger computational volume.
%delta_s = 0.005/refine % [m] spatial step
R = 1        % fraction of Courant limit. Must be <= 1
delta_t = R* delta_s/(c * sqrt(2)) 
              % [s] Time step size
sigma = 1.0e-10/pulse_compress
              % Controls spectral content of Gaussian derivative pulse - 
              % equals 1/omega_max
f_max = (1/sigma)/(2*pi) % Freq of largest significant spectral component, in Hz. Information only
m_offset = 4*sigma;  % Controls switch-on time
Peak = 1;      % Peak amplitude of E field
% Set parameters for PEC cylinder
radius = 0.03 % [m] radius of cylinder
N_centre_x = round(0.75*N_x)
N_centre_y = round(0.5*N_y)

% Check that the simulation specification is valid:

if ((N_centre_x-L)*delta_s <= radius)
   disp('Error in simulation data. Scattered/total field not entirely to the left of target')
   return
end

%------------------------------------------

% Set up material grid (free space to start)
C_Ex = ones(N_x+1,N_y+1)*delta_t/(eps_0*delta_s); 
C_Ey = ones(N_x+1,N_y+1)*delta_t/(eps_0*delta_s); 
D_Hz = ones(N_x,N_y)*delta_t/(mu_0*delta_s); 
% Now force the electric fields to zero inside (and on the surface of) the PEC
% Note that the indices of the centre are treated as per usual FDTD
% indices, i.e. the actual location is:
%x_c=(N_centre_x-1))*delta_s ; y_c=(N_centre_x-1))*delta_s
if cyl_present % Otherwise just leave it as free space
  for ii = 1:N_x+1
    for jj = 1:N_y+1
        if ( sqrt( ((ii-1/2-(N_centre_x-1))*delta_s)^2 +  ((jj-1-(N_centre_y-1))*delta_s)^2 ) <=  radius ) 
           C_Ex(ii,jj) = 0;
        end
        if ( sqrt( ((ii-1-(N_centre_x-1))*delta_s)^2 +  ((jj-1/2-(N_centre_y-1))*delta_s)^2 ) <=  radius ) 
           C_Ey(ii,jj) = 0;
        end
     end
  end
end 
%------------------------------------------

% Set up storage for time histories.
H_z_point1 = zeros(1,M);
E_y_point1 = zeros(1,M);
point1_x = N_x/4;
point1_y = N_y/2;
H_z_point2 = zeros(1,M);
E_y_point2 = zeros(1,M);
point2_x = round((N_x+L)/2);
point2_y = N_y/2;


% Produce a simple graphical output, showing the cylinder, scat/tot zone
% interface and the point at which the scattered field will be computed. 
mesh_pic=zeros(N_x+1,N_y+1);
for ii = 1:N_x+1
    for jj = 1:N_y+1
      if C_Ex(ii,jj) == 0
        mesh_pic(ii,jj) = N_x/2; % To get vertical scale in plot OK when plotted
      elseif ii==L
        mesh_pic(ii,jj) = N_x/4;
      elseif (ii==point1_x & jj==point1_y)
        mesh_pic(ii,jj) = N_x/2;
      end
    end 
end
mesh((1:1:N_y+1),(1:1:N_x+1),mesh_pic);
axis image;
title('Simulation region')
disp('Press any key to continue')
pause;


% First time step - Initialize values for H_z, E_x and E_y. See note below
% w.r.t matrix sizes.
H_z_nmin1 = zeros(N_x,N_y); %
E_x_nmin1 = zeros(N_x+1,N_y+1); %
E_y_nmin1 = zeros(N_x+1,N_y+1); %
% Pre-allocation
H_z_n = zeros(N_x,N_y); %
E_x_n = zeros(N_x+1,N_y+1); %
E_y_n = zeros(N_x+1,N_y+1); %
% Get CPU time
start_time=cputime;
%------------------------------------------

movie_count = 1;
movie_interval = refine*50,

report_time_interval = 50; 


% Time loop
for m = 2:M,
    
  % A note on dimensions. Due to the nature of the Yee cell, the E and H fields 
  % have slightly different dimensions. 
  % E_x: N_x+1 by N_y+1
  % E_y: N_x+1 by N_y+1
  % H_z: N_x   by N_y
  % Essentially, it is because the H field are in the interior of the cells, whereas the E fields are
  % on the exterior of the cells. 
  % This must be explicitly taken into account in the matrix sizes and the update equations.
    
  % Update H fields: 
  H_z_n(1:N_x,1:N_y) = H_z_nmin1(1:N_x,1:N_y) ...
      + D_Hz(1:N_x,1:N_y).*(  E_x_nmin1(1:N_x,2:N_y+1)   - E_x_nmin1(1:N_x,1:N_y) ...
                           + E_y_nmin1(1:N_x,1:N_y) - E_y_nmin1(2:N_x+1,1:N_y) ) ;
 % Drive a test line source - used to check basic operation
 % H_z_n(N_x/2,N_y/2) = gaussder((m-1)*delta_t,m_offset,sigma);

                       
  % Special update on scat/tot field boundary
  E_y_nmin1_inc = ones(1,N_y)*Peak*gaussder_norm((m-1)*delta_t - (L-1)*delta_s/c,m_offset,sigma) ;
  H_z_n(L,1:N_y) = H_z_nmin1(L,1:N_y) ...
      + D_Hz(L,1:N_y).*(  E_x_nmin1(L,2:N_y+1)     - E_x_nmin1(L,1:N_y) ...
                        + E_y_nmin1(L,1:N_y) + E_y_nmin1_inc(1:N_y) - E_y_nmin1(L+1,1:N_y)) ;
  % Update E fields: (note that latest H fields must be used!)
  E_x_n(1:N_x,2:N_y) = E_x_nmin1(1:N_x,2:N_y) ...
      + C_Ex(1:N_x,2:N_y).*(  H_z_n(1:N_x,2:N_y)  - H_z_n(1:N_x,1:N_y-1) ) ; % Correction in 1st index
  E_y_n(2:N_x,1:N_y) = E_y_nmin1(2:N_x,1:N_y) ...
      - C_Ey(2:N_x,1:N_y).*(  H_z_n(2:N_x,1:N_y)  - H_z_n(1:N_x-1,1:N_y) ) ; % Correction in 2nd index
  
  % Special update on scat/tot field boundary (only needed for Ey)
  H_z_n_inc = ones(1,N_y)*(Peak/eta_0)*gaussder_norm((m-1/2)*delta_t - (L-1/2)*delta_s/c,m_offset,sigma) ;
  E_y_n(L,1:N_y) = E_y_nmin1(L,1:N_y) ...
      - C_Ey(L,1:N_y).*(  H_z_n(L,1:N_y)  - H_z_n_inc(1:N_y) - H_z_n(L-1,1:N_y) ) ; % Correction in 2nd index
  
  
  % Impose ABC on sides - assumes free space cell on boundary.
  % Left/right boundaries:
  E_y_n(1,1:N_y+1)   = E_y_nmin1(1,1:N_y+1)*(1-c*delta_t/delta_s)   + c*delta_t/delta_s*E_y_nmin1(2,1:N_y+1); % Correction in 2nd index
  E_y_n(N_x+1,1:N_y+1) = E_y_nmin1(N_x+1,1:N_y+1)*(1-c*delta_t/delta_s) + c*delta_t/delta_s*E_y_nmin1(N_x,1:N_y+1); % Correction in 2nd index
  % Top/bottom boundaries:
  E_x_n(1:N_x+1,1)   = E_x_nmin1(1:N_x+1,1)*(1-c*delta_t/delta_s)   + c*delta_t/delta_s*E_x_nmin1(1:N_x+1,2); % Correction in 1st index
  E_x_n(1:N_x+1,N_y+1) = E_x_nmin1(1:N_x+1,N_y+1)*(1-c*delta_t/delta_s) + c*delta_t/delta_s*E_x_nmin1(1:N_x+1,N_y); % Correction in 1st index
  
  % Fix outer values of E_tangential as PEC:
  %E_y_n(1,:) = 0;
  %E_y_n(N_x,:) = 0;
  %E_x_n(:,1) = 0;
  %E_x_n(:,N_y) = 0;
  
  % Store data
  % Movie
  if mod(m,movie_interval) == 0 
    mesh(eta_0*H_z_n) % Normalize
    title(strcat('\eta_o H_z field at timestep ',num2str(m)))
    H_z_Movie(movie_count) = getframe;
    movie_count = movie_count +1;
  end
  % Time history
  H_z_point1(m) = H_z_n(point1_x,point1_y);
  H_z_point2(m) = H_z_n(point2_x,point2_y);
  E_y_point1(m) = E_y_n(point1_x,point1_y);
  E_y_point2(m) = E_y_n(point2_x,point2_y);

  % Update for next iteration
  H_z_nmin1 = H_z_n;
  E_x_nmin1 = E_x_n;
  E_y_nmin1 = E_y_n;
  
  % Output some indication of how far the code is:
  % disp('.')
  if(rem(m,report_time_interval)==0)
     m
  end 
end
% End of time stepping
%------------------------------------------
disp('CPU time in seconds for time stepping')
run_time = cputime-start_time % seconds
floprate = 11*N_x*N_y*(M-1)/run_time

% movie(H_z_Movie,1,4);

time=[0:M-1]*delta_t;
plot(time/1e-9,E_y_point1)
xlabel('Time [ns]')
ylabel('E_y [V/m]')


print -deps RFGaussDerSimABC
disp('Press any key to continue')
pause;

% Compute and store time histories of source pulse

for mm = 1:M
  t(mm) = delta_t*(mm-1); 
  GaussDerSource(mm) = gaussder_norm(t(mm),m_offset,sigma);
end

delta_f = 1/((M-1)*delta_t);
freq=delta_f*(0:M-1);

clf;
run_movie = input('Run movie? (y/n)','s')    % Movie not always wanted, and difficult to stop. 
if run_movie =='y' | run_movie =='Y'
   movie(H_z_Movie,1,1)
end

% Remove large variables from workspace before saving
clear C_Ex C_Ey D_Hz 
clear E_x_n E_x_nmin1 E_y_n E_y_nmin1 H_z_n H_z_nmin1 mesh_pic

save cyl_fdtd_data 

