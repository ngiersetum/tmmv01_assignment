%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following script executes a VLM to study the aerodynamic properties
% of a finite swept wing.
% 
% Assignment for the course TMMV01 Aerodynamics at Linköping University,
% HT2023.
% 
% Group members:
%   Jan Zoellner (janzo813)
%   Niklas Gierse (nikgi434)
%
% Chosen assignment:
%   B747 - 12 deg dihedral
% 
% Geometrical data of the wing taken from the -400 model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear 
close all

%% Figure's positions
% ----------------------------------------------------------------------- %
scrsz = get(0,'ScreenSize');
scW=scrsz(3); % Screen Width (px)
scH=scrsz(4); % Screen Height (px)
% figure('Name','NameFigure','Position',[x0 y0 width height])
figure('Name','Wing','Position',[1 scH/4 scW/2 scH/2.5])
figure('Name','Distributions','Position',[scW/2 scH/4 scW/2 scH/2.5])

% Plots configuration
font=15;  % font size
lw=2;     % linewidth
szp=100;  % scatter size plot

%% Mesh Configuration
% ----------------------------------------------------------------------- %
npx = 3; % Number of panels in the streamwise direction
npy = 4; % Number of panels in the SEMI-spanwise direction
Num_panels = npx*npy

%% Flight Conditions
% ----------------------------------------------------------------------- %

alpha = 1*pi/180; % Angle of Attack (rad)
Uinf  = 20;       % Wind Speed (m/s)
rho   = 1.225;    % Air Density (kg/m3)

%% Wing Geometry
% ----------------------------------------------------------------------- %

Ale = deg2rad(37.5);            % Leading edge sweep angle (rad)
phi = deg2rad(12) + (1e-10);  % Dihedral angle (rad)
s   = 32;            % SEMI-span (m)
cr  = 14.84;            % Root chord (m)
ct  = 3.70 + (1e-10);  % Tip chord (m)
S   = 510.96;            % Wing surface (m^2)
AR  = 6.97;            % Wing Aspect Ratio (-)

% IMPORTANT: (1e-10) added at Dihedral Angle and Tip Chord in order to
% avoid numerical issues.

%% Panels´ Geometrical Properties
% ----------------------------------------------------------------------- %

tic
for j = ...:...
    
    % Local y coordinate of each slice 
    
    % Local chord at the coordinate y 
    
    % Local dx

    for i = ...:...      
    % Point A of each panel [x1n,y1n,z1n]
    Panels(i,j).x1n = abs(y_A_xy)*tan(Ale) + (1/4)*dx_A + dx_A*(i-1);           
    Panels(i,j).y1n = ...;
    Panels(i,j).z1n = ...;    
    
    % Point B of each panel [x2n,y2n,z2n]            
    % Control Points Locations [xm,ym,zm]

    end
     
    
end

% Time enlapsed in the geometric loop
T_geom_loop=toc; % (s)

%% MAIN LOOP
% ----------------------------------------------------------------------- %

% Set waitbar
h = waitbar(0,'Initializing Solver...');

% Reset Panel Counter
Panel=0;

tic

for i = ...   

    for j = ... % Solve only the Starboard wing due to symmetry
        
        % Panel Counter
        Panel=Panel+1;
        
        % Progress
        waitbar(Panel/(npx*npy),h,sprintf('Calculating %4.2f %%...',Panel/(npx*npy)*100))
        
        % Control point where the downwash from all the vortex is calculated
        xm = ...;
        ym = ...
        zm = ...;
     
        %% Starboard/Right Wing Vortex
        
        % Reset Vortex Counter
        Vortex=0; 
        
        for ii = ... 
            for jj = ... % Right wing loop
                
                % Vortex Counter
                Vortex=Vortex+1;
                
        % Induced Velocities function          
               
            end
        end
        
        %% Port/Left Wing Vortex
        
        % Reset Vortex Counter
        Vortex=0; 
        
        for ii = ...
            for jj = ... % Left wing loop 
                
                % Vortex Counter
                Vortex=Vortex+1;
                
         % Induced Velocities function       
                
            end
        end
           
    end
         
end

% Time enlapsed in the MAIN loop
T_MAIN_loop=toc; % (s) 

% Close wait bar
close(h)

%% Combine contributions from the Port and Startboard wings

Um = ...;
Vm = ...;
Wm = ...;

%% Assemble coefficients Matrix 
K = ...;

%% Solve the Vortex Strengths

GAMMAS = ...;

% Order Vortex Strengths in a matrix (for the Startboard wing)
for j = ... 

    for i = ...    
       
       V_strength( , ) = ...;
           
    end
         
end

%% Lift Coefficient

CL  = ...;

%% Spanwise Lift Distribution
% ----------------------------------------------------------------------- %
syms y
% Linear chord symbolic function (m)
% Mean aerodynamic chord (m)

for j = ... % Solve only the Starboard wing due to symmetry  

    % Local y coordinate of each slice    
    % Local y coordinate in adimensional form    
    % Local chord at the coordinate y (linear chord function)
     
    % Columwise sumatory of CL_nth
          
    % Columwise Lift
    
    % Columwise Lift Coefficient 
    
end

 





