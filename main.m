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
numPanels = npx*npy;

%% Flight Conditions
% ----------------------------------------------------------------------- %

alpha = 1*pi/180; % Angle of Attack (rad)
Uinf  = 20;       % Wind Speed (m/s)
rho   = 1.225;    % Air Density (kg/m3)

%% Wing Geometry
% https://www.wikipedia.org/wiki/Boeing_747#747-400
% ----------------------------------------------------------------------- %

Ale = deg2rad(37);            % Leading edge sweep angle (rad)
phi = deg2rad(12) + (1e-10);  % Dihedral angle (rad)
s   = 32;            % SEMI-span (m)
cRoot  = 14.84;            % Root chord (m)
cTip  = 3.70 + (1e-10);  % Tip chord (m)
S   = 541.2;            % Wing surface (m^2)
AR  = 7.7;            % Wing Aspect Ratio (-)

% IMPORTANT: (1e-10) added at Dihedral Angle and Tip Chord in order to
% avoid numerical issues.

%% Panels´ Geometrical Properties
% ----------------------------------------------------------------------- %

tic
for j = 1:npy
    % Local y coordinate of each slice
    dy = s/npy;
    yA = dy .* (j-1);
    
    % Local chord at the coordinate y
    cA = cRoot - (cRoot-cTip).*(yA/s);
    cB = cRoot - (cRoot-cTip).*((yA+dy)/s);
    cC = (cA+cB)/2;
    
    % Local dx
    dxA = cA/npx;
    dxB = cB/npx;
    dxC = cC/npx;

    for i = 1:npx
        % Point A of each panel [x1n,y1n,z1n]
        panels(i,j).x1n = yA*tan(Ale) + (1/4)*dxA + dxA*(i-1);
        panels(i,j).y1n = yA;
        panels(i,j).z1n = yA*tan(phi);
    
        % Point B of each panel [x2n,y2n,z2n]
        panels(i,j).x2n = (yA+dy)*tan(Ale) + (1/4)*dxB + dxB*(i-1);
        panels(i,j).y2n = yA+dy;
        panels(i,j).z2n = (yA+dy)*tan(phi);

        % Control Points (Point C) Locations [xm,ym,zm]
        panels(i,j).xmn = (yA+(dy/2))*tan(Ale) + (3/4)*dxC + dxC*(i-1);
        panels(i,j).ymn = yA+(dy/2);
        panels(i,j).zmn = (yA+(dy/2))*tan(phi);

    end
end

% Time enlapsed in the geometric loop
T_geom_loop=toc; % (s)

%% Plot the wing points
mat = cell2mat(struct2cell(panels));

figure
scatter3(mat(1,:), mat(2,:), mat(3,:), "r");
hold on
scatter3(mat(4,:), mat(5,:), mat(6,:), "r");
scatter3(mat(7,:), mat(8,:), mat(9,:), "g");
xlim([0 40])
ylim([0 40])
zlim([0 40])

%% MAIN LOOP
% ----------------------------------------------------------------------- %

% Set waitbar
h = waitbar(0,'Initializing Solver...');

% Reset Panel Counter
panel=0;

tic

for i = 1:npx
    for j = 1:npy % Solve only the Starboard wing due to symmetry
        
        % Panel Counter
        panel=panel+1;
        
        % Progress
        waitbar(panel/(npx*npy),h,sprintf('Calculating %4.2f %%...',panel/(npx*npy)*100))
        
        % Control point where the downwash from all the vortex is calculated
        xm = panels(i,j).xmn;
        ym = panels(i,j).ymn;
        zm = panels(i,j).zmn;
     
        %% Starboard/Right Wing Vortex
        
        % Reset Vortex Counter
        vortex=0; 
        
        for ii = 1:npx 
            for jj = 1:npy % Right wing loop
                
                % Vortex Counter
                vortex=vortex+1;
                
                vCoeff = zeros(1,3);
                [vCoeff(1), vCoeff(2), vCoeff(3)] = calculateInducedVelocity(xm, ym, zm, panels(ii,jj).x1n, panels(ii,jj).y1n, panels(ii,jj).z1n, panels(ii,jj).x2n, panels(ii,jj).y2n, panels(ii,jj).z2n);
                vCoeff = transpose(vCoeff);
            end
        end
        
        %% Port/Left Wing Vortex
        
        % Reset Vortex Counter
        vortex=0; 
        
        for ii = 0:10
            for jj = 0:20 % Left wing loop 
                
                % Vortex Counter
                vortex=vortex+1;
                
                % Induced Velocities function       
                vCoeff = zeros(1,3);
                [vCoeff(1), vCoeff(2), vCoeff(3)] = calculateInducedVelocity(xm, ym, zm, panels(ii,jj).x1n, -panels(ii,jj).y1n, panels(ii,jj).z1n, panels(ii,jj).x2n, -panels(ii,jj).y2n, panels(ii,jj).z2n);
                vCoeff = transpose(vCoeff);
            end
        end
           
    end
         
end

% Time enlapsed in the MAIN loop
T_MAIN_loop=toc; % (s) 

% Close wait bar
close(h)

%% Combine contributions from the Port and Startboard wings

Um = 0;
Vm = 0;
Wm = 0;

%% Assemble coefficients Matrix 
K = 0;

%% Solve the Vortex Strengths

GAMMAS = 0;

% Order Vortex Strengths in a matrix (for the Startboard wing)
for j = 0:10 

    for i = 0:20    
       
       V_strength(i,j) = 0
           
    end
         
end

%% Lift Coefficient

CL  = 0;

%% Spanwise Lift Distribution
% ----------------------------------------------------------------------- %
syms y
% Linear chord symbolic function (m)
% Mean aerodynamic chord (m)

for j = 0:10 % Solve only the Starboard wing due to symmetry  

    % Local y coordinate of each slice    
    % Local y coordinate in adimensional form    
    % Local chord at the coordinate y (linear chord function)
     
    % Columwise sumatory of CL_nth
          
    % Columwise Lift
    
    % Columwise Lift Coefficient 
    
end

 





