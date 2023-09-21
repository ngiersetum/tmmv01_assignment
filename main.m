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
wingFig = figure('Name','Wing','Position',[1 scH/4 scW/2 scH/2.5]);
distFig = figure('Name','Distributions','Position',[scW/2 scH/4 scW/2 scH/2.5]);

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

Ale = deg2rad(37);              % Leading edge sweep angle (rad)
phi = deg2rad(12) + (1e-10);    % Dihedral angle (rad)
s   = 32;                       % SEMI-span (m)
cRoot  = 14.84;                 % Root chord (m)
cTip  = 3.70 + (1e-10);         % Tip chord (m)
S   = 541.2;                    % Wing surface (m^2)
AR  = 7.7;                      % Wing Aspect Ratio (-)

% IMPORTANT: (1e-10) added at Dihedral Angle and Tip Chord in order to
% avoid numerical issues.

%% Panels´ Geometrical Properties
% ----------------------------------------------------------------------- %
grid = zeros(1,3); % initialize vector of grid points for plotting the wing

tic
for j = 1:npy

    % Width of one slice
    dy = s/npy;

    % Local y coordinate of each slice (0 on inner edge)
    yA = dy .* (j-1);
    
    % Local chord at the coordinate y at both edges and the middle
    cA = cRoot - (cRoot-cTip).*(yA/s);
    cB = cRoot - (cRoot-cTip).*((yA+dy)/s);
    cC = (cA+cB)/2;
    
    % Local dx (length of one slice) at both edges and the middle
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

        % Save points for plotting the wing grid
        grid = [grid; yA*tan(Ale) + dxA*(i-1), yA, yA*tan(phi); 
            (yA+dy)*tan(Ale) + dxB*(i-1), yA+dy, (yA+dy)*tan(phi); 
            (yA+dy)*tan(Ale) + dxB*i, yA+dy, (yA+dy)*tan(phi); 
            yA*tan(Ale) + dxA*i, yA, yA*tan(phi); 
            nan, nan, nan]; % add NaN at the end to avoid connecting lines in the plot
    end
end

% Time enlapsed in the geometric loop
T_geom_loop=toc; % (s)

%% Plot the wing
figure(wingFig)

mat = cell2mat(struct2cell(panels));

% Vortex vertices and control points
scatter3(mat(1,:), mat(2,:), mat(3,:), 'r', 'filled');
hold on
scatter3(mat(1,:), -mat(2,:), mat(3,:), 'r', 'filled');
scatter3(mat(4,:), mat(5,:), mat(6,:), 'r', 'filled');
scatter3(mat(4,:), -mat(5,:), mat(6,:), 'r', 'filled');
scatter3(mat(7,:), mat(8,:), mat(9,:), 'g', 'filled');
scatter3(mat(7,:), -mat(8,:), mat(9,:), 'g', 'filled');

% Panel grid
plot3(grid(:,1), grid(:,2), grid(:,3), 'b');
plot3(grid(:,1), -grid(:,2), grid(:,3), 'b');

% Wing outline (starboard)
outline = [0 0 0; s*tan(Ale) s s*tan(phi); s*tan(Ale)+cTip s s*tan(phi); cRoot 0 0; 0 0 0];
plot3(outline(:,1), outline(:,2), outline(:,3), 'k');
plot3(outline(:,1), -outline(:,2), outline(:,3), 'k');


xlim([-s+15 s+25])
ylim([-s-5 s+5])
zlim([-s-5 s+5])


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
                coefficientsRight(panel, vortex).u = vCoeff(1);
                coefficientsRight(panel, vortex).v = vCoeff(2);
                coefficientsRight(panel, vortex).w = vCoeff(3);
            end
        end
        
        %% Port/Left Wing Vortex
        
        % Reset Vortex Counter
        vortex=0; 
        
        for ii = 1:npx
            for jj = 1:npy % Left wing loop 
                
                % Vortex Counter
                vortex=vortex+1;
                
                % Induced Velocities function       
                vCoeff = zeros(1,3);
                [vCoeff(1), vCoeff(2), vCoeff(3)] = calculateInducedVelocity(xm, ym, zm, panels(ii,jj).x1n, -panels(ii,jj).y1n, panels(ii,jj).z1n, panels(ii,jj).x2n, -panels(ii,jj).y2n, panels(ii,jj).z2n);
                coefficientsLeft(panel, vortex).u = vCoeff(1);
                coefficientsLeft(panel, vortex).v = vCoeff(2);
                coefficientsLeft(panel, vortex).w = vCoeff(3);
            end
        end
           
    end
         
end

% Time enlapsed in the MAIN loop
T_MAIN_loop=toc; % (s) 

% Close wait bar
close(h)

%% Combine contributions from the Port and Startboard wings

vMatRight = cell2mat(struct2cell(coefficientsRight));
vMatLeft = cell2mat(struct2cell(coefficientsLeft));

Um = vMatRight(1,:,:) + vMatLeft(1,:,:);
Vm = vMatRight(2,:,:) + vMatLeft(2,:,:);
Wm = vMatRight(3,:,:) + vMatLeft(3,:,:);

%% Assemble coefficients Matrix

% normalV = [0; sin(phi); cos(phi)];
normalU = ones(size(Um));
normalV = ones(size(Vm));
normalW = ones(size(Wm));

normalU = normalU .* 0;
normalV = normalV .* sin(phi);
normalW = normalW .* cos(phi);

K = dot([Um; Vm; Wm], [normalU; normalV; normalW]);
Kmat = reshape(K, (npx*npy), []);  % make npx*npy x npx*npy matrix out of 1xNxM tensor

%% Define right hand side

rhs = ones(npx*npy, 1);

rhs = rhs .* (-Uinf*sin(alpha)*cos(phi));

%% Solve the Vortex Strengths

gammas = linsolve(Kmat, rhs);

% Order Vortex Strengths in a matrix (for the Startboard wing)
for j = 1:npy
    for i = 1:npx
       strengths(i,j) = gammas(npx*(j-1) + i);
    end
end

%% Lift Coefficient - WE ARE HERE

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

 





