%% MAE 423 Final Project: Cylinder
%  Delaney Miller and Jason Mulderrig
%  17 January 2016

% Numerically simulate flow with heat transfer over a circular cylinder
% at Re = 200.

%% Part (a)
% Select fluid, determine scale for cylinder, and establish and
% initialize stream function array with boundary conditions. Find the
% initial inviscid steady-state solution by solving Laplacian. Plot the
% streamlines in the computational domain (contours for constant Psi).

% Fluid properties (air at 300 K)
nu    = 1.568e-5;  % kinematic viscosity
alpha = 1.9e-5;    % thermal diffusivity
k = 0.024;         % thermal conductivity, W/(m*K)

% Cylinder properties
radius = 25;  % in grid points 
x0 = 100;     % x-coordinate of center
y0 = 101;     % y=coordinate of center

Re = 200;  % Reynolds number based on cylinder diameter
u0 = .01;    % free stream velocity (m/s)

% Grid spacing
m = 500;   % grid points in x
n = 201;   % grid points in y (extra point that can be center)
% for Re_h < 10, diffusion processes are adequately modeled
Re_h = Re/(2*radius);
h1 = Re_h*nu/u0;  % grid spacing based on nu
h2 = Re_h*alpha/u0; % grid spacing based on alpha
h = min(h1, h2);  % pick minimum grid spacing

% Define points in cylinder
[x, y] = defineCylinder(n,m,x0,y0,radius);

% Streamfunction boundary conditions
psi = zeros(n,m);  % zero everywhere
psi(1,:) =  -u0*h*(n-1)/2;  % top boundary is a streamline
psi(n,:) =  -psi(1,:);     % bottom boundary is a streamline
n_mid = 101;      % midpoint (where psi = 0)
for j = 1:n
    psi(j,:) = u0*h*(j - n_mid);
end

psi = reassignCylinder(psi,0,x,y);

% Initialize streamfunction distribution (inviscid) at t = 0;
F = 1.5;  % over-relaxation factor
c = 0;    % convergence flag
eps = zeros(n,m);

while c == 0
    psi_temp = psi;
    % Note: by using Gauss-Seidel iteration (replacing values in stream-
    % function array), the streamline at n = 101 will not be exactly 0.
    
    for j = 2:n-1
        for i = 2:m-1
            if sqrt((i-x0)^2 + (j-y0)^2) <= radius
                psi(j,i) = 0;
%             elseif j == 101
%                 psi(j,i) = 0;
            else
                psi(j,i) = psi(j,i) + (F/4)*(psi(j,i+1) + psi(j,i-1) + ...
                    psi(j+1,i) + psi(j-1,i) - 4*psi(j,i));
            end
            if psi_temp(j,i) == 0
                eps(j,i) = 0;
            else
                eps(j,i) = abs((psi(j,i) - psi_temp(j,i))/psi_temp(j,i));
            end
        end
    end
    
    % test convergence
    max_eps = max(max(eps));
    if max_eps < 0.01  % convergence
        c = 1;
    end
    % continue while loop if not converged
end

% Plot streamlines using contour map
figure(1)
contour(psi, 50)
axis equal
title('1a) Streamlines at t = 0')

%% INITIALIZE TEMPERATURE
% Initialize temperature distribution (400K at cylinder, 300K elsewhere)

temp_i = 300; % bulk flow temperature
temp_c = 400; % cylinder temperature
temp = temp_i * ones(n,m);

% temperature boundary conditions (if different from bulk flow)
temp(1:n, 1) = 300;  % inflow
temp(1, 1:m) = 300;  % top free lid
temp(n, 1:m) = 300;  % bottom free lid

% assign cylinder temperature
temp = reassignCylinder(temp, temp_c,x, y);

%% Part (b)
% Develop subroutine to 'turn on' no-slip condition to generate
% vorticity at cylinder boundary. Simulate the subsequent unsteady
% evolution of the flow. Create a short movie of the vorticity field
% illustrating initial evolution (several hundred time steps).

% Create arrays to store points in interior of cylinder
[x_int, y_int] = cylinderInterior(n,m,x0,y0,radius);
circ = cat(2,x,y);
circ_int = cat(2,x_int,y_int);

% Create arrays to store points on boundary of cylinder
circ_bound = setdiff(circ, circ_int,'rows');
x_bound = circ_bound(:,1);
y_bound = circ_bound(:,2);

% Initialize bulk vorticity (w = 0) at t = 0
w = zeros(n,m);  % w = 0 everwhere at t = 0 (inviscid)

% Initialize convergence error array
eps_t = zeros(n,m);

% Pick time step
steps = 0; 
time_steps = 40000;
dt = 0.5*h/(3*u0);  % we use 3*u0 as the approximate maximum velocity

% Store vorticity data for movie
frame_spacing = 50;  % save every 10th vorticity frame
w_array = zeros(n,m,time_steps/frame_spacing);
w_max = 0;
w_min = 0;

% STORE TEMPERATURE DATA FOR MOVIE
temp_array = zeros(n,m,time_steps/frame_spacing);
temp_max = temp_c;
temp_min = temp_i;

% STORE HEAT TRANSFER DATA
q_dt = zeros(1,time_steps);
q_sum = zeros(1,time_steps);
q_total = 0;
perimeter = length(x_bound);  % perimeter of cylinder (length of boundary)

% Create for loop to simulate unsteady flow
for t = 1:time_steps
    
    % Set vorticity at solid boundary of cylinder (no slip)
    for b = 1:length(x_bound)
        iw = x_bound(b);
        jw = y_bound(b);
        % sum nearest neighbors
        w(jw,iw) = -2*(psi(jw+1, iw) + psi(jw-1, iw) ...
            + psi(jw, iw-1) + psi(jw, iw+1)) / h^2;  
    end
    
    % find the velocities u and v based on psi
    u = u0*ones(n,m);  % free stream velocity
    v = zeros(n,m);
    for b = 2:m-1  % can ignore boundaries because w = 0 there
        for a = 2:n-1
            u(a,b) = (psi(a+1,b) - psi(a-1,b))/(2*h);
            v(a,b) = (psi(a,b-1) - psi(a,b+1))/(2*h);
        end
    end
    
    % advance vorticity AND TEMPERATURE in convective fluid (for t+1)
    w_prev = w;
    temp_prev = temp;
    for j = 2:n-1
        for i = 2:m-1
            if sqrt((i-x0)^2 + (j-y0)^2) > radius
                % Calculate Laplacian of vorticity
                w_lp = (w_prev(j+1,i) + w_prev(j-1,i) + w_prev(j,i+1) + ...
                    w_prev(j,i-1) - 4*w_prev(j,i))/h^2;
                
                % CALCULATE LAPLACIAN OF TEMPERATURE
                temp_lp = (temp_prev(j+1,i) + temp_prev(j-1,i) + ...
                    temp_prev(j,i+1) + temp_prev(j,i-1) - ...
                    4*temp_prev(j,i)) / h^2;

                % Find advective terms (vorticity and TEMPERATURE)
                if u(j,i) < 0
                    duw = u(j,i+1)*w_prev(j,i+1) - u(j,i)*w_prev(j,i);
                    dut = u(j,i+1)*temp_prev(j,i+1) - u(j,i)*temp_prev(j,i);
                else
                    duw = u(j,i)*w_prev(j,i) - u(j,i-1)*w_prev(j,i-1);
                    dut = u(j,i)*temp_prev(j,i) - u(j,i-1)*temp_prev(j,i-1);
                end
                
                if v(j,i) < 0
                    dvw = v(j+1,i)*w_prev(j+1,i) - v(j,i)*w_prev(j,i);
                    dvt = v(j+1,i)*temp_prev(j+1,i) - v(j,i)*temp_prev(j,i);
                else
                    dvw = v(j,i)*w_prev(j,i) - v(j-1,i)*w_prev(j-1,i);
                    dvt = v(j,i)*temp_prev(j,i) - v(j-1,i)*temp_prev(j-1,i);
                end
            
                % Find new vorticity
                w(j,i) = w_prev(j,i) + dt*(-duw/h - dvw/h + nu*w_lp);
                
                % FIND NEW TEMPERATURE
                temp(j,i) = temp_prev(j,i) + ...
                    dt*(-dut/h - dvt/h + alpha*temp_lp);
            else
                w(j,i) = w_prev(j,i);
                temp(j,i) = temp_prev(j,i);
            end
        end
    end
    % Make sure that vorticity inside cylinder is 0
    % w = reassignCylinder(w,0,x_int,y_int);
    
    % Save vorticity array
    if rem(t,frame_spacing) == 0
        w_array(:,:,t/frame_spacing) = w;
        temp_array(:,:,t/frame_spacing) = temp;
    end
    w_max = max(max(w));
    w_min = min(min(w));
    temp_max = max(max(temp));
    temp_min = min(min(temp));

    % Find the new value of psi at this time step (n+1)
    % Advance psi by numerically solving the Poisson equation
    c = 0;
    psi_prev = psi;
    count = 0;
    while c == 0
        count = count + 1;
        if count == 100
            break;
        end
        psi_t = psi;
        for j = 2:n-1
            for i = 2:m-1
                 if sqrt((i-x0)^2 + (j-y0)^2) <= radius
                    psi(j,i) = 0;
                  else
                    psi(j,i) = psi(j,i) + (F/4)*(psi(j,i+1) + psi(j,i-1) + ...
                    psi(j+1,i) + psi(j-1,i) + 4*h^2*w(j,i) - 4*psi(j,i));
                end
                % Find convergence error
                if psi_t(j,i) == 0
                    eps_t(j,i) = 0;
                else
                    eps_t(j,i) = abs((psi(j,i) - psi_t(j,i))/psi_t(j,i));
                end
            end
        end
        
        % test convergence
        max_eps_t = max(max(eps_t));
        if max_eps_t < 0.05   % convergence
            c = 1;
        end
    % continue while loop if not converged
    end
    
    % outflow boundary condition (assume vorticity convects out)
    for e = 1:n
        psi(e,m) = 2*psi_prev(e,m-1) - psi_prev(e,m-2);
        w(e,m) = 2*w_prev(e,m-1) - w_prev(e,m-2);
        temp(e,m) = temp_prev(e,m-1);
    end
    
    % Make sure psi = 0 for cylinder and cylinder boundary
    % psi = reassignCylinder(psi,0,x,y);
    
    % CALCULATE HEAT TRANSFER
    q_t = 0;
    for q = 1:perimeter
        % Node to the left
        if x_bound(q) >= 75 && x_bound(q) <= 83 && ...
                y_bound(q) >= 84 && y_bound(q) <= 118
            q_t = q_t + k * (temp(y_bound(q), x_bound(q)-1)...
                - temp(y_bound(q), x_bound(q))) / h;
        end
        % Node to the right
        if x_bound(q) >= 117 && x_bound(q) <= 125 && ...
                y_bound(q) >= 84 && y_bound(q) <= 118
            q_t = q_t + k * (temp(y_bound(q), x_bound(q)+1)...
                - temp(y_bound(q), x_bound(q))) / h;
        end
        % Node to the top
        if x_bound(q) >= 83 && x_bound(q) <= 117 && ...
                y_bound(q) >= 76 && y_bound(q) <= 83
            q_t = q_t + k * (temp(y_bound(q)-1, x_bound(q))...
                - temp(y_bound(q),x_bound(q))) / h;
        end
        % Node to the bottom
        if x_bound(q) >= 83 && x_bound(q) <= 117 && ...
                y_bound(q) >= 119 && y_bound(q) <= 126
            q_t = q_t + k * (temp(y_bound(q)+1, x_bound(q))...
                - temp(y_bound(q), x_bound(q))) / h;
        end
    end
    
    q_dt(1,t) = q_t; % in W/m^2
    q_total = q_total + dt*q_t;  % in J/m^2
    q_sum(1,t) = q_total; % in J/m^2
    
    % Next time step
    steps = steps + 1;
    
end

%% CALCULATE AND PLOT HEAT TRANSFER

% Create time array
time = 1:time_steps;
time = time*dt;

% Plot the heat transfer rate from the cylinder at a given time
figure()
plot(time,q_dt);
xlabel('Time elapsed (seconds)')
ylabel('Heat Transfer Rate (W/m^2)')
title('Heat transfer rate from cylinder over time')

% Plot the total heat transfer from cylinder over time
figure()
plot(time,q_sum);
xlabel('Time elapsed (seconds)')
ylabel('Total heat transfer from cylinder (J/m^2)')
title('Total heat transfer from cylinder over time')

%% Make a short movie of vorticity field

map = colormap('jet');
frames = time_steps/frame_spacing;
w_min2 = w_min/10;
w_max2 = w_max/10;
for z = 1:frames
    w_frame = w_array(:,:,z);
    index = fix((w_frame-w_min2)/(w_max2-w_min2)*64) + 1;
    index(index<1) = 1;
    index(index>64) = 64;
    Frames(z) = im2frame(index,map);
end

vorticity = VideoWriter('vorticity1_16.avi');
vorticity.FrameRate = 30;  % make sure to change fps in part c
open(vorticity)
writeVideo(vorticity,Frames)
close(vorticity)

%% TEMPERATURE VIDEO

for z = 1:frames
    temp_frame = temp_array(:,:,z);
    index = fix((temp_frame - temp_min)/(temp_max - temp_min)*64) + 1;
    index(index < 1) = 1;
    index(index > 64) = 64;
    Frames_Temp(z) = im2frame(index,map);
end

tempVideo = VideoWriter('temperature1_16.avi');
tempVideo.FrameRate = 30;  % make sure to change fps in part c
open(tempVideo)
writeVideo(tempVideo,Frames_Temp)
close(tempVideo)

%% Part (c)
% Estimate the Strouhal number of the developed vortex street and
% compare to experimental results for circular cylinder at Re = 200.
% Expect Strouhal number of approximately 0.2.

freq_obs = 10/53; % # of vortices shed / time elapsed in video
fps = 5;  % frames per second in video
tsf = frame_spacing;   % time steps per frame
freq = freq_obs / (fps * tsf * dt);
d = 2 * radius * h;  % diameter of cylinder (m)
St = freq*d/u0;  % Strouhal number estimate

% Our results yielded an estimated Strouhal number of 0.2264, which is
% a litlte larger than the experimental result of 0.2.