
%%
clc 
close all
clearvars

% Add matlab_functions to path
addpath('../advect_functions');

% Add Coleen specific functions to path
addpath(genpath('coleen_functions'));

% Neighbors
%load(['../../data/all_neighbors_2_360x200.mat'])
load(['../../data/index_neighborhood_100x50_2.mat']);

%% Time / draw
days = 30;             % sim time
num_hours = 24;     % dt
iters = ceil(days * 24 / num_hours);             

%% System Parameters
% Fish Swimming Speed
% Here BATS Tuna Speed 86-260km/day ~ 1m/s to 3m/s
% Biomass is approximately 1gww/m^3 ~ 
% Tuna consume ~5% of their mass per day
fish_speed = 1;       % speed in meters/second
%
fish_init = 1.0;        % fish biomass in g / m^3
food_init = 1000;       % food biomass in g / m^3
%
% Currents
u_max = 0.0;             % current speed in meters/second
v_max = 0.0;
% Time Discretization
dt   = num_hours/24;    % dt in days
dt_s = dt*24*60*60;     % dt in seconds
day  = 24*60*60;        % seconds in a day ( for dimensionalizing )
%
epsilon = 10e-15;      % tiny value, 
%
%% Grid setup
% Number of CELLS in each direction
m = 100;
n = 50;
% Domain's start and end locations ( distances )
a = 0;
b = 1000;
c = 0.0;
d = 500.0;
% Spatial step sizes
x = linspace(0, 1, m).^1.1;  % adjust exponent for nonlinearity
x_vec = a + (b - a) * x;
y = linspace(0, 1, n).^1.1;
y_vec = c + (d - c) * y;

y_vec = y_vec - (max(y_vec)/2);

% Grid
%[X, Y] = meshgrid( a : dx : b, c : dy : d );
[X, Y] = meshgrid( x_vec, y_vec);
Z = zeros(size(X));
mesh(X,Y,Z);

%dx = (b-a)/(m-1);     % kilometers
%dy = (d-c)/(n-1);     % kilometers
dx = diff(X, 1, 2);
dy = diff(Y, 1, 1);
dx = [dx dx(:,end) ];
dy = [dy ; dy(end,:)];

[ny, nx] = size(X);
cell_area = zeros(ny, nx);

for i = 1:ny-1
    for j = 1:nx-1
        % Four corners of the cell
        xcorn = [X(i,j)   X(i,j+1) X(i+1,j+1) X(i+1,j)];
        ycorn = [Y(i,j)   Y(i,j+1) Y(i+1,j+1) Y(i+1,j)];

        % Shoelace formula for polygon area
        cell_area(i,j) = 0.5 * abs( sum(xcorn .* circshift(ycorn,-1)) ...
                                  - sum(ycorn .* circshift(xcorn,-1)) );
    end
end
cell_area(end,:) = cell_area(end-1,:);
cell_area(:,end) = cell_area(:,end-1);

dx_m = dx * 1000.0;   % dx in meters
dy_m = dy * 1000.0;   % dy in meters

% Find neighbors
%neighborhood = find_all_neighbors_by_index( m, n );
%
%% Initial Density Generation for Fish
Fish = fish_init * ones(n, m);
Fish = zeroCornerMatrix(Fish);
Fish_init = Fish;
%
x0_fish = 300; % km
Fish = makeFood('humpfish', Fish_init, X, Y, x0_fish, 0, 1, 0);

%% Initial Food Generation
Food = food_init * ones(n, m);
Food = makeFood('gradient', Food, X, Y, 0, 0, 1, 0);
Food = fliplr(Food); % Flip gadient left to right
Food = Food .* food_init;

%% Island Masking
mask = ones(size(Food));
Fish = mask .* Fish;
Food = mask .* Food; 

%% Initial Current Field
current = zeros(n,m,2);
current(:,:,1) = -2;
current(:,:,2) = 0;
%% Convert to day from seconds
% Convert to m per day
current = current .* day;
% Fish speed to m per day
fish_speed = fish_speed .* day;

%% Mask Values
%current(:,:,1) = squeeze(current(:,:,1)) .* mask;
%current(:,:,2) = squeeze(current(:,:,2)) .* mask;
apparent_current = current;

Food = Food .* mask;
Fish = Fish .* mask;

%% Iteration Starts
for it = 1 : iters
    fish_old = sum(Fish(:));
     Fish = AdvectPredator( Fish,...
            Food, ...
            current, ...
            dt, ...
            dx_m, ...
            dy_m, ...
            neighborhood, ...
            fish_speed, ...
            mask, ...
            cell_area, ...
            m, ...
            n);
    fish_new = sum(Fish(:));
    fprintf('Iteration: %5d/%5d, diff=%6.5e \n', it, iters, fish_new-fish_old);
    
    %% Plot Fish Concentration
    titleString = ['Fish and food at step ', num2str(it)];

    plotFishFood(2,X,Y,Fish,Food,current,it,dt,titleString);
    % figure(23);
    % ax1 = gca;
    % img1 = pcolor(X, Y, Fish);
    % img1.FaceColor='interp';
    % img1.EdgeColor='none';
    % cmap1 = cmocean('dens');
    % colormap(ax1,cmap1); cb1 = colorbar;
    % ylabel(cb1, '$gww/m^3$','FontSize',14,'interpreter','latex');
    % set(cb1, 'TickLabelInterpreter','latex');
    % axis tight
    % pbaspect([2 1 1])
    % xlabel(ax1,'x [km]', 'interpreter','latex');
    % ylabel(ax1,'y [km]', 'interpreter','latex');
    % title(ax1,'Fish Concentration in $gww/m^3$', 'interpreter','latex');
    % clim([0. 0.1]);
    % drawnow;
    %caxis([0 0.1]);
    %
end
