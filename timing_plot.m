set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

N = 64; % number of grid point along one direction

%%% Setup grid in x-y
x = linspace(0,1,N+1); % type 1 grid
dx = x(2)-x(1);
% convert to type 2 grid
x = x(1:end-1) + dx/2.0;

y = x;
dy = dx;
h = dx;

[X,Y] = meshgrid(x,y); % make 2D grid

u = readmatrix("adaptive_timing_750s_cutoff99_m3.csv");
v = readmatrix("matrix_times_500s_cutoff99_m3.csv");
length(v);

times = linspace(0, 750, 76);
plot(times, u);
hold on
plot(times, v);
hold off

ylabel("Computation time (s)")
xlabel("Simulation time (s)")
title("Runtime comparison versus original method")
xlim([0, 750]);

legend("Narrow banding, m=4, $\lambda$=0.9", "Finite difference Laplacian matrix", 'Location', 'southeast')

figure(2)

points = readmatrix("sampled_points.csv");

times = linspace(0, 750, 9600);
plot(times, points);
ylabel("Number of points")
xlabel("Simulation time (s)")
title("Number of points used per iteration, m=4, $\lambda$=0.9")
