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


size(row)

tiledlayout(1, 3);
nexttile

epsilon=dx;

colormap parula;

rng(0);
u = epsilon*randn(N*N,1);
u = u - mean(u);

surf(X, Y, reshape(u, N, N));
view(0, 90)
title("Random initial condition");

nexttile

u = readmatrix("data/AllenCahn50s_soln.csv");

surf(X, Y, reshape(u, N, N));
view(0, 90)
title("The clench")
colorbar('southoutside', 'Ticks', [-0.98, -0.5, 0, 0.5, 0.98], ...
    'TickLabels', {'-1', '-0.5', '0', '0.5', '1'})
set(gca,'ytick',[])

nexttile

u = readmatrix("gif/AllenCahnStep75.csv");

surf(X, Y, reshape(u, N, N));
view(0, 90)
title("Final state")
set(gca,'ytick',[])