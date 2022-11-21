%%% Setup grid in x-y
N = 64; % number of grid point along one direction
x = linspace(0,1,N+1); % type 1 grid
dx = x(2)-x(1);
% convert to type 2 grid
x = x(1:end-1) + dx/2.0;

y = x;
dy = dx;
h = dx;


[X,Y] = meshgrid(x,y); % make 2D grid

u = readmatrix("gif/AllenCahnStep50.csv");
surf(X, Y, u);



for i = 1:75
    
    clf
    filename = strcat("gif/AllenCahnStep", num2str(i), ".csv");
    u = readmatrix(filename);
    zlim([-1.1 1.1])
    img = surf(X,Y,reshape(u,N,N));
    set(img,'edgecolor','none')
    title(['t = ' num2str(i*10)])
    view(0, 90);
    colorbar('Ticks', [-0.97, -0.5, 0, 0.5, 0.97], ...
    'TickLabels', {'-1', '-0.5', '0', '0.5', '1'})
    camlight
    
    if i == 1
        gif('allencahn.gif');
    else
        gif
    end

end