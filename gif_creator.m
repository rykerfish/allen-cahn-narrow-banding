u = readmatrix("gif/AllenCahnStep1.csv");
v = readmatrix("gif/AllenCahnStep50.csv");

sum(abs(u-v), 'all')
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



for i = 1:50
    
    clf
    filename = strcat("gif/AllenCahnStep", num2str(i), ".csv");
    u = readmatrix(filename);
    img = surf(X,Y,reshape(u,N,N));
    set(img,'edgecolor','none')
    zlim([-1.1 1.1])
    title(['t = ' num2str(i*10)])
    view(0, 90);
    colorbar
    camlight
    
    if i == 1
        gif('allencahn.gif');
    else
        gif
    end

end