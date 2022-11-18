set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

N = 64; % number of grid point along one direction
u = readmatrix("data/AllenCahn50s_soln.csv");
u = reshape(u, N, N);

%%% Setup grid in x-y
x = linspace(0,1,N+1); % type 1 grid
dx = x(2)-x(1);
% convert to type 2 grid
x = x(1:end-1) + dx/2.0;

y = x;
dy = dx;
h = dx;

[X,Y] = meshgrid(x,y); % make 2D grid

cutoff = 0.85;
[row, col] = find(abs(u) < cutoff);

m = 2;
[row, col] = m_nearest_neighbors(m, row, col, N);

size(row)

scaled_row = row/N;
scaled_col = col/N;

tiledlayout(1, 2);
nexttile

colormap('bone');
img = surf(X,Y,reshape(u,N,N)-1.1);
set(img,'edgecolor','none')
zlim([-2.1 1.1])
view(0, 90);
camlight
drawnow

nexttile

points = plot(scaled_row', scaled_col', "LineStyle","none", 'Marker','.', 'Color', [220/255 72/255 72/255]);
rotate(points, [0.5, 0.5, 0], 180)
points.MarkerSize = 20;

hold on


colormap('bone');
img = surf(X,Y,reshape(u,N,N)-1.1);
set(img,'edgecolor','none')
zlim([-2.1 1.1])
title('Sampled points with m=2, $\lambda=0.85$')
view(0, 90);
camlight
drawnow

% m: number of surrounding neighbors to add
% row, col: indicies of cells being targeted
% N: number of cells in the total grid
function [new_row, new_col] = m_nearest_neighbors(m, row, col, N)

    % preallocate for speed
    new_row = zeros(N*N, 1);
    new_col = zeros(N*N, 1);

    index = 1; % maintain index of next empty spot in arrays
    for k = 1:length(row)
        for i = -m:m
            for j = -m:m
                cell_row = row(k) + i;
                cell_col = col(k) + j;
                
                % don't add cell if its beyond the boundaries
                if cell_row < 1 || cell_col < 1 || cell_row > N || cell_col > N
                    continue
                end

                new_row(index) = cell_row;
                new_col(index) = cell_col;
                index = index + 1;
            end
        end
    end

    % leaves only unique pairings so we don't double count cells and cuts
    % off the zero entry leftover from preallocating the array
    A = [new_row new_col];
    A = unique(A, 'rows', 'sorted');
    A = A(2:end, :);

    new_row = A(:, 1);
    new_col = A(:, 2);

end