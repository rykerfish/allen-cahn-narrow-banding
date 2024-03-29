set(0,'defaulttextInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');
set(0,'defaultLineLineWidth',3);
set(0,'defaultAxesFontSize',35)

% clf
% clearvars

tFinal = 750;
timing_gap = 10;

adaptive_start_time = tic;
adaptive_times = zeros(1, tFinal/timing_gap);
t_index = 1;

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
    
    
    %%% 
    epsilon = dx;
    CFL = 0.001;
    dt = CFL*(dx);
    Nsteps = round(tFinal/dt);
    t_plot = 640000;
    
    u = readmatrix("data/AllenCahn50s_soln.csv");
    u = reshape(u, [N,N]);
    
    F = @(u_0) exp(dt) ./ ( sqrt( (1./(u_0.^2) + exp(2*dt) - 1)));
    
%     tplot = 10000;
%     tot_mass = zeros(1, round(Nsteps/tplot));
%     time = zeros(1, round(Nsteps/tplot));
    u_0 = u; % initial condition
    plot_index = 1;
    
    mass_index = 1;
    num_points_sampled = zeros(1, tFinal/timing_gap);
    cutoff = 0.90; % value we narrow band around
    sample_index = 1;
    for n = 1:Nsteps
        u_prev = u;
    
    
        % m is how far we extend the narrow bands in a square of cells, i.e. 
        % we use all cells in a (2m+1)*(2m+1) square around the center cell
        m = 4;
        % recompute the m nearest neighbors every this many time steps. this
        % seems to be okay since the interface moves so slowly, and saves a ton
        % of time since the nearest_neighbors function is a bit slow. can lower
        % it for better accuracy if you run into problems
        recompute_thresh = 5000; 
        if mod(n, recompute_thresh) == 1 %
            % get indicies of narrow band location
            [row, col] = find(abs(u) < cutoff);
            [row, col] = m_nearest_neighbors(m, row, col, N);
            num_points_sampled(sample_index) = length(row);
            sample_index = sample_index + 1;
        end
    
        u_star = u;
        for j = 1:length(row)
            x_i = row(j);
            y_i = col(j);
            u_center = u(x_i, y_i);
            if x_i == 1 % left edge
                u_west = u_center;
            else
                u_west = u(x_i-1, y_i);
            end
            if x_i == N % right edge
                u_east = u_center;
            else
                u_east=u(x_i+1, y_i);
            end
            if y_i == 1 % bottom edge
                u_south = u_center;
            else
                u_south = u(x_i,y_i-1);
            end
            if y_i == N % top edge
                u_north = u_center;
            else
                u_north=u(x_i, y_i+1);
            end
    
            % finite difference laplacian w 5 point stencil
            u_lap = (u_west + u_east + u_south + u_north - 4*u_center)*(1/h^2);
            u_star(x_i, y_i) = u(x_i, y_i) + epsilon^2*dt*u_lap; % forward euler
        end
    
        u_star2 = sign(u_star).*F(u_star);
    
        beta = (1/dt) * ( (sum(u_0 - u_star2, 'all')) / sum( (1/2)*u_star2 .* (u_star2-sqrt(2)) , 'all')  );
    
        u = u_star2 + dt*beta * (1/2)*(u_star2.*(u_star2 - sqrt(2)));
    
        if mod(n,t_plot) == 0
            name = strcat("gif/AllenCahnStep", num2str(plot_index), ".csv");
            plot_index = plot_index + 1;
            writematrix(u, name);
            u_old = u;
        end

        if n*dt > t_index*timing_gap
            current_time = tic;
            adaptive_times(t_index) = toc(adaptive_start_time);
            t_index = t_index + 1;
        end
    end
    adaptive_times(t_index) = toc(adaptive_start_time);

    adaptive_times = [0 adaptive_times]
    
writematrix(u, "AllenCahn750s_adaptive_soln.csv")
writematrix(num_points_sampled, "sampled_points.csv");

writematrix(adaptive_times, "adaptive_timing_750s_cutoff99_m3.csv");

matrix_start_time = tic;
matrix_times = zeros(1, tFinal/timing_gap);
t_index = 1;

    %%% Setup grid in x-y
    N = 64; % number of grid point along one direction
    x = linspace(0,1,N+1); % type 1 grid
    dx = x(2)-x(1);
    % convert to type 2 grid
    x = x(1:end-1) + dx/2.0;
    
    y = x;
    dy = dx;
    
    [X,Y] = meshgrid(x,y); % make 2D grid
    
    %%% make finite difference laplacian
    
    e = (1/dx^2)*ones(N,1);
    D2 = spdiags([e -2*e e], -1:1, N, N);
    
    % neumann boundary conditions
    D2(1,1) = -(1/dx^2);
    D2(end,end) = -(1/dx^2);
    
    
    I_n = speye(N);
    
    D_full = kron(I_n, D2) + kron(D2, I_n);
    
    %%% 
    epsilon = dx;
    CFL = 0.001;
    dt = CFL*(dx);
    Nsteps = round(tFinal/dt);
    t_plot = 10000;
    
    Lm = speye(N*N) - (dt)*epsilon^2*D_full;
    
    dLm = decomposition(Lm);
    
    
    u = readmatrix("data/AllenCahn50s_soln.csv");
    
    F = @(u_0) exp(dt) ./ ( sqrt( (1./(u_0.^2) + exp(2*dt) - 1)));
    
    u_0 = u; % initial condition
    
    i = 2;
    for n = 1:Nsteps
        u_star = dLm\u;
        u_star2 = sign(u_star).*F(u_star);
    
        beta = (1/dt) * ( (sum(u_0 - u_star2)) / sum( (1/2)*u_star2 .* (u_star2-sqrt(2)) )  );
    
        u = u_star2 + dt*beta * (1/2)*(u_star2.*(u_star2 - sqrt(2)));
    
        % mass plotting code
%         if mod(n, t_plot) == 0
%             tot_mass(i) = sum(u)*(1/dx^2);
%             time(i) = n*dt;
%             i = i + 1;
%         end
%     
%         if mod(n,t_plot) == 0
%            h = surf(X,Y,reshape(u,N,N));
%            set(h,'edgecolor','none')
%            zlim([-1.1 1.1])
%            title(['t = ' num2str(n*dt)])
%            view(0, 90);
%            colorbar
%            camlight
%            drawnow
%         end

        if n*dt > t_index*timing_gap
            current_time = tic;
            matrix_times(t_index) = toc(matrix_start_time);
            t_index = t_index + 1;
        end

    end
    matrix_times(t_index) = toc(matrix_start_time);

    matrix_times = [0 matrix_times]

    writematrix(u, "AllenCahn500s_matrix_soln.csv")

    writematrix(matrix_times, "matrix_times_500s_cutoff99_m3.csv");

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