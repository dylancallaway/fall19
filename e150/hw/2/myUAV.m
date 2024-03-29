function [cost, Tstar, Lstar, Mstar] = myUAV(plot_flag, N_m, N_t, N_o, r_i, O_j, T_j, W_mt, W_mo, W_mm, w_t1, w_t2, w_o1, w_o2, w_m1, w_m2, a_1, a_2, b_1, b_2, c_1, c_2)
% Simulation
% fprintf("Starting UAV swarm simulation...\n");
tic;

% Static variables
ndims = 3;

% From Table 1
A_i = 1; % m^2. agent characteristic area
C_di = .25; % agent coefficient of drag
m_i = 10; % kg. agent mass
magF_pi = 200; % N, prop force magnitude
v_a = 0; %m/s. air velocity
p_a = 1.225; % kg/m^3. air density
dt = .2; % s. time step size
t_f = 60; % s. maximum task time
time_steps = 0:dt:t_f; % time step iterator
num_steps = length(time_steps); % time_step index iterator

% From Table 2
agent_sight = 5; % m. target mapping distance
crash_range = 2; % m. agent collision distance

% Initial conditions
t = 1;
v_i = zeros(3, N_m);
if isempty(O_j)
    O_j = [NaN; NaN; NaN];
end
map_done = 0;
all_agents_crashed = 0;

while ( t < num_steps )
    %% Target mapping
    for d = 1:ndims
        dmt_ij_temp(:, :, d) = r_i(d, :) - T_j(d, :)'; % expanded subtraction
    end
    dmt_ij = vecnorm(dmt_ij_temp, 2, 3); % euclidean distance. cols are agents, rows are targets
    
    % uncap_und is a N_m x N_t matrix where each col is the distance from
    % the jth target to each of the agents. 1 x N_t vector where each col is logical
    % wether or not that target should be kept (1 has been mapped)
    cap_ind = any(logical(dmt_ij <= agent_sight)', 1); % check if target has been mapped
    
    % T_j is 3 x N_t here. Once a target is mapped, all its rows are NaN
    T_j(:, cap_ind) = NaN; % NaN targets that have been mapped
    
    % Check map success criteria
    if all(isnan(T_j))
        map_done = 1;
        break;
    end
    
    %% Agent agent crashing
    for d = 1:ndims
        dmm_ij_temp(:, :, d) = r_i(d, :) - r_i(d, :)'; % expanded subtraction
    end
    dmm_ij = vecnorm(dmm_ij_temp, 2, 3); % euclidean distance between ith agent and jth agent
    dmm_ij(dmm_ij == 0) = Inf; % set agent self interaction to Inf so there are no false crashes due to 0 distance
    
    % row vector where each col is an agent that
    % crashed
    crash_ind_mm = any(logical(dmm_ij <= crash_range), 1); % get indices of all agents within crash range (N_m x N_m)
    r_i(:, crash_ind_mm) = NaN; % NaN agents that have crashed
    
    % check agents dead criteria
    if all(isnan(r_i))
        all_agents_dead = 1;
        break;
    end
    
    %% Agent obstacle crashing
    for d = 1:ndims
        dmo_ij_temp(:, :, d) = r_i(d, :) - O_j(d, :)'; % expanded subtraction
    end
    dmo_ij = vecnorm(dmo_ij_temp, 2, 3); % euclidean distance between ith agent and jth obstacle
    % row vector where each col is an agent that
    % crashed
    crash_ind_mo = any(logical(dmo_ij <= crash_range), 1); % get indices of all agents within crash range (N_o x N_m)
    r_i(:, crash_ind_mo) = NaN; % NaN agents that have crashed
    
    % check agents dead criteria
    if all(isnan(r_i))
        all_agents_dead = 1;
        break;
    end
    
    %% Agent world crashing
    crash_ind_mw = any((r_i < [-150; -150; -60]) | (r_i > [150; 150; 60]), 1); % Check if agents are out of world bounds
    r_i(:, crash_ind_mw) = NaN; % NaN the agents that went OOB
    
    if all(isnan(r_i))
        all_agents_dead = 1;
        break;
    end
    
    %% Agent target interaction
    for d = 1:ndims
        % T_j is N_t x 1 here
        % r_i is 1 x N_m here
        % nmt_ij_temp is N_t x N_m (after ndims is N_t x N_m x ndims)
        nmt_ij_temp(:, :, d) = T_j(d, :)' - r_i(d, :);
    end
    
    % cols are agents, rows are targets
    nmt_ij = nmt_ij_temp ./ vecnorm(nmt_ij_temp, 2, 3); % unit vector between ith agent and jth target
    attraction = w_t1.*exp(-a_1.*dmt_ij);
    repulsion = w_t2.*exp(-a_2.*dmt_ij);
    nmt_bar_ij = (attraction - repulsion).*nmt_ij; % cols are ith agent, rows are jth target, dims are [x, y, z]
    nmt_bar_ij(isnan(nmt_bar_ij)) = 0; % sum NaN always returns NaN, set NaNs to 0 for sum
    Nmt_i = sum(nmt_bar_ij, 1); % Agent total interaction vector with all targets
    
    % interaction vector for the ith agent and all targets
    Nmt_i = reshape(Nmt_i, N_m, 3)'; % reshape so each col is an agent
    
    %% Agent agent interaction
    for d = 1:ndims
        % r_i is N_m x 1 here
        % r_i is 1 x N_m here
        % nmm_ij_temp is N_m x N_m (after ndims is N_m x N_m x ndims)
        nmm_ij_temp(:, :, d) = r_i(d, :)' - r_i(d, :);
    end
    nmm_ij = nmm_ij_temp ./ vecnorm(nmm_ij_temp, 2, 3); % unit vector between ith agent and jth agent
    attraction = w_m1.*exp(-c_1.*dmm_ij);
    repulsion = w_m2.*exp(-c_2.*dmm_ij);
    nmm_bar_ij = (attraction - repulsion).*nmm_ij; % cols are ith agent, rows are jth agent, dims are [x, y, z]
    nmm_bar_ij(isnan(nmm_bar_ij)) = 0; % set NaNs to 0 for sum
    Nmm_i = sum(nmm_bar_ij, 1); % Agent total interaction vector with all other agents
    
    % interaction vector for the ith agent and all other agents
    Nmm_i = reshape(Nmm_i, N_m, 3)'; % reshape so each col is an agent
    
    %% Agent obstacle interaction
    
    % SEEMS LIKE ALL ARE CRASHING WHEN ONE IS CRASHING
    
    for d = 1:ndims
        % O_j is N_o x 1 here
        % r_i is 1 x N_m here
        % nmm_ij_temp is O_j x N_m (after ndims is O_j x N_m x ndims)
        nmo_ij_temp(:, :, d) = O_j(d, :)' - r_i(d, :);
    end
    nmo_ij = nmo_ij_temp ./ vecnorm(nmo_ij_temp, 2, 3); % unit vector between ith agent and jth agent
    attraction = w_o1.*exp(-b_1.*dmo_ij);
    repulsion = w_o2.*exp(-b_2.*dmo_ij);
    nmo_bar_ij = (attraction - repulsion).*nmo_ij; % cols are ith agent, rows are jth agent, dims are [x, y, z]
    nmo_bar_ij(isnan(nmo_bar_ij)) = 0; % set NaNs to 0 for sum
    Nmo_i = sum(nmo_bar_ij, 1); % Agent total interaction vector with all other agents
    
    % interaction vector for the ith agent and all obstacles
    Nmo_i = reshape(Nmo_i, N_m, 3)'; % reshape so each col is an agent
    
    %% Total interaction vector
    Ntot_i = W_mt .* Nmt_i + W_mo .* Nmo_i + W_mm .* Nmm_i;
    
    %% Propulsive force vector
    n_i = Ntot_i ./ vecnorm(Ntot_i, 2);
    
    %%
    % Each col is the ith agent
    F_pi = magF_pi .* n_i; % propulsive force on each agent
    v_i(isnan(v_i)) = 0;
    F_di = .5 .* p_a .* C_di .* A_i .* vecnorm(v_a - v_i) .* (v_a - v_i); % drag force on each agent
    F_t = F_pi + F_di; % sum forces (neglect gravity)
    a_i = F_t ./ m_i; % calculate current acceleration of each agent
    
    %% Plotting
    if (plot_flag)
        if t == 1 || mod(t, 40) == 0 || t == 200
            scatter3(r_i(1, :), r_i(2, :), r_i(3, :), 'b', 'filled');
            hold on
            scatter3(T_j(1, :), T_j(2, :), T_j(3, :), 'g', 'filled');
            scatter3(O_j(1, :), O_j(2, :), O_j(3, :), 'r', 'filled');
            hold off
            set(gca,'XLim',[-150 150],'YLim',[-150 150],'ZLim',[-60 60])
            xlabel("X")
            ylabel("Y")
            zlabel("Z")
            legend("UAVs", "Targets", "Obstacles")
            drawnow
            pause(0.02)
        end
    end
    
    v_i = v_i + a_i .* dt; % calculate velocity at next time step
    r_i = r_i + v_i .* dt; % calculate position at next time step
    
    t = t + 1;
end

% get final time taken
t_final = time_steps(t);

% From Table 3
w_1 = 70; % weight of mapping in net cost
w_2 = 10; % weight of time usage in net cost
w_3 = 20; % weight of agent losses in net cost

% Create cost function (eqn 27)
Mstar = 1 - sum(isnan(T_j), 2) / N_t; % fraction of unmapped targets
Tstar = w_2 * t_final / t_f; % fraction of used time out of total time
Lstar = sum(isnan(r_i), 2) / N_m; % fraction of crashed agents out of initial agents
Mstar = w_1.*Mstar(1, 1);
Lstar = w_3.*Lstar(1, 1);
cost = Mstar + Tstar + Lstar; % cost function
% fprintf("Unmapped: %.2f\nTime used: %.2f\nCrashed: %.2f\n\n", Mstar, Tstar, Lstar);

% fprintf("Done with UAV swarm simulation!\n");
end

