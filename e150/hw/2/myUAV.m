function [cost] = myUAV(N_m, N_t, N_o, r_i, O_j, T_j, W_mt, W_mo, W_mm, w_t1, w_t2, w_o1, w_o2, w_m1, w_m2, a_1, a_2, b_1, b_2, c_1, c_2)
% Simulation
fprintf("Starting UAV swarm simulation...\n");
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
    %     dmt_ij = [];
    %     dmt_ij_temp = [];
    %     nmt_ij = [];
    %     nmt_bar_ij = [];
    %     dmm_ij = [];
    %     dmm_ij_temp = [];
    %     nmm_ij = [];
    %     nmm_bar_ij = [];
    %     dmo_ij_temp = [];
    %     dmo_ij = [];
    %     nmo_ij = [];
    %     nmo_bar_ij = [];
    %     Nmt_i = [];
    %     Nmo_i = [];
    %     Nmm_i = [];
    %     Ntot_i = [];
    %     n_i = [];
    %     nmt_ij_temp = [];
    %     nmo_ij_temp = [];
    %     nmm_ij_temp = [];
    %     attraction = [];
    %     repulsion = [];
    
    %% Target mapping
    for d = 1:ndims
        dmt_ij_temp(:, :, d) = r_i(d, :) - T_j(d, :)'; % expanded subtraction
    end
    dmt_ij = vecnorm(dmt_ij_temp, 2, 3); % euclidean distance between ith agent and jth target
    
    % uncap_und is a N_m x N_t matrix where each col is the distance from
    % the jth target to each of the agents. 1 x N_t vector where each col is logical
    % wether or not that target should be kept (1 has been mapped)
    cap_ind = all(logical(dmt_ij <= agent_sight)', 1); % check if target has been mapped
    
    % T_j is 3 x N_t here. Once a target is mapped, all its rows are NaN
    T_j(:, cap_ind) = NaN; % NaN targets that have been mapped
    %     dmt_ij(all(cap_ind, 1)', :) = NaN; % remove distances to mapped targets
    
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
    crash_ind_mm = all(logical(dmm_ij <= crash_range), 1); % get indices of all agents within crash range (N_m x N_m)
    r_i(:, crash_ind_mm) = NaN; % NaN agents that have crashed
    %     dmm_ij(:, all(crash_ind_mm, 1)) = NaN; % remove distances to crashed agents
    
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
    crash_ind_mo = all(logical(dmo_ij <= crash_range), 1); % get indices of all agents within crash range (N_o x N_m)
    r_i(:, crash_ind_mo) = NaN; % NaN agents that have crashed
    %     dmo_ij(:, all(crash_ind_mo, 1)) = NaN; % remove distances to crashed agents
    
    % check agents dead criteria
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
    nmt_ij = nmt_ij_temp ./ vecnorm(nmt_ij_temp, 2, 3); % unit vector between ith agent and jth target
    attraction = w_t1.*exp(-a_1.*dmt_ij);
    repulsion = w_t2.*exp(-a_2.*dmt_ij);
    nmt_bar_ij = (attraction - repulsion).*nmt_ij; % cols are ith agent, rows are jth target, dims are [x, y, z]
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
    F_di = .5 .* p_a .* C_di .* A_i .* norm(v_a - v_i) .* (v_a - v_i); % drag force on each agent
    F_t = F_pi + F_di; % sum forces (neglect gravity)
    a_i = F_t ./ m_i; % calculate current acceleration of each agent
    
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
Tstar = t_final / t_f; % fraction of used time out of total time
Lstar = 1 - sum(isnan(r_i), 2) / N_m; % fraction of crashed agents out of initial agents
cost = w_1*Mstar(1, 1) + w_2*Tstar + w_3*Lstar(1, 1); % cost function

toc;
fprintf("Done with UAV swarm simulation!\n")
end

