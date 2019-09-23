function [cost] = myUAV(N_m, N_t, N_o, r_i, O_j, T_j, W_mt, W_mo, W_mm, w_t1, w_t2, w_o1, w_o2, w_m1, w_m2, a_1, a_2, b_1, b_2, c_1, c_2)
% Simulation
fprintf("Starting UAV swarm simulation...\n");
tic;

% Static variables
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

% From Table 3
w_1 = 70; % weight of mapping in net cost
w_2 = 10; % weight of time usage in net cost
w_3 = 20; % weight of agent losses in net cost

% Initial conditions
map_done = 0;
all_agents_dead = 0;
t = 1;
v_i = zeros(3, N_m);
ndims = 3;

while ( t < num_steps )
    i = 1;
    
    dmt_ij = [];
    nmt_ij = [];
    nmt_bar_ij = [];
    dmm_ij = [];
    nmm_ij = [];
    nmm_bar_ij = [];
    dmo_ij = [];
    nmo_ij = [];
    nmo_bar_ij = [];
    Nmt_i = [];
    Nmo_i = [];
    Nmm_i = [];
    Ntot_i = [];
    n_i = [];
    
    % Agent target interaction
    for d = 1:ndims
        dmt_ij_temp(:, :, d) = r_i(d, :) - T_j(d, :)'; % expanded subtraction
    end
    dmt_ij = vecnorm(dmt_ij_temp, 2, 3); % euclidean distance between ith agent and jth target
    uncap_ind = ~logical(dmt_ij <= agent_sight)'; % check if target can be mapped
    T_j = T_j(:, any(uncap_ind, 1)); % Keep only targets that have not been mapped
    for d = 1:ndims
        nmt_ij_temp(:, :, d) = T_j(d, :)' - r_i(d, :);
    end
    nmt_ij = nmt_ij_temp ./ vecnorm(nmt_ij_temp, 2, 3); % unit vector between ith agent and jth target
    attraction = w_t1.*exp(-a_1.*dmt_ij);
    repulsion = w_t2.*exp(-a_2.*dmt_ij);
    nmt_bar_ij = (attraction - repulsion).*nmt_ij; % cols are ith agent, rows are jth target, dims are [x, y, z]
    Nmt_i = sum(nmt_bar_ij, 1); % Agent total interaction vector with all targets
    Nmt_i1 = reshape(Nmt_i, 15, 3)';
    
    while ( i <= size(r_i, 2) )
        agents_dead = 0;
        
        
        % check if all targets are mapped
        if size(T_j, 2) == 0
            map_done = 1;
            break;
        end
        
        % Agent agent interaction
        j = 1;
        while ( j <= size(r_i, 2) )
            if i == j
                nmm_bar_ij{j, i} = [0; 0; 0];
                j = j + 1;
            else
                dmm_ij(j, i) = norm( r_i(:, i) - r_i(:, j) ); % euclidean distance between ith agent and jth agent
                if dmm_ij(j, i) <= crash_range
                    if i > j
                        r_i(:, i) = [];
                        v_i(:, i) = [];
                        r_i(:, j) = [];
                        v_i(:, j) = [];
                    else
                        r_i(:, j) = [];
                        v_i(:, j) = [];
                        r_i(:, i) = [];
                        v_i(:, i) = [];
                    end
                    agents_dead = 2;
                    break;
                else
                    nmm_ij{j, i} = ( r_i(:, j) - r_i(:, i) ) ./ ( norm(r_i(:, j) - r_i(:, i)) ); % unit vector between ith agent and jth agent
                    nmm_bar_ij{j, i} = ( w_m1 .* exp(-c_1.*dmm_ij(j, i)) - w_m2 .* exp(-c_2.*dmm_ij(j, i)) ) .* nmm_ij{j, i}; % interaction vector between ith agent and jth agent
                    j = j + 1;
                end
            end
        end
        
        if agents_dead == 2
            % check if any agents are still alive
            if size(r_i, 2) == 0
                all_agents_dead = 1;
                break;
            end
            continue;
        end
        
        % Agent obstacle interaction
        j = 1;
        while ( j <= size(O_j, 2) )
            dmo_ij(j, i) = norm( r_i(:, i) - O_j(:, j) ); % euclidean distance between ith agent and jth obstacle
            if dmo_ij(j, i) <= crash_range
                r_i(:, i) = [];
                v_i(:, i) = [];
                agents_dead = 1;
                break;
            else
                nmo_ij{j, i} = ( O_j(:, j) - r_i(:, i) ) ./ ( norm(O_j(:, j) - r_i(:, i)) ); % unit vector between ith agent and jth obstacle
                nmo_bar_ij{j, i} = ( w_o1 .* exp(-b_1.*dmo_ij(j, i)) - w_o2 .* exp(-b_2.*dmo_ij(j, i)) ) .* nmo_ij{j, i}; % interaction vector between ith agent and jth obstacle
                j = j + 1;
            end
        end
        
        if agents_dead == 1
            % check if any agents are still alive
            if size(r_i, 2) == 0
                all_agents_dead = 1;
                break;
            end
            continue;
        end
        
        % Each col is the interaction vector of the ith agent
        
        
        % Agent obstacle total interaction vector
        Nmo_i(:, i) = sum( [nmo_bar_ij{:}], 2 );
        
        % Agent agent total interaction vector
        Nmm_i(:, i) = sum( [nmm_bar_ij{:}], 2 );
        
        % Total interaction for ith agent
        Ntot_i(:, i) = W_mt .* Nmt_i(:, i) + W_mo .* Nmo_i(:, i) + W_mm .* Nmm_i(:, i);
        
        % Propulsive force direction vector
        n_i(:, i) = Ntot_i(:, i) ./ norm(Ntot_i(:, i));
        
        i = i + 1;
    end
    
    % check stop criteria
    if map_done || all_agents_dead
        break;
    end
    
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

% Create cost function (eqn 27)
Mstar = size(T_j, 2) / N_t; % fraction of unmapped targets
Tstar = t_final / t_f; % fraction of used time out of total time
Lstar = (N_m - size(r_i, 2)) / N_m; % fraction of crashed agents out of initial agents
cost = w_1*Mstar + w_2*Tstar + w_3*Lstar; % cost function

toc;
fprintf("Done with UAV swarm simulation!\n")
end

