function [PI, Orig, Lambda, best_design] = myGenetic(f, interval, parents, TOL_GA, G, S, dv)
% Input parameters
% parents: scalar, number of design strings to preserve and breed
% TOL_GA: scalar, acceptable cost functon theshold to stop evolution
% G: scalar, maximum total generations
% S: scalar, number of total design strings per generation
% dv: scalar, number of design variables per string
% f: a function that can take in either a scalar or a column vector of
% design parameters
% Pretty sure there should be a search interval specified here
% Interval is now a dv x 2 matrix where each row is the interval for the mth
% design variable
%
%
% Outputs
% PI: G x S array where PI(g,s) is the cost of the s'th ranked design
% in the g'th generation
% Orig: G x S array where Orig(g,:) repesents the indices of each
% sorted entry from before they were sorted.
% Lambda: S x dv array of the most recent generations design strings

if mod(parents,2) ~= 0
    error('Number of parents must be even.')
end
if parents > S
    error('Number of parents to keep must be less than or equal to population size.')
end

K = parents/2;
% pop = interval(1) + (interval(2)-interval(1)) .* rand(dv, S); % Creates a dv x S matrix where each column is a design
pop = [];
for d = 1:dv
    pop_temp = interval(d,1) + (interval(d,2)-interval(d,1)) .* rand(1, S);
    pop(d,:) = pop_temp;
end
raw_cost = zeros(1,S); % Preallocate generation cost vector
Tstar = zeros(1, S);
Mstar = zeros(1, S);
Lstar = zeros(1, S);
best_parents = zeros(dv, parents); % Preallocate best parents
PI = zeros(G,S); % Preallocate PI
Orig = zeros(G,S); % Preallocate Orig

t0 = tic;
for g = 1:G % For each generation
    if g == 1
        for s = 1:S % For each design vector (columns of pop)
            func_vals = num2cell(pop(:,s)); % Temp storage to pass into function as cell array
            [raw_cost(s)] =  f(func_vals{:}); % Calculate cost of each design in the population
        end
    else
        % Iterate only through entries AFTER the parents
        for s = parents+1:S-parents % For each design vector (columns of pop)
            func_vals = num2cell(pop(:,s)); % Temp storage to pass into function as cell array
            [raw_cost(s)] =  f(func_vals{:}); % Calculate cost of each design in the population
        end
        raw_cost(1:parents) = sorted_cost(1:parents);
    end
    
    [sorted_cost, sorted_indices] = sort(raw_cost, 2); % Sort costs and get indices pre-sort
    
    PI(g,:) = sorted_cost; % Assign correct output names per project
    Orig(g,:) = sorted_indices;
    
    best_indices = sorted_indices(1:parents); % Get pre-sort indices of best parents to find the best design vector
    best_parents(:, 1:parents) = pop(:, best_indices); % Find best parents. Each column is a parent. each row is a design variable
    best_design = best_parents(:, 1)';
    children = zeros(dv, parents); % Preallocate and reset to zero the children array
    for k = 2*(1:K) % Even only iterator through every other parent
        for i = 1:-1:0 % Full iterator to make 2 children from each pair of 2 parents
            phi = rand(1); % Random weight [0, 1]
            child = best_parents(:,k-1).*phi + best_parents(:,k).*(1-phi); % Create child from 2 best parents
            children(:,k-i) = child; % Add 2 children from 2 best parents to array of all children
        end
    end
    
    % Add parents, children, and random to a new generation. Creates a dv x S matrix where each column is a design
    pop_temp3 = [];
    for d = 1:dv
        pop_temp2 = interval(d,1) + (interval(d,2)-interval(d,1)) .* rand(1, S-2*parents);
        pop_temp3(d,:) = pop_temp2;
    end
    pop = [best_parents, children, pop_temp3];
    
    %     % If the best design is within the tolerance for finding the global min
    %     % (only works if global min is 0)
    %     if sorted_cost(1)<TOL_GA
    %         break
    %     end
    
    if sorted_cost(1)<TOL_GA
        break;
    end
    
    myProgressBar(toc(t0), g, G, 'Genetic Algorithm');
end

Lambda = pop'; % % Assign correct output names per project
PI(g+1:end,:) = [];
Orig(g+1:end,:) = [];

end
