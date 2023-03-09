% -------------------------------------------------------------------------
% Refernce: Grasshopper Optimization Algorithm (GOA) source codes demo V1.0
% Author and programmer: Seyedali Mirjalili
% Main paper: S. Saremi, S. Mirjalili, A. Lewis
% Grasshopper Optimisation Algorithm: Theory and Application
% DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.01.004
% Redistributions by Wu Runcong 2023/03/09.
% -------------------------------------------------------------------------
function [TargetFitness, TargetPosition, Convergence_curve, Trajectories, fitness_history, position_history] = GOA(N, Max_iter, dim, des)

tic; fw = waitbar(0, 'GOA initialization.');

ub = [1, 2.4, 15];      % [n0 ratio loop]
lb = [0.1, 1.2, 3];

if size(ub, 1) == 1
    ub = ones(dim, 1) .* ub';
    lb = ones(dim, 1) .* lb';
end

flag = 0;
if (rem(dim, 2) ~= 0)	% this algorithm should be run with a even number of variables. This line is to handle odd number of variables
    dim = dim + 1;
    ub = [ub; 1];
    lb = [lb; -1];
    flag = 1;
end

% Initialize the population of grasshoppers
% GrassHopperPositions = initialization(N, dim, ub, lb);
if size(ub, 1) == 1
    GrassHopperPositions = rand(N, dim).*(ub - lb) + lb;
elseif size(ub, 1) > 1
    GrassHopperPositions = zeros(N, dim);
    for i = 1 : dim
        high = ub(i);
        low = lb(i);
        GrassHopperPositions(:, i) = rand(1, N) .* (high - low) + low;
    end
end
GrassHopperFitness = zeros(1, N);

fitness_history = zeros(N, Max_iter);
position_history = zeros(N, Max_iter, dim);
Convergence_curve = zeros(1, Max_iter);
Trajectories = zeros(N, Max_iter);

cMax = 1;
cMin = 0.00004;

% Calculate the fitness of initial grasshoppers
for i = 1 : size(GrassHopperPositions, 1)
    GrassHopperFitness(1, i) = Fitness(GrassHopperPositions(i, :), des);
    
    fitness_history(i, 1) = GrassHopperFitness(1, i);
    position_history(i, 1, :) = GrassHopperPositions(i, :);
    Trajectories(:, 1) = GrassHopperPositions(:, 1);
end

[sorted_fitness, sorted_indexes] = sort(GrassHopperFitness);

% Find the best grasshopper (target) in the first population
Sorted_grasshopper = zeros(N, dim);
for newindex = 1 : N
    Sorted_grasshopper(newindex, :) = GrassHopperPositions(sorted_indexes(newindex), :);
end

TargetPosition = Sorted_grasshopper(1, :);
TargetFitness = sorted_fitness(1);
Convergence_curve(1) = TargetFitness;

% Main loop
GrassHopperPositions_temp = zeros(N, dim);
for l = 2 : Max_iter    % Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
    c = cMax - l * ((cMax - cMin) / Max_iter);	% Eq. (2.8) in the paper
    
    for i = 1 : size(GrassHopperPositions, 1)
        temp = GrassHopperPositions';
        for k = 1 : 2 : dim
            S_i = zeros(2, 1);
            for j = 1 : N
                if i ~= j
                    Dist = distance(temp(k : k+1, j), temp(k : k+1, i));            % Calculate the distance between two grasshoppers
                    
                    r_ij_vec = (temp(k : k+1, j) - temp(k : k+1, i)) / (Dist+eps);	% xj-xi/dij in Eq. (2.7)
                    xj_xi = 2 + rem(Dist, 2);                                       % |xjd - xid| in Eq. (2.7) 
                    
                    s_ij = ((ub(k:k+1) - lb(k:k+1)) * c / 2) * S_func(xj_xi) .* r_ij_vec;	% The first part inside the big bracket in Eq. (2.7)
                    S_i = S_i + s_ij;
                end
            end
            S_i_total(k:k+1, :) = S_i;
        end
        
        X_new = c * S_i_total'+ (TargetPosition);   % Eq. (2.7) in the paper      
        GrassHopperPositions_temp(i, :) = X_new';         
    end
    % GrassHopperPositions
    GrassHopperPositions = GrassHopperPositions_temp;
    
    for i = 1 : size(GrassHopperPositions, 1)
        % Relocate grasshoppers that go outside the search space 
        Tp = GrassHopperPositions(i, :) > ub';
        Tm = GrassHopperPositions(i, :) < lb';
        GrassHopperPositions(i, :) = (GrassHopperPositions(i, :) .* (~(Tp + Tm))) + ub' .* Tp + lb' .* Tm;
        
        % Calculating the objective values for all grasshoppers
        GrassHopperFitness(1, i) = Fitness(GrassHopperPositions(i, :), des);
        
        fitness_history(i, l) = GrassHopperFitness(1, i);
        position_history(i, l, :) = GrassHopperPositions(i, :);
        
        Trajectories(:, l) = GrassHopperPositions(:, 1);
        
        % Update the target
        if GrassHopperFitness(1, i) < TargetFitness
            TargetPosition = GrassHopperPositions(i, :);
            TargetFitness = GrassHopperFitness(1, i);
        end
        
        CurrentTimes = (l - 1) * N + i;
        waitbar(CurrentTimes / (N * Max_iter), fw, ...
            ['GOA loop process - ' num2str(CurrentTimes) ' ... ' num2str(floor(toc)) ' s']);
    end
        
    Convergence_curve(l) = TargetFitness;
end

if (flag == 1)
    TargetPosition = TargetPosition(1 : dim-1);
end

close(fw);
