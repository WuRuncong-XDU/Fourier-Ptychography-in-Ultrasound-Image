% -------------------------------------------------------------------------
% This script reference by: Chameleon Swarm Algorithm (CSA) source codes version 1.0
% Malik Sh. Braik,
% Chameleon Swarm Algorithm: A Bio-inspired Optimizer for Solving Engineering Design Problems
% Expert Systems with Applications
% DOI: https://doi.org/10.1016/j.eswa.2021.114685
% Redistributions by Wu Runcong 2023/03/09.
% -------------------------------------------------------------------------
function [fmin0, gPosition, cg_curve] = CSA(SearchAgent, Max_iter, dim, des)
tic; fw = waitbar(0, 'CSA initial');

ub = [1, 2.4, 15];          % [n0 ratio loop]
lb = [0.1, 1.2, 3];
if size(ub, 2) == 1
    ub = ones(1, dim) * ub;
    lb = ones(1, dim) * lb;
end

%% Convergence curve
cg_curve = zeros(1, Max_iter);

%% Initialize the first population of search agents
% number for both u and l
if size(ub, 2) == 1
    u_new = ones(1, dim) * ub;
    l_new = ones(1, dim) * lb;
else
    u_new = ub;
    l_new = lb;
end

% If each variable has a different l and u
chameleonPositions = zeros(SearchAgent, dim);
for i = 1 : dim
    u_i = u_new(i);
    l_i = l_new(i);
    chameleonPositions(:, i) = rand(SearchAgent, 1) .* (u_i - l_i) + l_i;
end

%% Evaluate the fitness of the initial population
fit = zeros(SearchAgent, 1);
for i = 1 : SearchAgent
    fit(i, 1) = Fitness(chameleonPositions(i, :), des);
end

%% Initalize the parameters of CSA
fitness = fit;	% Initial fitness of the random positions

[fmin0, index] = min(fit);

chameleonBestPosition = chameleonPositions;	% Best position initialization
gPosition = chameleonPositions(index, :);	% initial global position

v = 0.1 * chameleonBestPosition;            % initial velocity
v0 = 0.0 * v;

%% Start CSA
% Main parameters of CSA
rho = 1.0;
% p1 = 2.0;
% p2 = 2.0;
c1 = 2.0;
c2 = 1.80;
gamma = 2.0;
alpha = 4.0;
beta = 3.0;

%% Start CSA
for t = 1 : Max_iter
    a = 2590 * (1 - exp(-log(t)));
    omega = (1-(t/Max_iter))^(rho*sqrt(t/Max_iter));
    p1 = 2 * exp(-2 * (t / Max_iter)^2);
    p2 = 2 / (1 + exp((-t + Max_iter/2) / 100));
    
    mu = gamma * exp(-(alpha * t / Max_iter)^beta);
    
    ch = ceil(SearchAgent *  rand(1, SearchAgent));
    %% Update the position of CSA (Exploration)
    for i = 1 : SearchAgent
        if rand >= 0.1
            chameleonPositions(i, :) = chameleonPositions(i, :) + p1*(chameleonBestPosition(ch(i), :)...
                - chameleonPositions(i, :)) * rand() + p2 * (gPosition - chameleonPositions(i, :)) * rand();
        else
            for j = 1:dim
                chameleonPositions(i, j) = gPosition(j) + mu * ((ub(j) - lb(j)) * rand+lb(j)) * sign(rand - 0.50) ;
            end
        end
    end
    
    %% Rotation of the chameleons - Update the position of CSA (Exploitation)
    % Rotation 180 degrees in both direction or 180 in each direction
    % chameleonPositions = rotation(chameleonPositions, searchAgents, dim);
    
    %% Chameleon velocity updates and find a food source
    for i = 1 : SearchAgent
        v(i, :) = omega * v(i, :) + c1 * (chameleonBestPosition(i, :) - chameleonPositions(i, :)) * rand...
            + c2 * (gPosition - chameleonPositions(i, :)) * rand;
        chameleonPositions(i, :) = chameleonPositions(i, :) + (v(i, :).^2 - v0(i, :).^2) / (2 * a);
    end
    v0 = v;
    
    %% handling boundary violations
    for i = 1 : SearchAgent
        if chameleonPositions(i, :) < lb
            chameleonPositions(i, :) = lb;
        elseif chameleonPositions(i, :) > ub
            chameleonPositions(i, :) = ub;
        end
    end
    
    %% Relocation of chameleon positions (Randomization)
    for i = 1 : SearchAgent
        ub_ = sign(chameleonPositions(i, :) - ub) > 0;
        lb_ = sign(chameleonPositions(i, :) - lb) < 0;
        
        chameleonPositions(i, :) = (chameleonPositions(i, :) .* (~xor(lb_, ub_))) + ub .* ub_ +lb .* lb_;  % *2
        
        fit(i, 1) = Fitness(chameleonPositions(i, :), des);             % fobj (chameleonPositions(i, :))
        
        if fit(i) < fitness(i)
            chameleonBestPosition(i, :) = chameleonPositions(i, :);     % Update the best positions
            fitness(i) = fit(i);                                        % Update the fitness
        end
        waitbar(((t-1)*SearchAgent+i) / (Max_iter * SearchAgent), fw, ...
            ['CSA loop process - ' num2str((t-1)*SearchAgent+i) ' ... ' num2str(floor(toc)) ' s']);
    end
    
    
    %% Evaluate the new positions
    [fmin, index] = min(fitness);                       % finding out the best positions
    
    % Updating gPosition and best fitness
    if fmin < fmin0
        gPosition = chameleonBestPosition(index, :);	% Update the global best positions
        fmin0 = fmin;
    end
    
    %% Visualize the results
    cg_curve(t) = fmin0;	% Best found value until iteration t
    
end

% ngPosition = find(fitness == min(fitness));
% g_best = chameleonBestPosition(ngPosition(1), :);     % Solutin of the problem, the same as 'gPosition'.
% fmin0 = Fitness(g_best, des);
close(fw);
end