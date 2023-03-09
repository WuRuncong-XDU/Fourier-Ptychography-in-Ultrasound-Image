% -------------------------------------------------------------------------
% This script reference by: Marine Predators Algorithm source code
% programming: Afshin Faramarzi & Seyedali Mirjalili
% paper:
% A. Faramarzi, M. Heidarinejad, S. Mirjalili, A.H. Gandomi,
% Marine Predators Algorithm: A Nature-inspired Metaheuristic
% Expert Systems with Applications
% DOI: doi.org/10.1016/j.eswa.2020.113377
% Redistributions by Wu Runcong 2023/03/09.
% -------------------------------------------------------------------------
function [Top_predator_fit, Top_predator_pos, Convergence_curve] = MPA(SearchAgents, Max_iter, dim, des)
tic; fw = waitbar(0, 'MPA initialization.');

ub = [1, 2.4, 15];          % [n0 ratio loop]
lb = [0.1, 1.2, 3];
if size(ub, 2) == 1
    ub = ones(1, dim) * ub;
    lb = ones(1, dim) * lb;
end

Top_predator_pos = zeros(1, dim);
Top_predator_fit = inf;

Convergence_curve = zeros(1, Max_iter);
stepsize = zeros(SearchAgents, dim);
fitness = inf(SearchAgents, 1);

% Initialize first population of 'prey'
Boundary_no = size(ub, 2);
if Boundary_no == 1
    Prey = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
elseif Boundary_no > 1
    Prey = zeros(SearchAgents, dim);
    for i = 1 : dim
        ub_i = ub(i);
        lb_i = lb(i);
        Prey(:, i) = rand(SearchAgents, 1) .* (ub_i - lb_i) + lb_i;
    end
end

Xmin = repmat(ones(1, dim).*lb, SearchAgents, 1);
Xmax = repmat(ones(1, dim).*ub, SearchAgents, 1);

Iter = 0;
FADs = 0.2;
P = 0.5;

while Iter < Max_iter
    % ------------------- Detecting top predator -----------------
    for i = 1 : size(Prey, 1)
        
        Flag4ub = Prey(i, :) > ub;
        Flag4lb = Prey(i, :) < lb;
        Prey(i, :) = (Prey(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        
        fitness(i, 1) = Fitness(Prey(i,:), des);
        
        if fitness(i, 1) < Top_predator_fit
            Top_predator_fit = fitness(i, 1);
            Top_predator_pos = Prey(i, :);
        end
        
        waitbar((Iter*SearchAgents+i) / (Max_iter * SearchAgents), fw, ...
            ['MPA loop process - ' num2str(Iter*SearchAgents+i) ' ... ' num2str(floor(toc)) ' s']);
    end
    
    % ------------------- Marine Memory saving -------------------
    if Iter == 0
        fit_old = fitness;
        Prey_old = Prey;
    end
    
    Inx = (fit_old < fitness);
    Indx = repmat(Inx, 1, dim);
    Prey = Indx .* Prey_old + ~Indx .* Prey;
    fitness = Inx .* fit_old + ~Inx .* fitness;
    
    fit_old = fitness;
    Prey_old = Prey;
    
    % ------------------------------------------------------------
    Elite = repmat(Top_predator_pos, SearchAgents, 1);	% (Eq. 10)
    CF = (1 - Iter / Max_iter)^(2 * Iter / Max_iter);
    
    RL = 0.05 * levy(SearchAgents, dim, 1.5);            % Levy random number vector
    RB = randn(SearchAgents, dim);                       % Brownian random number vector
    
    for i = 1 : size(Prey, 1)
        for j = 1 : size(Prey, 2)
            R = rand();
            % ------------------ Phase 1 (Eq.12) -------------------
            if Iter < Max_iter / 3
                stepsize(i, j) = RB(i, j) * (Elite(i, j) - RB(i, j) * Prey(i, j));
                Prey(i, j) = Prey(i, j) + P * R * stepsize(i, j);
                
            % --------------- Phase 2 (Eqs. 13 & 14) ----------------
            elseif Iter > Max_iter / 3 && Iter < 2 * Max_iter / 3
                if i > size(Prey, 1) / 2
                    stepsize(i, j) = RB(i, j) * (RB(i, j) * Elite(i, j) - Prey(i, j));
                    Prey(i, j) = Elite(i, j) + P * CF * stepsize(i, j);
                else
                    stepsize(i, j) = RL(i, j) * (Elite(i, j) - RL(i, j) * Prey(i, j));
                    Prey(i, j) = Prey(i, j) + P * R * stepsize(i, j);
                end
                
            % ----------------- Phase 3 (Eq. 15) -------------------
            else
                stepsize(i, j) = RL(i, j) * (RL(i, j) * Elite(i, j) - Prey(i, j));
                Prey(i, j) = Elite(i, j) + P * CF * stepsize(i, j);
                
            end
        end
    end
    
    % ------------------ Detecting top predator ------------------
    for i = 1 : size(Prey, 1)
        Flag4ub = Prey(i, :) > ub;
        Flag4lb = Prey(i, :) < lb;
        Prey(i, :) = (Prey(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        
        fitness(i, 1) = Fitness(Prey(i,:), des);
        
        if fitness(i, 1) < Top_predator_fit
            Top_predator_fit = fitness(i, 1);
            Top_predator_pos = Prey(i, :);
        end
    end
    
    % ---------------------- Marine Memory saving ----------------
    if Iter == 0
        fit_old = fitness;
        Prey_old = Prey;
    end
    
    Inx = (fit_old < fitness);
    Indx = repmat(Inx, 1, dim);
    Prey = Indx .* Prey_old + ~Indx .* Prey;
    fitness = Inx .* fit_old + ~Inx .* fitness;
    
    fit_old = fitness;
    Prey_old = Prey;
    
    % ---------- Eddy formation and FADs?effect (Eq 16) -----------
    if rand() < FADs
        U = rand(SearchAgents, dim) < FADs;
        Prey = Prey + CF * ((Xmin + rand(SearchAgents, dim) .* (Xmax - Xmin)) .* U);
        
    else
        r = rand();
        Rs = size(Prey, 1);
        stepsize = (FADs * (1 - r) + r) * (Prey(randperm(Rs), :) - Prey(randperm(Rs), :));
        Prey = Prey + stepsize;
    end
    
    Iter = Iter + 1;
    Convergence_curve(Iter) = Top_predator_fit;
end

close(fw);