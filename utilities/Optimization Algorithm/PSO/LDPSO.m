function [BestFitness, BestPosition, PSOBestCurve] = LDPSO(c1, c2, wmax, wmin, MaxDT, D, N, xmax, xmin, vmax, vmin, des, xInitial)
% 迭代次数
irNum = N + N * MaxDT;
fw = waitbar(0, 'Process ... ');
tic

% c1学习因子1; c2学习因子2; w惯性权重; MaxDT最大迭代次数
% D搜索空间维数（未知数个数; N初始化群体个体数目;）

% -------------- 初始化种群的个体(可以在这里限定位置和速度的范围) --------------
x = zeros(N, D);
v = zeros(N, D);
for i = 1 : N
    for j = 1 : D
       x(i, j) = rand * (xmax(j) - xmin(j)) + xmin(j);	% 随机初始化位置
    end
    v(i, 1:D) = rand * (vmax - vmin) + vmin;            % 随机初始化速度
end
x(:, 1:end-1) = roundn(x(:, 1:end-1), -2);
x(:, end) = round(x(:, end));
if exist('xInitial', 'var')
    for i = 1 : length(xInitial)
        x(i, :) = xInitial{i};
    end
end

% -------------- 先计算各个粒子的适应度，并初始化Pi和Pg --------------
p = zeros(1, N);    % 局部最优解（输出）
y = zeros(size(x)); % 局部最优解（参数）
for i = 1 : N
    p(i) = Fitness(x(i, :), des);
    y(i,: ) = x(i, :);
    waitbar(i / irNum, fw, ['Loop process1 - ' num2str(i) ' ... ' num2str(floor(toc)) ' s']);
end
BestPosition = x(find(min(p) == p, 1), :);          % BestPosition为全局最优

% -------------- 进入主要循环，按照公式依次迭代，直到满足精度要求 --------------
w = [wmax zeros(1, MaxDT)];
PSOBestCurve = [min(p); zeros(MaxDT-1, 1)];         % 计算每次迭代最优解的收敛值
for t = 1 : MaxDT   
    for i = 1 : N                                   % BestPosition更新就不能并行
        v(i, :) = w(t) * v(i, :) + c1 * rand * (y(i, :) - x(i, :)) + ...
            c2 * rand * (BestPosition - x(i, :));
        for j = 1 : D
            v(i, j) = max(min(v(i, j), vmax(j)), vmin(j));
        end    
        x(i, :) = x(i, :) + v(i, :);
        x = [roundn(x(:, 1:end-1), -5) round(x(:, end))];	% 小数为由2到5
        for j = 1 : D
            x(i, j) = max(min(x(i, j), xmax(j)), xmin(j));
        end
        temp = Fitness(x(i, :), des);
        if temp < p(i)                              % 更新局部最优解 
            p(i) = temp;
            y(i, :) = x(i, :);
        end
        if p(i) < PSOBestCurve(t)                   % 更新全局最优解（不能用temp比较，因为第一次的存在）
            PSOBestCurve(t) = p(i);
            BestPosition = y(i, :);
        end
        waitbar((N+(t-1)*N+i) / irNum, fw, ['Loop process2 - ' num2str((t-1)*N+i) ' ... ' num2str(floor(toc)) ' s']);
    end
    w(t+1) = wmin + (MaxDT - t) * (wmax - wmin) / MaxDT;	% 惯性权重w从wmax逐渐下降到wmin
    if t < MaxDT
        PSOBestCurve(t+1) = PSOBestCurve(t);        % 初始化最优值
    end
    writePbest(PSOBestCurve(t), t, BestPosition);   % Record system every Pbest.
end

BestFitness = PSOBestCurve(end);

waitbar(1, fw, 'Iteration done !');
close(fw);
