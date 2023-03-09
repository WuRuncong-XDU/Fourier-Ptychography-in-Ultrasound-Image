function [BestFitness, BestPosition, PSOBestCurve] = PSO(N, MaxIT, Dim, Des)
% 初始化位置
% xInitial = {[0.32 1.74 5], [0.227 1.973 4], [0.259 2.308 15], [0.292  1.507  15]};
xInitial = {};

% 初始化PSO参数
c1 = 2;
c2 = 2;                     % 学习因子
wmax = 0.9;                 % 惯性权重
wmin = 0.4;
D = Dim;                    % 搜索空间维数
xmax = [1, 2.4, 15];        % [n0 ratio loop]
xmin = [0.1, 1.2, 3];       % 确定搜索范围[0.1~1 1.2~2.4 2~10]
vmax = (xmax - xmin) / 2;   % vmax(3)=-6~6有点高了
vmin = -vmax;               % 搜索速度

% Main program of PSO
tic
[BestFitness, BestPosition, PSOBestCurve] = LDPSO(c1, c2, wmax, wmin, MaxIT, D, N, xmax, xmin, vmax, vmin, Des, xInitial);
toc

end
