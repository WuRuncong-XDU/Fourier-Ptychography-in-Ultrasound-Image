function [BestFitness, BestPosition, PSOBestCurve] = PSO(N, MaxIT, Dim, Des)
% ��ʼ��λ��
% xInitial = {[0.32 1.74 5], [0.227 1.973 4], [0.259 2.308 15], [0.292  1.507  15]};
xInitial = {};

% ��ʼ��PSO����
c1 = 2;
c2 = 2;                     % ѧϰ����
wmax = 0.9;                 % ����Ȩ��
wmin = 0.4;
D = Dim;                    % �����ռ�ά��
xmax = [1, 2.4, 15];        % [n0 ratio loop]
xmin = [0.1, 1.2, 3];       % ȷ��������Χ[0.1~1 1.2~2.4 2~10]
vmax = (xmax - xmin) / 2;   % vmax(3)=-6~6�е����
vmin = -vmax;               % �����ٶ�

% Main program of PSO
tic
[BestFitness, BestPosition, PSOBestCurve] = LDPSO(c1, c2, wmax, wmin, MaxIT, D, N, xmax, xmin, vmax, vmin, Des, xInitial);
toc

end
