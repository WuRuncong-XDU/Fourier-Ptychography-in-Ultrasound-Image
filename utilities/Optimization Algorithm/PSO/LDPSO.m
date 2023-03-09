function [BestFitness, BestPosition, PSOBestCurve] = LDPSO(c1, c2, wmax, wmin, MaxDT, D, N, xmax, xmin, vmax, vmin, des, xInitial)
% ��������
irNum = N + N * MaxDT;
fw = waitbar(0, 'Process ... ');
tic

% c1ѧϰ����1; c2ѧϰ����2; w����Ȩ��; MaxDT����������
% D�����ռ�ά����δ֪������; N��ʼ��Ⱥ�������Ŀ;��

% -------------- ��ʼ����Ⱥ�ĸ���(�����������޶�λ�ú��ٶȵķ�Χ) --------------
x = zeros(N, D);
v = zeros(N, D);
for i = 1 : N
    for j = 1 : D
       x(i, j) = rand * (xmax(j) - xmin(j)) + xmin(j);	% �����ʼ��λ��
    end
    v(i, 1:D) = rand * (vmax - vmin) + vmin;            % �����ʼ���ٶ�
end
x(:, 1:end-1) = roundn(x(:, 1:end-1), -2);
x(:, end) = round(x(:, end));
if exist('xInitial', 'var')
    for i = 1 : length(xInitial)
        x(i, :) = xInitial{i};
    end
end

% -------------- �ȼ���������ӵ���Ӧ�ȣ�����ʼ��Pi��Pg --------------
p = zeros(1, N);    % �ֲ����Ž⣨�����
y = zeros(size(x)); % �ֲ����Ž⣨������
for i = 1 : N
    p(i) = Fitness(x(i, :), des);
    y(i,: ) = x(i, :);
    waitbar(i / irNum, fw, ['Loop process1 - ' num2str(i) ' ... ' num2str(floor(toc)) ' s']);
end
BestPosition = x(find(min(p) == p, 1), :);          % BestPositionΪȫ������

% -------------- ������Ҫѭ�������չ�ʽ���ε�����ֱ�����㾫��Ҫ�� --------------
w = [wmax zeros(1, MaxDT)];
PSOBestCurve = [min(p); zeros(MaxDT-1, 1)];         % ����ÿ�ε������Ž������ֵ
for t = 1 : MaxDT   
    for i = 1 : N                                   % BestPosition���¾Ͳ��ܲ���
        v(i, :) = w(t) * v(i, :) + c1 * rand * (y(i, :) - x(i, :)) + ...
            c2 * rand * (BestPosition - x(i, :));
        for j = 1 : D
            v(i, j) = max(min(v(i, j), vmax(j)), vmin(j));
        end    
        x(i, :) = x(i, :) + v(i, :);
        x = [roundn(x(:, 1:end-1), -5) round(x(:, end))];	% С��Ϊ��2��5
        for j = 1 : D
            x(i, j) = max(min(x(i, j), xmax(j)), xmin(j));
        end
        temp = Fitness(x(i, :), des);
        if temp < p(i)                              % ���¾ֲ����Ž� 
            p(i) = temp;
            y(i, :) = x(i, :);
        end
        if p(i) < PSOBestCurve(t)                   % ����ȫ�����Ž⣨������temp�Ƚϣ���Ϊ��һ�εĴ��ڣ�
            PSOBestCurve(t) = p(i);
            BestPosition = y(i, :);
        end
        waitbar((N+(t-1)*N+i) / irNum, fw, ['Loop process2 - ' num2str((t-1)*N+i) ' ... ' num2str(floor(toc)) ' s']);
    end
    w(t+1) = wmin + (MaxDT - t) * (wmax - wmin) / MaxDT;	% ����Ȩ��w��wmax���½���wmin
    if t < MaxDT
        PSOBestCurve(t+1) = PSOBestCurve(t);        % ��ʼ������ֵ
    end
    writePbest(PSOBestCurve(t), t, BestPosition);   % Record system every Pbest.
end

BestFitness = PSOBestCurve(end);

waitbar(1, fw, 'Iteration done !');
close(fw);
