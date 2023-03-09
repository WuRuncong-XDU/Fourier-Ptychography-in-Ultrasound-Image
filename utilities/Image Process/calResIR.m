function output = calResIR(lineData, delta, limit, range)
range = round(range / 2);
interpTime = 20;
lineData = interp1(1:length(lineData), lineData, linspace(1, length(lineData), length(lineData) * interpTime), 'pchip');
lineRD = diff(lineData);	% 计算导数
indexR = 0;
maxIndexR = zeros(1, 1);
for i = 2:length(lineRD)
    if lineRD(i-1) > 0 && lineRD(i) <= 0 && lineData(i) > limit
        indexR = indexR + 1;
        maxIndexR(indexR) = i;
    end
end
res = 0;
try
    for j = 1:length(maxIndexR)
        r = lineData(maxIndexR(j)-range*interpTime:maxIndexR(j)+range*interpTime);
        indexR = find(max(0, r - 0.5*max(r)) > 0);
        res = res + (indexR(end) - indexR(1)) * delta / interpTime * 10^3;    % unit mm
    end
    output = res / length(maxIndexR); % 平均分辨率
catch message
    % figure;
    % plot(lineData);
    % fprintf(['calResIR.m - ' message.message, '\n']);
    output = 1;
    return
end
end