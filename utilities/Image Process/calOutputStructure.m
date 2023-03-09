% For 30MHz image
function output = calOutputStructure(outputObject, spsize, depth)
%% 平均横向分辨率
interpTime = 20;
lineR = outputObject(288, 200:650);
lineR = interp1(1:length(lineR), lineR, linspace(1, length(lineR), length(lineR) * interpTime), 'pchip');
lineRD = diff(lineR);	% 计算导数
indexR = 0;
maxIndexR = zeros(1, 6);
for i = 2:length(lineRD)
    if lineRD(i-1) > 0 && lineRD(i) <= 0 && lineR(i) > 180
        indexR = indexR + 1;
        maxIndexR(indexR) = i;
    end
end
resR = 0;
for j = 1:length(maxIndexR)
    r = lineR(maxIndexR(j)-8*interpTime:maxIndexR(j)+8*interpTime);
    indexR = find(max(0, r - 0.5*max(r)) > 0);
    resR = resR + (indexR(end) - indexR(1)) * spsize / interpTime * 10^3;    % unit mm
end
output.resR = resR / length(maxIndexR); % 平均横向分辨率

%% 平均纵向分辨率（可以不计算）
lineH = outputObject(200:384, 391);
lineH = interp1(1:length(lineH), lineH, linspace(1, length(lineH), length(lineH) * interpTime), 'pchip');
lineHD = diff(lineH); % 计算导数
indexH = 0;
maxIndexH = zeros(1, 3);
for i = 2:length(lineHD)
    if lineHD(i-1) > 0 && lineHD(i) <= 0 && lineH(i) > 180
        indexH = indexH + 1;
        maxIndexH(indexH) = i;
    end
end
resH = 0;
pixelH = depth / size(outputObject, 1);  % unit [mm]
for j = 1:length(maxIndexH)
    h = lineH(maxIndexH(j)-8*interpTime:maxIndexH(j)+8*interpTime);
    indexH = find(max(0, h - 0.5*max(h)) > 0);
    resH = resH + (indexH(end) - indexH(1)) * pixelH / interpTime * 10^3;    % unit mm
end
output.resH = resH / length(maxIndexH); % 平均纵向分辨率

%% Contrast (ISLR)
iRegion = outputObject(273:303, 376:406);
[row, col] = find(iRegion == max(iRegion, [], 'all'), 1);
sidelobes = iRegion;
sidelobes(row-5:row+5, col-5:col+5) = 0;
% cal contrast value
output.contrast = 10 * log10(sum(sidelobes.^2, 'all') / sum(iRegion.^2, 'all'));   % less than 0

