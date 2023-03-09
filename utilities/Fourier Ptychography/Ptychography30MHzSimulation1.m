%% Ptychography for HFL40 30 MHz transducer, Verasonics standard simulation data.
function [Output, OutputImage, RefImage] = Ptychography30MHzSimulation1(inputParams)
%% 输入参数
n0 = inputParams(1);        % Replace refraction
ratio = inputParams(2);     % 优化参数2, 1.2~2.4
loop = inputParams(3);      % 优化参数3

%% setup the parameters for the coherent imaging system
global data;
imSeqLowRes = data;
[depthPts, radiusPts, number] = size(data);

% HFL40 : 31.25MHz, 1540m/s
spacing = 0.069 / 1000;         % Pitch in [m]
waveLength = 1540 / 31.25e6;    % 4.928e-05
k0 = 2 * pi / waveLength;
spsize = spacing * 256 / radiusPts;     % 根据图片像素决定
psize = spsize / ratio;         % Final pixel size of the reconstrucio1 （必须比1大）

depthMm = 8e-3;                 % 深度8mm
pixelH = depthMm / depthPts;	% 纵向分辨率
startP = 5 * waveLength;        % lambda
h = (0 : depthPts-1) * pixelH + startP;
% h = (1 : depthPts) / 4 * waveLength;  % image data have been calibration not match 4* acq-rate
r = spsize * radiusPts / 2;
NA = r ./ sqrt(r^2 + h.^2);

%% 将声源比作阵列式的光源
theta = [-4 -3 -2 -1 0 1 2 3 4] * ((25*pi/180) / (9-1));	% 直接得到各个角度
kx_relative = -sin(theta);	% Create kx, ky wavevectors

n = floor(radiusPts * (spsize / psize));	% Image size of the final output
kx = k0 * kx_relative;
dkx = 2 * pi / (psize * n);
cutoffFrequency = n0 * NA * k0;             % correcte cutoff frequency
kmax = pi / spsize;
kxm = linspace(-kmax, kmax, radiusPts);

%% recover the high resolution image
objectRecover = ones(size(imSeqLowRes, 1), n);	% initial guess of the object
% fftshift对二维矩阵一起转换，所以只能分层进行处理
objectRecoverFT = zeros(size(objectRecover));

parallelMode = 'on';
if strcmp(parallelMode, 'off')
    seq = gseq(number);                 % 从零开始，左右（负正角度）交替的顺序
    for rowIndex = 1:size(imSeqLowRes, 1)
        objectRecoverFT(rowIndex, :) = fftshift(fft(objectRecover(rowIndex, :)));
        CTF = (kxm.^2 < cutoffFrequency(rowIndex)^2);
        for tt = 1 : loop
            for i3 = 1 : number
                i2 = seq(i3);
                kxc = round((n+1)/2 + kx(1, i2)/dkx);
                kxl = round(kxc - (radiusPts-1)/2);
                kxh = round(kxc + (radiusPts-1)/2);
                lowResFT = (radiusPts/n)^2 * objectRecoverFT(rowIndex, kxl:kxh) .* CTF;
                im_lowRes = ifft(ifftshift(lowResFT));
                im_lowRes = (n/radiusPts)^2 * imSeqLowRes(rowIndex, :, i2) .* exp(1i .* angle(im_lowRes));
                lowResFT = fftshift(fft(im_lowRes)) .* CTF;
                objectRecoverFT(rowIndex, kxl:kxh) = (1-CTF) .* objectRecoverFT(rowIndex, kxl:kxh) + lowResFT;
            end
        end
        objectRecover(rowIndex, :) = ifft(ifftshift(objectRecoverFT(rowIndex, :)));
    end
else
    parfor rowIndex = 1 : size(imSeqLowRes, 1)
        tmpobjectRecoverFT = fftshift(fft(objectRecover(rowIndex, :)));
        CTF = (kxm.^2 < cutoffFrequency(rowIndex)^2);
        tmpseq = gseq(number);
        tmpkx = kx;
        tmpimSeqLowRes = imSeqLowRes;
        for tt = 1 : loop
            for i3 = 1 : number
                i2 = tmpseq(i3);
                kxc = round((n+1)/2 + tmpkx(1, i2)/dkx);
                kxl = round(kxc - (radiusPts-1)/2);
                kxh = round(kxc + (radiusPts-1)/2);
                lowResFT = (radiusPts/n)^2 * tmpobjectRecoverFT(kxl:kxh) .* CTF;
                im_lowRes = ifft(ifftshift(lowResFT));
                im_lowRes = (n/radiusPts)^2 * tmpimSeqLowRes(rowIndex, :, i2) .* exp(1i .* angle(im_lowRes));
                lowResFT = fftshift(fft(im_lowRes)) .* CTF;
                tmpobjectRecoverFT(kxl:kxh) = (1-CTF) .* tmpobjectRecoverFT(kxl:kxh) + lowResFT;
            end
        end
        objectRecover(rowIndex, :) = ifft(ifftshift(tmpobjectRecoverFT(:)));
    end
end

%% 图像处理 & 校正
% Output Image
OutputImage = dynamic(abs(objectRecover), 45, 2.5);
OutputImage = max(0, imresize(OutputImage, [depthPts, radiusPts]));

% Reference Image
RefImage = imSeqLowRes(:, :, 5);
RefImage = dynamic(abs(RefImage), 45, 2.5);
RefImage = max(0, imresize(RefImage, [depthPts, radiusPts]));

%% Cal output structure
Output = calOutputStructure(OutputImage, spsize, depthMm);
