function output = dynamic(input, dr, gain, mode, gamma)
if ~exist('mode', 'var')
    mode = 'sim';
end
input = (input - min(input, [], 'all')) / (max(input, [], 'all')-min(input, [], 'all')) * 2^16;
vbd = max(input, [], 'all');
maxdB = 20 * log10(vbd);
logValue = min(maxdB, max(0, 20 * log10(input)) + gain);
if strcmp(mode, 'sim')
    gamma = 4.5;
else
    if ~exist('gamma', 'var')
        gamma = 5;
    end
end
output = 255 * (max(0, logValue - maxdB + dr) / dr).^ gamma;
end
