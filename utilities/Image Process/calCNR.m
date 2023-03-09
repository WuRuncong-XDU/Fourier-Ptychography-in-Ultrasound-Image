% Calculate CNR for HFL40 simulation
function contrast = calCNR(input, mode)
if ~exist('mode', 'var')
    mode = 0;
end
if mode == 0
    iRegion = input(273:303, 376:406);
else
    iRegion = input;
end
[row, col] = find(iRegion == max(iRegion, [], 'all'), 1);
sidelobes = iRegion;
sidelobes(row-6:row+6, col-6:col+6) = 0;
% cal contrast value
contrast = 10 * log10(sum(sidelobes.^2, 'all') / sum(iRegion.^2, 'all'));   % less than 0
