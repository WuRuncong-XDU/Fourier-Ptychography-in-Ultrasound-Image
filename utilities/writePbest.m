function writePbest(Pbest, t, pg)
filepath = '.\file\PbestRecord.txt';
if isfile('.\file\PbestRecord.txt')
    fid = fopen(filepath, 'w');
else
    fid = fopen(filepath, 'a');
end
fprintf(fid, '%d\n', t);
fprintf(fid, '%.3f  ', pg);
fprintf(fid, '\n%f\n\n', Pbest);
fclose(fid);
