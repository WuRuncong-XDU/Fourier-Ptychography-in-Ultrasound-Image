function writePg(pg, Pbest)
filepath = '.\file\PgRecord.txt';
if isfile(filepath)
    fid = fopen(filepath, 'w');
else
    fid = fopen(filepath, 'a');
end
fprintf(fid, '%.3f  ', pg);
fprintf(fid, '\n%f\n\n', Pbest);
fclose(fid);
