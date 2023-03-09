function seq = gseq(number)
seq = zeros(1, number);
seq(1) = round(number / 2);
seq(2) = seq(1) - 1;
for i = 3:number
    if mod(i, 2) == 1
        seq(i) = seq(i-2) + 1;
    else
        seq(i) = seq(i-2) - 1;
    end
end
end