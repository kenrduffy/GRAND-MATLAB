function [H, k] = H_eBCH(n,t)
if floor(log2(n)) ~= log2(n)
    error('n must be a power of two.');
end
m = log2(n);

PriPoly = primpoly(m);

nz = gf([0:2^m-1],m,PriPoly);
H = ones(n,1);
for i = 1:t
    temp = nz.^(2.*i-1);
    temp = double(temp.x)';
    temp = dec2bin(temp);
    H = [H str2arr(temp)];
end
H = H';

k = size(H,2) - size(H,1);
end

function out = str2arr(in)
[a,b] = size(in);
out = -1.*ones(a,b);
for i = 1:a
    for j = 1:b
        out(i,j) = str2num(in(i,j));
    end
end
end
