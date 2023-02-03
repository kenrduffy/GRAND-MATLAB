function [chat_list, s_list, NT] = ORBGRAND(llr, H, L, Tmax)
n = length(llr);
chat_list = []; 
s_list = [];
NT = 1;
[~, il] = sort(abs(llr));
[~, invil] = sort(il);

c_HD = HardDec(llr);

if sum(mod(H * c_HD,2)) == 0
    chat_list = [chat_list, c_HD]; 
    s_list = [s_list, 0];
end

TEP = zeros(n,1);
while (NT < Tmax && (length(s_list)<L) )
    TEP = NextTEP(TEP, n);
    chat = mod(c_HD+TEP(invil),2);
    NT = NT + 1;
    if sum(mod(H * chat,2)) == 0
        chat_list = [chat_list, chat];
        s_list = [s_list, findLW(TEP, n)];
    end
end
end

function cHD = HardDec(L)
cHD = zeros(length(L),1);
for i = 1:length(L)
    if L(i) < 0
        cHD(i) = 1;
    end
end
end


function LW = findLW(TEP, n)
LW = 0;
for i = 1:n
   if TEP(i) == 1
       LW = LW + i;
   end
end
end

function  flag = isLastTEP(TEP, n)
LW = findLW(TEP, n);
if LW == 0
    flag = 1;
    return
end
for i = n:-1:1
    if i*(i-1) < 2*LW
        if TEP(i) == 0
            flag = 0;
            return
        end
        LW = LW - i;
        if LW == 0
            flag = 1;
            return
        end
    end
end
flag = 0;
end

function  TEP = findLastTEP(n, LW)
TEP = zeros(n,1);
for i = n:-1:1
    if i <= LW
        TEP(i) = 1;
        LW = LW - i;
        if LW == 0
            return;
        end
    end
end
end

function  EP = NextTEP(TEP, n)
LW = findLW(TEP, n);
if isLastTEP(TEP, n) == 1
    EP = findLastTEP(n, LW+1);
    return
else
    while(1)
        lm = find(TEP(3:n) == 1,1) + 2;
        r = findLW([zeros(lm,1); TEP(lm+1:n)], n);
        if lm*(lm-1)/2 + r < LW
            TEP(1:lm) = 0;
        else
            EP = [findLastTEP(lm-1, LW-r); 0; TEP(lm+1:n)];
            return
        end
    end
end
end