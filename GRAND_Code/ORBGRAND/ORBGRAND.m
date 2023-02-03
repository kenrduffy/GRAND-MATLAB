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

ErrorPattern = zeros(n,1);
while (NT < Tmax && (length(s_list)<L) )
    ErrorPattern = NextErrorPattern(ErrorPattern, n);
    chat = mod(c_HD+ErrorPattern(invil),2);
    NT = NT + 1;
    if sum(mod(H * chat,2)) == 0
        chat_list = [chat_list, chat];
        s_list = [s_list, LogisticWeight(ErrorPattern, n)];
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

function chat = AddTEP(cHD,pattern)
chat = zeros(length(cHD),1);
for i = 1:length(cHD)
    if pattern(i) == true
        chat(i) = 1 - cHD(i);
    else
        chat(i) = cHD(i);
    end
end
end

function LW = LogisticWeight(ErrorPattern, n)
LW = 0;
for i = 1:n
   if ErrorPattern(i) == 1
       LW = LW + i;
   end
end
end

function  flag = IsLast(ErrorPattern, n)
LW = LogisticWeight(ErrorPattern, n);
if LW == 0
    flag = 1;
    return
end
for i = n:-1:1
    if i*(i-1) < 2*LW
        if ErrorPattern(i) == 0
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
end

function  ErrorPattern = MaxIntegerPartition(n, LW)
ErrorPattern = zeros(n,1);
for i = n:-1:1
    if i <= LW
        ErrorPattern(i) = 1;
        LW = LW - i;
        if LW == 0
            return;
        end
    end
end
end

function  EP = NextErrorPattern(ErrorPattern, n)
LW = LogisticWeight(ErrorPattern, n);
if IsLast(ErrorPattern, n) == 1
    EP = MaxIntegerPartition(n, LW+1);
    return
else
    while(1)
        lm = find(ErrorPattern(3:n) == 1,1) + 2;
        %r = LogisticWeight(ErrorPattern(lm+1:n), n-lm);
        r = LogisticWeight([zeros(lm,1); ErrorPattern(lm+1:n)], n);
        if lm*(lm-1)/2 + r < LW
            ErrorPattern(1:lm) = 0;
        else
            EP = [MaxIntegerPartition(lm-1, LW-r); 0; ErrorPattern(lm+1:n)];
            %EP = [MaxIntegerPartition(LW-r, lm-1), 0, ErrorPattern(lm+1:n)];
            return
        end
    end
end
end