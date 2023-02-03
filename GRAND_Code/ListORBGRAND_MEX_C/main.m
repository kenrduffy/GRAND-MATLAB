clear;clc;

%% define H of eBCH
n = 64;
t = 2;

[H, k] = H_eBCH(n,t);

L = 2;
%% Monte Carlo
EbN0dB = 2:0.25:5;
EsN0dB = EbN0dB + 10*log10(2*k/n);
NoErrors = 100;
maxIt = 10^6;

Tmax = 10^6;

%% all zero codeword
BLER = zeros(1,length(EsN0dB));
AvgN = zeros(1,length(EsN0dB));
for sp=1:length(EsN0dB)
    BlockError = 0;
    sumN = 0;   
    ntx = 0;
    scal = sqrt(10^(EsN0dB(sp)/10));

    t_c = 0;
    t_matlab = 0;

    while (BlockError < NoErrors && ntx < maxIt)
        ntx = ntx+1;
        
        c = zeros(n,1);
        x = (1-2*c)*scal;
        y = x + randn([n 1]);
        llr = 2*scal*y;
        
        % c implementation
        tic
        [chat, score, NT] = orbGRAND(llr, uint8(reshape(H',[],1)), uint64(L), uint64(Tmax));
        t_c = t_c + toc;

        % matlab implementation
        tic
        [chat1, score1, NT1] = ORBGRAND(llr, H, L, Tmax);
        t_matlab = t_matlab + toc;

        sumN = sumN + NT;

        if (~isequal(NT,NT1))
            disp("error: matlab and c implementation doesnt match")
            
        end
        
        if (~isequal(c,chat(:,1)))
            BlockError = BlockError + 1;
        end
    end
    AvgN(sp) = sumN/ntx;
    BLER(sp) = BlockError/ntx;

    disp(['C/MEX  implementation: ' num2str(t_c/ntx) ' seconds']);
    disp(['Matlab implementation: ' num2str(t_matlab/ntx) ' seconds']);
    
    disp([num2str(EbN0dB(sp)) ' dB: BLER        = ' num2str(BLER(sp))]);    
    disp([num2str(EbN0dB(sp)) ' dB: AvgN        = ' num2str(AvgN(sp))]);
    %save(['results/Polar-' num2str(n) '-' num2str(k) '-CRC-' num2str(l_crc) '-' hexPoly '-L-' num2str(L) '.mat'],'EsN0dB','FER','FER_ML');
end
