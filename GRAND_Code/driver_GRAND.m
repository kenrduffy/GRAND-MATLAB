% Ken R. Duffy, 2018-2022.
% This drives the main simulation.

% Guessing Random Additive Noise Decoding (GRAND)
% All code is subject to license
% GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf

% GRAND algorithms efficiently decode codes where n-k is moderate. 
% For this inefficient MATLAB implementation, obtaining the full
% performance of codes with n-k>16 may be time-consuming.

% GRAND
% K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random 
% additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 
% 4023–4040, 2019.

% Ordered Reliability Bits GRAND (ORBGRAND)
% Basic ORBGRAND as introduced in
% K. R. Duffy, “Ordered reliability bits guessing random additive noise 
% decoding," in IEEE ICASSP, 2021, pp. 8268–8272 
% and implemented using the Landslide algorithm introduced in 
% ORBGRAND
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding,” IEEE Trans. Signal Process., vol. 70, 
% pp. 4528-4542, 2022.

% 1-line ORBGRAND from
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding,” IEEE Trans. Signal Process., vol. 70, 
% pp. 4528-4542, 2022.
% and implemented here using materials introduced in 
% K. Galligan, M. Médard, K. R. Duffy, "Block turbo decoding with ORBGRAND"
% arXiv:2207.11149, 2022.

clear;

% Chose the decoder from: 
% 'GRAND' (hard detection); 
% 'ORBGRAND' (soft detection); 
% 'ORBGRAND1' (soft detection);
DECODER='ORBGRAND1';

% Sim range in Eb/N0
ebn0=4:0.5:6;

% Run simulation until this many errors observed for each Eb/N0 value
err_thresh = 50;

% Modulation schemes available using MATLAB's toolbox, which will be used
% in a complex-valued channel
modlist = {'pi/2-BPSK','BPSK','QPSK','16QAM','64QAM','256QAM'}; 
bpsList = [1 1 2 4 6 8];
% Pick the modulation
modulation = 'BPSK';
% Determine the number of bits per symbol in the modulation
nmodbits = bpsList(strcmpi(modlist,modulation));

% Pick the code
% RLC, PAC, CAPOLAR, BCH or CRC. 
code_class = 'RLC';

% Random Linear Code. 
if isequal(code_class,'RLC')
    n=128;
    k=116;
    
    % If a RLC was previously used, reload it and reuse the code
    filename = ['../RESULTS/' DECODER '_' code_class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
    if exist(filename, 'file') == 2
        load(filename,'code');
        G=code.G;
        H=code.H;
    % If not, make a random generator
    else
        % Make a random parity check matrix
        [G,H] = make_RLC(k,n,0.5);
        code.G=G;
        code.H=H;
    end

% Polar-assisted convolutional codes, as introduced by 
% E. Arikan, "From sequential decoding to channel polarization 
% and back again", arXiv:1908.09594, 2019.
elseif isequal(code_class,'PAC')
    n=128; % Must be power of 2.
    k=116;
    % Generator for convolutional code, taken from arXiv:1908.09594.
    conv_code = [1 1 1 0 1 1 0 1];
    [G, H] = make_pac_code(n,k,conv_code);
    code.G = G;
    code.H = H;
    % Record what convolution was used.
    code.conv_code = conv_code;

% CRC-Assisted Polar code. As GRAND algorithms are universal, rather than 
% use the Polar bits for error correction and the CRC bits for error 
% detection, as the product of linear codes is linear, instead GRAND uses 
% all bits for error correction.
elseif isequal(code_class,'CAPOLAR')
    % 5G New Radio uplink code
    filename = '../CODES/capolar_k_116_n_128_ul.mat';
    code = open(filename);
    G = code.G;
    H = code.H;
    % k,n is the code length
    [k,n]=size(G);

% Bose–Chaudhuri–Hocquenghem code. Normally only used with hard 
% detection decoding, but soft version of GRAND can decode too.
elseif isequal(code_class,'BCH')
    filename = '../CODES/BCH_k_113_n_127.mat';
    code = open(filename);
    code.class=code_class;
    G = code.G;
    H = code.H;
    % k,n is the code length
    [k,n]=size(G); 

% Cyclic Redundancy Check (CRC) code. Normally only used for error
% detection, but GRAND algorithms can decode with hard or soft information.
elseif isequal(code_class,'CRC')
    k=116;
    % Polynomial from https://users.ece.cmu.edu/~koopman/crc/ which has a
    % repository of high quality CRCs curated and maintained by 
    % Philip Koopman, Carnegie Mellon University.
    hex_poly = '0x8f3'; % 12 bit
    poly=koopman2matlab(hex_poly); % Convert from Koopman notation
    % Record the polynomial
    code.poly = poly;
    % Set up the CRC check using MATLAB toolbox
    [G,H,n] = make_CRC_GH(k,poly);
    code.G=G;
    code.H=H;
end

% Reminder that this implementation is not optimized or parallelized, so 
% larger n-k may take time to run.
if n-k>16
    disp('GRAND and ORBGRAND are readily highly parallelizable in software and hardware.')
    disp('The MATLAB implementations here are for instruction and are NOT parallelized.')
    disp(['If GRAND algorithms find an error, they do so after ' ...
        'approximately 2^{n-k} code-book queries'])
    disp('and so 2^{n-k} is an upper bound on complexity in terms of query numbers.')
    disp(['With n-k=' num2str(n-k) ' expect slow decoding in noisy channels with this inefficient implementation.'])
end

% Record the code_class in the data construct
code.class = code_class;

%GRAND parameters
max_query=inf; % Can reduce to an abandonment value.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert Eb/N0 to SNR
snr_db = ebn0+10*log10(k/n)+10*log10(nmodbits);

% Pad modulation with meaningless bits if necessary.
mod_pad=zeros(1,mod(n,nmodbits));

%For each SNR record the following
num_decoded = zeros(1,length(snr_db)); % Number of packets decoded
num_demod_errs = num_decoded; % Number of demodulations that are in error
num_demod_bit_errs = num_decoded; % Number of erroneous demodulated bits
num_errs = num_decoded; % Number of erroneous decodings
num_bit_errs = num_decoded; % Number of erroneous bits
num_aband = num_decoded; % Number of abandoned decodings
num_queries = num_decoded; % Total number of code-book queries

for ii=1:length(snr_db)
    
    % Noise variance
    sigma2 = 1/(10^(0.1*snr_db(ii))); 
    % Using MATLAB's channel function
    awgnchan = comm.AWGNChannel('NoiseMethod','Variance','Variance',sigma2);

    % Keep decoding pacekts until you see err_thresh erroneous decodings
    while num_errs(ii)<err_thresh
        num_decoded(ii)=num_decoded(ii)+1;
                
        % Encode
        % Uniform random information word
        u = binornd(ones(1,k),0.5);
        % Codeword
        c = mod(u*G,2);
	    % Modulate
	    modOut = nrSymbolModulate([c mod_pad]',modulation);
        % Add White Gaussian noise
        rSig = awgnchan(modOut);
        % Soft demodulate
        y_soft = nrSymbolDemodulate(rSig,modulation,sigma2);
    	y_soft = y_soft(1:n)';
        % Hard demodulate
        y_demod = (y_soft<0);

        % Count the times the demodulation is in error
        if ~isequal(y_demod, c)
            num_demod_errs(ii) = num_demod_errs(ii)+1;
            % Count the demodulated bit errors
            num_demod_bit_errs(ii) = num_demod_bit_errs(ii)+sum(abs(y_demod-c));
        end
        
        % Decode with GRAND (hard detection)
        if isequal(DECODER,'GRAND')
            [y_decoded,~,n_guess,abandoned] = bin_GRAND(H,max_query,y_demod); 
        % Decode with basic ORBGRAND (soft detection)
        elseif isequal(DECODER,'ORBGRAND')
            [y_decoded,~,n_guess,abandoned] = bin_ORBGRAND(H,max_query,y_soft); 
        % Decode with 1-line ORBGRAND (soft detection)
        elseif isequal(DECODER,'ORBGRAND1')
            [y_decoded,~,n_guess,abandoned] = bin_ORBGRAND1(H,max_query,y_soft); 
        else
            disp('DECODER must be GRAND or ORBGRAND or ORBGRAND1')
            return;
        end

        % Total number of queries made at this SNR
        num_queries(ii) = num_queries(ii) + n_guess;

        % If there is an error in the decoding
        if ~isequal(y_decoded, c)
            % Increment the number of errors observed at this SNR
            num_errs(ii)=num_errs(ii)+1;
            % Count the bit errors.
            if ~abandoned
                num_bit_errs(ii) = num_bit_errs(ii) + sum(abs(y_decoded-c));
            else
                % If we've abandoned, then use the demodulated sequence 
                % to determine how many bit errors there were.
                num_bit_errs(ii) = num_bit_errs(ii) + sum(abs(y_demod-c));
                num_aband(ii) = num_aband(ii)+1;
            end
            % Report to the terminal every time a new error is observed
            disp(['n=' num2str(n) ', k=' num2str(k) ', R=' num2str(k/n,'%.2f') ', ' num2str(ebn0(ii)) ' dB, Dec=' num2str(num_decoded(ii)) ', Err=' num2str(num_errs(ii)) ', BLER=' num2str(num_errs(ii)/num_decoded(ii)) ', BER=' num2str(num_bit_errs(ii)/(num_decoded(ii)*n)) ', EG=' num2str(num_queries(ii)/num_decoded(ii),'%.0f')]); 
        end
    end
end


% Save the results.

if isequal(code.class,'CRC')
    filename = ['../RESULTS/' DECODER '_' code.class '_' hex_poly '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
elseif isequal(code.class,'PAC')
    filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(bin2dec(num2str(conv_code))) '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
else
    filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_' num2str(nmodbits) '.mat'];
end

% Store the data in a way where results from future simulations with the 
% same code can be appended.

if exist(filename, 'file') == 2
    load(filename,'code')
    for ii=1:length(snr_db)
        these = find(code.snr==snr_db(ii));
        if isempty(these)
            code.snr(end+1) = snr_db(ii);
            code.num_decoded(end+1)=num_decoded(ii);
            code.num_demod_errs(end+1)=num_demod_errs(ii);
            code.num_demod_bit_errs(end+1)=num_demod_bit_errs(ii);
            code.num_errs(end+1)=num_errs(ii);
            code.num_bit_errs(end+1)=num_bit_errs(ii);
            code.num_aband(end+1)=num_aband(ii);
            code.num_queries(end+1)=num_queries(ii);
        else
            code.num_decoded(these)=code.num_decoded(these)+num_decoded(ii);
            code.num_demod_errs(these)=code.num_demod_errs(these)+num_demod_errs(ii);
            code.num_demod_bit_errs(these)=code.num_demod_bit_errs(these)+num_demod_bit_errs(ii);
            code.num_errs(these)=code.num_errs(these)+num_errs(ii);
            code.num_bit_errs(these)=code.num_bit_errs(these)+num_bit_errs(ii);
            code.num_aband(these)=code.num_aband(these)+num_aband(ii);
            code.num_queries(these)=code.num_queries(these)+num_queries(ii);
        end
    end     
else
    code.n = n;
    code.k = k;
    code.modulation = modulation;
    code.nmodbits = nmodbits;
    code.G = G;
    code.snr = snr_db;
    code.num_decoded=num_decoded;
    code.num_demod_errs=num_demod_errs;
    code.num_demod_bit_errs=num_demod_bit_errs;
    code.num_errs=num_errs;
    code.num_bit_errs=num_bit_errs;
    code.num_aband=num_aband;
    code.num_queries=num_queries;
end

% Sort according to increasing snr;
[~,ord]=sort(code.snr);
code.snr = code.snr(ord);
code.num_decoded = code.num_decoded(ord);
code.num_demod_errs = code.num_demod_errs(ord);
code.num_demod_bit_errs = code.num_demod_bit_errs(ord);
code.num_errs = code.num_errs(ord);
code.num_bit_errs = code.num_bit_errs(ord);
code.num_aband = code.num_aband(ord);
code.num_queries = code.num_queries(ord);

% Summary statistics   
code.ebn0 = code.snr-10*log10(k/n)-10*log10(nmodbits); %Eb/N0
code.R = code.k/code.n; % Code rate
code.BLERdemod = code.num_demod_errs./code.num_decoded; % Demod BLER
code.BERdemod = code.num_demod_bit_errs./(code.num_decoded*code.n); % Demod BER
code.BLER = code.num_errs./code.num_decoded; % Decoded BLER
code.BER = code.num_bit_errs./(code.num_decoded*code.n); % Decoded BER
code.EG = code.num_queries./code.num_decoded; % Average number of code-book queries
code.max_query = max_query;

% Code
code.G=G;
code.H=H;

save(filename,'code')

disp(['Decoder=' DECODER ', code=' code.class ', ' num2str(sum(num_decoded)) ' packets decoded.'])


