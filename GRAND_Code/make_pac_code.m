% Create a Polar-assisted convolutional code, following the instruction
% in E. Arikan, "From sequential decoding to channel polarization 
% and back again", arXiv:1908.09594, 2019.

function [G, H] = make_pac_code(N, K, c)

    %Inputs: codeword length as integer N, 
    %        number of data bits as integer K,
    %        convolution polynomial c as array of ints, 

    %Outputs: K by N binary generator matrix G
    %         N-K by N binary parity check matrix H

    %generate rate-profiling vector A using score metric

    [~, indices] = sort(sum(de2bi([0:(N-1)]),2), 'descend');
    A = sort(indices(1:K)');

    %rate profiling matrix given vector A of indices
    rate_prof_mat = zeros(K,N);
    for i = 1:K
        rate_prof_mat(i, A(1,i)) = 1;
    end

    %Generate convolutional precoding matrix given polynomial c, 

    %T is convolution matrix
    T = zeros(N,N);
    for i = 1:length(c)
        k = 0;
        for j = 1:(N+1-i)
            T(j,i+k) = c(1,i);
            k = k+1;
        end
    end

    %polar coding matrix, the Kronecker power of polar transform to
    %power of log2(N)
    P_1 = [1 0; 1 1];
    P_n = P_1;
    for i = 1:(log2(N)-1)
        P_n = kron(P_n, P_1);
    end

    %generator matrix G is product of these matrices, modulo 2
    G = mod(rate_prof_mat*T*P_n,2);
    
    %need to get G in standard form through row ops of modulo 2
    stand_G = G;
    i = 1;
    j = 1;

    %for each column find non-zero element, swap rows, subtract from other
    %rows
    while (i<=K)&&(j<=N)
        row_offset = find(stand_G(i:K,j),1);

        if isempty(row_offset)
            j=j+1;
        else
            pivot_row = row_offset + i - 1;

            %swap rows
            stand_G([i pivot_row],j:N) = stand_G([pivot_row i],j:N);
            for ii = [1:i-1 i+1:K]
                stand_G(ii,j:N) = mod(stand_G(ii,j:N) - stand_G(ii,j)*stand_G(i,j:N),2);
            end
            i = i+1;
            j = j+1;
        end
    end
    
    %only valid if G is full rank
    if rank(stand_G) == K

        %if not in standard form still, swap columns and find a column
        %swapped parity check matrix, and switch columns back to get H
        if ~isequal(stand_G(:,1:K),eye(K))
            swap_cols = [1:N];
            col_swap_G = stand_G;
            for i = 1:K
                if col_swap_G(i,i) == 0
                    for j = i+1:N
                        if col_swap_G(i,j) == 1
                            col_swap_G(:,[i j]) = col_swap_G(:, [j i]);
                            swap_cols(1, [i j]) = swap_cols(1, [j i]);
                            break
                        end
                    end
                end
            end

            %make makeshift H then swaps cols back to get our true H
            col_swap_H = mod(gen2par(col_swap_G),2);
            H(:, swap_cols) = col_swap_H;
           
        %already in standard form
        else
            H = mod(gen2par(stand_G),2);
        end
    
    else
        error('Needs full rank G - try different c array.')
    end
end
