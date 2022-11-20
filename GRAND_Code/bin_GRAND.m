%  GRAND (hard detection, BSC)

% Inputs:
%   n               - code length
%   H               - Parity check matrix or CRC function
%   max_query       - Maximum number of code-book queries to abandonment
%   y_demod         - Channel hard output
%
% Outputs:
%   y_decoded       - Decoded codeword
%   putative_noise  - noise
%   n_guesses       - Number of guesses performed
%   abandoned       - 1 if abandoned, 0 if found a codeword

function [y_decoded,putative_noise,n_guesses,abandoned] = bin_GRAND(H,max_query,y_demod)

    n=size(H,2);
    n_guesses = 0;
    abandoned = 0;
    
    [err_loc_vec, err_vec, ~] = gen_next_err(n); %Generate initial vectors

    while n_guesses < max_query
        [decoded,putative_noise,ng]=bin_syn_check(H,y_demod,err_vec);
        %How many guesses have been made
        n_guesses = n_guesses+ng;
        if (decoded == 1)
            y_decoded = mod(y_demod-putative_noise,2);
            return;
        end
        [err_loc_vec, err_vec, ~] = gen_next_err(n, err_loc_vec); %Generate next error vector
    end

    % If abandoned
    y_decoded = -1*ones(size(y_demod));
    abandoned = 1;


end

function [decoded,noise_found,n_guesses]=bin_syn_check(H,y_demod,err_vec)

    decoded = 0;
    noise_found = NaN(1,size(err_vec,2));
    % How many guesses have we made
    n_guesses=0;
    while (decoded == 0 && n_guesses<size(err_vec,1))
        % How much work have we done
        n_guesses=n_guesses+1;
        t = mod(y_demod-err_vec(n_guesses,:),2);
        if (mod(H*t',2) == zeros(size(H,1),1))
            decoded = 1;
            noise_found = err_vec(n_guesses,:);
        end   
    end

end

%
%Description: This function generates the new error location vector,
%given a previous error location vector and mask vector of possible error locations
%
%Inputs:
%   n           - Code length
%   err_loc_vec - Previous error location vector
%
%Outputs:
%   err_loc_vec - New error location vector. [] if a new error location vector cannot be generated
%   err_vec     - Binary error vector that corresponds to err_loc_vec. Zero vector if no error could be generated
%

function [err_loc_vec, err_vec, weight] = gen_next_err(n, err_loc_vec)


    mask_vec = 1:n;
    err_vec = zeros(1,n);

    if ~exist('err_loc_vec', 'var')
        %Initialize zero error location vector
        err_loc_vec = [];
        weight = 0;
        return;
    end


    %Get next error location vector
    err_loc_vec = increase_error(mask_vec, err_loc_vec);
    if isequal(err_loc_vec, []) %Could not generate error
        weight=inf;
        return;
    end

    %Generate error
    weight = length(err_loc_vec);
    err_vec(err_loc_vec) = 1;


    end

    function err_loc_vec = increase_error(mask_vec, err_loc_vec)
    %
    %Description: This function generates the next error location vector given the previous one
    %

    max_possible_err_loc = length(mask_vec);
    %Try to remain in the same weight
    success = false;
    for ii=length(err_loc_vec):-1:1
        new_err_loc_vec = err_loc_vec;
    %     success = true;
        if new_err_loc_vec(ii)==max(mask_vec) %Cannot increase this error
            continue;
        end
        success = true;
        ind_in_mask = find(mask_vec==new_err_loc_vec(ii)); %Find the index in mask_vec such that mask_vec(index)==err_loc_vec(index)
        new_err_loc_vec(ii) = mask_vec(ind_in_mask+1);
        for jj=ii+1:length(err_loc_vec)
            new_ind = ind_in_mask+jj-ii+1;
            if new_ind<=max_possible_err_loc
                new_err_loc_vec(jj) = mask_vec(new_ind);
            else
                success = false;
                break;
            end
        end
        if success==true
            err_loc_vec = new_err_loc_vec;
            break;
        end
    end

    %If need to go to the next weight
    if success==false
        if length(err_loc_vec)==max_possible_err_loc
            err_loc_vec = [];
            return;
        end
        for ii=1:length(err_loc_vec)+1
            err_loc_vec(ii) = mask_vec(ii);
        end
    end

end