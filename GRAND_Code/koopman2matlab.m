% Convert a polynomial in Koopman notation into one suitable for MATLAB's
% comms package
function poly = koopman2matlab(k_poly)

    poly= dec2bin(hex2dec(k_poly))-'0';
    poly = [poly 1];
    
end
