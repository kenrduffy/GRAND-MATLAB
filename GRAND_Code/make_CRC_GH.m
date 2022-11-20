% Convert a generator polynomial into [G,H] binary code matrices
% Requires MATLAB's comms toolbox


function [G,H,n] = make_CRC_GH(k,poly)

    GCRC = comm.CRCGenerator(poly);
    n = length(GCRC(zeros(k,1)));
    
    u=[zeros(1,k-1) 1];
    
    G = zeros(k,n);
    
    for ii=1:k
        u=circshift(u,1);
        G(ii,:) = (GCRC(u')'==1);
    end
    
    H =[G(:,k+1:end)' eye(n-k)];
    
end
