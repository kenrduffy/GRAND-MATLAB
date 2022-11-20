% Landslide algorithm for generating integer partitions from 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits 
% guessing random additive noise decoding," IEEE Trans. Signal 
% Process., vol. 70, pp. 4528-4542, 2022.


% W is the target logistic weight, w is the Hamming weight and n is the 
% length of the code.

function z = landslide(W,w,n)

    W1=W-w*(w+1)/2;
    n1=n-w;
    % Create the first integer partition
    jj=1;
    % Start with empty vector and breaking at first index
    u = zeros(1,w);
    k=1;
    u = mountain_build(u,k,w,W1,n1);
    z(jj,:)=u;
    % Evaluate drops
    d=circshift(u,-1)-u;
    d(w)=0;
    % Evaluate accumuated drops
    D = cumsum(d,'reverse');
    % Each loop generates a new integer partition
    while D(1)>=2
        % Find the last index with an accumulated drop >=2
        k=find(D>=2,1,'last');
        % Increase its index by one.
        u(k)=u(k)+1;
        u = mountain_build(u,k,w,W1,n1);
        % Record the partition
        jj=jj+1;
        z(jj,:)=u;
        % Evaluate drops
        d=circshift(u,-1)-u;
        d(w)=0;
        % Evaluate accumulated drops
        D = cumsum(d,'reverse');
    end
    z = z + repmat([1:w],size(z,1),1);
end

function u = mountain_build(u,k,w,W1,n1)
    u(k+1:w) = u(k)*ones(1,w-k);
    W2 = W1-sum(u);
    q = floor(W2/(n1-u(k)));
    r = W2-q*(n1-u(k));
    if q ~= 0
        u(w-q+1:w)=n1*ones(1,q);
    end
    if w-q>0
    	u(w-q)=u(w-q)+r;
    end
end
