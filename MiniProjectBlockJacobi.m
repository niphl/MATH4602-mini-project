function numIters = MiniProjectBlockJacobi(n, d, epsilon)
    numIters = 0;
    x = zeros(n^2,1);
    temp = zeros(n^2,1);
    % cholesky factorization of A_n
    ldiag = zeros(n,1); % diagonal elements
    llow = zeros(n,1); % elements below the diagonal (starts with dummy 0)
    ldiag(1) = sqrt(2^d);
    for i = 2:n
        llow(i) = -1/ldiag(i-1);
        ldiag(i) = sqrt(2^d-llow(i)^2);
    end
    % this b is valid only for n >= 3
    b = zeros(n^2,1);
    for i = 1:n^2
        b(i) = 2^d - 4 + (mod(i,n) == 1 || mod(i,n) == 0) + (i<=n || i > n^2 - n);
    end
    while norm(x-ones(n^2,1),2) > epsilon
        % computing b-(A-D)*x and storing it in x
        for i = 1:n^2
            if i <= n
                temp(i) = b(i) + x(i+n);
            elseif i > n^2 - n
                temp(i) = b(i) + x(i-n);
            else
                temp(i) = b(i) + x(i-n) + x(i+n);
            end
        end
        x = temp;
        % computing D \ x using the cholesky factorization computed earlier
        % forward substitution
        for i = 1:n
            temp(i*n-n+1) = x(i*n+1-n)/ldiag(1);
            for j = 2:n
                temp(i*n+j-n) = (x(i*n+j-n)-temp(i*n+j-n-1)*llow(j))/ldiag(j);
            end
        end
        x = temp;
        % backward substitution
        for i = 1:n
            temp(i*n) = x(i*n)/ldiag(n);
            for j = n-1:-1:1
                temp(i*n+j-n) = (x(i*n+j-n)-temp(i*n+j-n+1)*llow(j+1))/ldiag(j);
            end
        end
        x = temp;
        numIters = numIters + 1;
    end
end
