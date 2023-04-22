function numIters = MiniProjectGaussSeidel(n, d, epsilon)
    numIters = 0;
    x = zeros(n^2,1);
    temp = zeros(n^2,1);
    % this b is valid only for n >= 3
    b = zeros(n^2,1);
    for i = 1:n^2
        b(i) = 2^d - 4 + (mod(i,n) == 1 || mod(i,n) == 0) + (i<=n || i > n^2 - n);
    end
    while norm(x-ones(n^2,1),2) > epsilon
        % computing b-U*x and saving it as x
        for i = 1:n^2
            if i == n^2
                temp(i) = b(i);
            elseif i > n^2 - n
                temp(i) = b(i) + x(i+1);
            else
                temp(i) = b(i) + x(i+1)*(mod(i,n)~=0) + x(i+n);
            end
        end
        x = temp;
        % computing L \ x
        for i = 1:n^2
            if i == 1
                temp(i) = x(i)/2^d;
            elseif i <= n
                temp(i) = (x(i)+temp(i-1))/2^d;
            else
                temp(i) = (x(i)+temp(i-1)*(mod(i,n)~=1)+temp(i-n))/2^d;
            end
        end
        x = temp;
        numIters = numIters + 1;
    end
end
