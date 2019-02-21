function [H, D1] = D1_scratch1(n)
    H = eye(n);
    H(1, 1) = 0.5;
    H(end, end) = 0.5;
    Q = zeros(n);
    for i=1:n-2
        Q(i + 1, i) = -0.5;
        Q(i + 1, i + 2) = 0.5;
    end
    Q(1, 1) = -0.5;
    Q(1, 2) = 0.5;
    Q(end, end-1) = -0.5;
    Q(end, end) = 0.5;
    D1 = inv(H)*Q;
end