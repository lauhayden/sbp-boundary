% 1st derivative, p=2 SBP operator
% From SBP SAT Review paper

function [H,D1]=D1_4(n)
h11 = 17.d0/48.d0; h22 = 59.d0/48.d0; h33 = 43.d0/48.d0; h44 = 49.d0/48.d0;
q11 = -1.d0/2.d0; q12 = 59.d0/96.d0; q13 = -1.d0/12.d0; q14 = -1.d0/32.d0;
q23 = 59.d0/96.d0; q24 = 0.d0;q34 = 59.d0/96.d0;
%setting up H and D1
H = eye(n,n);
H(1,1) = h11;H(n,n) = h11;H(2,2) = h22;H(n-1,n-1) = h22;
H(3,3) = h33;H(n-2,n-2) = h33;H(4,4) = h44;H(n-3,n-3) = h44;
Q = zeros(n,n);
Q(1,1) = q11;Q(1,2) = q12;Q(1,3) = q13; Q(1,4) = q14;
Q(2,1) = -q12; Q(2,3) = q23; Q(2,4) = q24;Q(3,1) = -q13;Q(3,2) = -q23;
Q(3,4) = q34; Q(3,5) = -1.d0/12.d0;Q(4,1) = -q14;Q(4,2) = -q24;
Q(4,3) = -q34; Q(4,5) = 2.d0/3.d0; Q(4,6) = -1.d0/12.d0;
%internal nodes
for i =5:n-4
    Q(i,i-2:i+2) = ...
    [1.d0/12.d0,-2.d0/3.d0,0.d0,2.d0/3.d0,-1.d0/12.d0];
end
%bottom portion of the matrix
for i = 1:6
    for j = 1:6
        Q(n-(i-1),n-(j-1)) = -Q(i,j);
    end
end
D1 = inv(H)*Q;