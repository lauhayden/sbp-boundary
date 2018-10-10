% 1st derivative, p=3 SBP operator
% From SBP SAT Review paper

function [H,D1]=D1_6(n,q56)
% optimized value comment out to use other values
q56 =5591070156686698065364559.d0/7931626489314500743872000.d0;
h11 = 0.13649d5 / 0.43200d5; h22 = 0.12013d5 / 0.8640d4;
h33 = 0.2711d4 / 0.4320d4; h44 = 0.5359d4 / 0.4320d4;
h55 = 0.7877d4 / 0.8640d4; h66 = 0.43801d5 / 0.43200d5;
q11 = -1.d0/2.d0;q12 = -0.953d3 / 0.16200d5 + q56;
q13 = 0.715489d6 / 0.259200d6 - (4.d0 * q56);
q14 = -0.62639d5 / 0.14400d5 + (6.d0 * q56);
q15 = 0.147127d6 / 0.51840d5 - (4.d0 * q56);
q16 = -0.89387d5 / 0.129600d6 + q56;
q23 = -0.57139d5 / 0.8640d4 + (10.d0 * q56);
q24 = 0.745733d6 / 0.51840d5 - (20.d0 * q56);
q25 = -0.18343d5 / 0.1728d4 + (15.d0 * q56);
q26 = 0.240569d6 / 0.86400d5 - (4.d0 * q56);
q34 = -0.176839d6 / 0.12960d5 + (20.d0 * q56);
q35 = 0.242111d6 / 0.17280d5 - (20.d0 * q56);
q36 = -0.182261d6 / 0.43200d5 + (6.d0 * q56);
q45 = -0.165041d6 / 0.25920d5 + (10.d0 * q56);
q46 = 0.710473d6 / 0.259200d6 - (4.d0 * q56);
q47 = 1.d0/6.d1;q57 = -3.D0/2.d1; q58 = 1.d0/6.d1;
%setting up H and D1
H = eye(n,n);
H(1,1) = h11;H(n,n) = h11;H(2,2) = h22;H(n-1,n-1) = h22;
H(3,3) = h33;H(n-2,n-2) = h33;H(4,4) = h44;H(n-3,n-3) = h44;
H(5,5) = h55;H(n-4,n-4) = h55;H(6,6) = h66;H(n-5,n-5) = h66;
Q = zeros(n,n);
Q(1,1) = q11;Q(1,2) = q12;Q(1,3) = q13; Q(1,4) = q14;
Q(1,5) = q15; Q(1,6) = q16;Q(2,1) = -q12; Q(2,3) = q23;
Q(2,4) = q24; Q(2,5) = q25;Q(2,6) = q26;Q(3,1) = -q13;Q(3,2) = -q23;
Q(3,4) = q34; Q(3,5) = q35;Q(3,6) = q36;
Q(4,1) = -q14;Q(4,2) = -q24;Q(4,3) = -q34; Q(4,5) = q45;
Q(4,6) = q46; Q(4,7) = q47;Q(5,1) = -q15; Q(5,2) = -q25;
Q(5,3) = -q35; Q(5,4) = -q45;Q(5,6) = q56; Q(5,7) = q57;
Q(5,8) = q58;Q(6,1) = -q16; Q(6,2) = -q26;Q(6,3) = -q36; Q(6,4) = -q46;
Q(6,5) = -q56; Q(6,7) = 3.d0/4.d0;Q(6,8) = -3.d0/2.d1; Q(6,9) = 1.d0/6.d1;
%internal nodes
for i =7:n-6
    Q(i,i-3:i+3) = ...
    [-1.d0/6.d1, 3.d0/2.d1, -3.d0/4.d0, 0.d0, 3.d0/4.d0, -3.d0/2.d1, 1.d0/6.d1];
end
%bottom portion of the matrix
for i = 1:9
    for j = 1:9
        Q(n-(i-1),n-(j-1)) = -Q(i,j);
    end
end
D1 = inv(H)*Q;