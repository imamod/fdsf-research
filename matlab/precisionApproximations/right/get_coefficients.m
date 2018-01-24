function [a, b] = get_coefficients(A, B, y0, z, N, k, baseSize)
for i = 1:size(A,2)
    for j = 1:size(A,2)
        if (j > 0 && j < N + 1)
            A(i,j) = (y0(i))^(2*(j-1));
        elseif (j > N && j < baseSize)
            A(i,j) = -z(i)*y0(i)^(2*(j-N-1));
        elseif (j == baseSize)
            A(i,j) = -z(i)*y0(i)^(2*N-2) + (y0(i)^(2*N));
        end
    end
end
% disp(A); 
% disp(B');
E = A\B';
% disp(E);

for j = 1:length(E)
    if (j > 0 && j <= N)
        a(j) = E(j);
    elseif (j > N && j <= baseSize)
        b(j-N) = E(j);
    end
end
a(N+1) = b(N) - (k+1)*(pi^2)/(3*k);