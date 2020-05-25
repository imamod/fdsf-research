% param ����� ������� �����
% param ����� �������������� �����
% param ������ �������� �������
% return ������������������ �������� ������� � �����
function [U] = rationalApprox(x, X, U_prec)
    % ������� ������� ������ ���� ����� ����� �����
    baseGridSize = length(x);
    % N - ����� ������ � ���������
    % M - ����� ������ � �����������
    % ���� ����� ����� ������, �� ������� ���������� � ��������� � �����������,
    % ����� � ����������� �� 1 ������, ��� � ���������
    if (mod(baseGridSize - 1, 2) == 0)
        N = (baseGridSize + 1)/2;
        M = N;
        disp('N = ');disp(N);
        disp('M = ');disp(M);
    else
        N = baseGridSize/2;
        M = N + 1;
        disp('N = ');disp(N);
        disp('M = ');disp(M);
    end
    approxSize = M + N - 1;
    
    % ������ ������� A � B
    B = transpose(U_prec);
    A = zeros(approxSize,approxSize);
    for i = 1:baseGridSize
        for j = 1:approxSize
            if (j < N + 1)
                A(i,j) = x(i)^(j-1);
            elseif (j <= approxSize)
                A(i,j) = -U_prec(i)*x(i)^(j-N);
            end
        end
    end
   % disp(A); 
   % disp(B);
    E = A\B;
  %  disp(E);

    a = ones(1,N); b = ones(1,M);
    for j = 1:length(E)
        if (j < N + 1)
            a(j) = E(j);
        elseif (j <= approxSize)
            b(j-N+1) = E(j);
        end
    end
    
    disp('--------------------------------');
    disp('������������ �:'); disp(a');
    disp('----------------------------');
    disp('������������ b:'); disp(b');
    disp('----------------------------');

    gridSize = length(X);
    U = zeros(1,gridSize);
    for j = 1:gridSize
        S1 = 0; S2 = 0; 
        for n = 1:N
            S1 = S1 + a(n)*X(j)^(n-1);
        end
        for m = 1:M
            S2 = S2 + b(m)*X(j)^(m-1);
        end
        U(j) = S1/S2;
    end

end
