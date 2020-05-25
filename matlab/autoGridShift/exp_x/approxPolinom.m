% param ����� ������� �����
% param ����� �������������� �����
% param ������ �������� �������
% return ������������������ �������� ������� � �����
function [U] = approxPolinom(x, X, U_prec)
    
    N = length(x);

    % ���������������� ������� exp(x)
    % U_prec = zeros(1, N + 1);
    B = transpose(U_prec);
    A = zeros(N, N);
    for j = 1:N
        A(:,j) = x.^(j-1);
    end

    %disp('---B---');
    %disp(B);

    %disp('---A---');
    %disp(A);
    a = A\B;
    %disp('---������������---');
    %disp(a);

    U = [];
    for i=1:length(X)
        u = 0;
        for j = 1:N
            u = u + a(j)*X(i)^(j-1);
        end
        U = [U, u];
    end
end
