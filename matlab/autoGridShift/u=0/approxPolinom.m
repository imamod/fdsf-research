% param ����� ������� �����
% param ����� �������������� �����
% return ������������������ �������� ������� � �����
function [U] = approxPolinom(x, X)
    
    N = length(x)-1;

    % ���������������� ������� = 0
    U_prec = zeros(1, N + 1);

    B = transpose(U_prec - x.^(N + 1));
    A = zeros(N+1, N+1);
    for j = 0:N
        A(:,j+1) = x.^j;
    end

    %disp('---B---');
    %disp(B);

    %disp('---A---');
    %disp(A);
    a = A\B;
    %disp('---������������---');
    %disp(a);

    % ����������� �� ���������
    a(N+2)=1;

    U = [];
    for i=1:length(X)
        u = X(i)^(N+1);
        for j = 1:N+1
            u = u + a(j)*X(i)^(j-1);
        end
        U = [U, u];
    end
end
