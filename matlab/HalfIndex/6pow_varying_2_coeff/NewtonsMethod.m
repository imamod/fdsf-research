function res = NewtonsMethod(x, k)
    a0 = 3 * (k + 7.0 / 8); % ��������� �����������
    res = NewtonInternal( x, a0, k, 0.001);
end

function f = func( x, a, k)
      f = a - 3 * (k + 7.0 / 8) * (1 + exp(x - a / 3));
end

%  ������ �����������
function df = dfunc(x, a, k)
   df = 1 + (k + 7.0 / 8) * exp(x - a / 3);
end

%  ����� ������� 
function res = NewtonInternal(x, a0, k, eps)

    a1 = a0 - func(x, a0, k) ./ dfunc(x, a0, k);
    while (abs(a0 - a1) > eps) 
        a0 = a1;
        a1 = a0 - func(x, a0, k) ./ dfunc(x, a0, k);
    end

    res = a1;
end
