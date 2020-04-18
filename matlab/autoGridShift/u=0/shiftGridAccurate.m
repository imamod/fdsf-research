% @param  ������ �����
% @param  ������ �����������
% @return ������ ��������� ����� 
function [X] = shiftGridAccurate(x, d)
  N = length(x);
  % �������� �������� �������� �����
  v = zeros(1,N);
  for i = 1:length(d)-1
      % ���� � ������ ����������� ���� ������ ������ �����, ���� ���������
      if d(i)*d(i+1) > 0
          return;
      end
      % ��������� �������� �������� �����
      v_pm = (d(i+1) + d(i))/(d(i+1) - d(i));
      % ����������� ��������� ��������, ������� �������
      MUL_CONST = 1/(2*sqrt(3));
      % �������� ��� ���������� ������������ TODO: ������� ����� ������
      % X_i1mi = x(i+1) - x(i);
      % X_i2mi1 = x(i+2) - x(i+1);
      % nom = X_i1mi^2 + X_i1mi*X_i2mi1 + X_i2mi1^2;
      % denom = (x(i+2) - x(i))*(2*X_i1mi^2 + 5*X_i1mi*X_i2mi1 + 2*(X_i2mi1^2));
      % MUL_CONST = (nom^(3/2))/denom;
      v(i+1) = (x(i+2) - x(i))*v_pm*MUL_CONST;
  end
  %disp('������� ���������'); disp(v);
  
  % ������� ������� tau
  TAU_BIG = 10^10;
  % ������ tau
  tau=[];
  for i = 1:length(v)-1
    denom = v(i) - v(i+1);
    if denom < 0
        nextTau = TAU_BIG;
    else
        nextTau = (x(i+1) - x(i))/denom;
    end
    tau = [tau, nextTau];
  end
  disp(length(tau));
  
  % ������� ������������ tau
  tau_min = min(tau);
  % ����������� ����������� tau, ������� �������
  A = 0.2; % ������ 0.2 ������ ������ �� �����, ������ 0.1 ���� �� �����, ������ ��� ������� ������� �� �����
  % ��������� �� ����������� �����������  
  tau_shift = tau_min*A;
  constTau = 1; 
  if tau_shift > constTau
    tau_shift = constTau;
  end
  disp('tau_shift: '); disp(tau_shift);
  X = [];
  % tau_max = 1;
  for i = 1:N-1
      %%% �������� ����������� ��������� ��������� (������� �������,
      %%% ������������� ����������) -  ����������
      %if i < N-1
      %    tau_shift = min([A*tau(i), A*tau(i+1), tau_max]);
      %else
      %    tau_shift = min([A*tau(i), tau_max]);
      %end
      %disp('tau_shift: '); disp(tau_shift);
      X = [X, x(i) + v(i)*tau_shift];
  end
  X(N) = x(end);  
end
