% @param  ������ �����
% @param  ������ �����������
% @return ������ ��������� ����� 
function [x_out] = shiftGridAccurate(x, d)
  N = length(x);
  % �������� �������� �������� �����
  v = zeros(1,N);
  for i = 1:length(d)-1
      % ���� � ������ ����������� ���� ������ ������ �����, ���� ���������
      if d(i)*d(i+1) > 0
          x_out = [];
          return;
      end
      % ��������� �������� �������� �����
      v_pm = (d(i+1) + d(i))/(d(i+1) - d(i));
      % ����������� ��������� ��������, ������� �������
      MUL_CONST = 1/(2*sqrt(3));
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
  x_out = [];
  % tau_max = 1;
  for i = 1:N-1
      x_out = [x_out, x(i) + v(i)*tau_shift];
  end
  x_out(N) = x(end);
end
