% @param  вектор узлов
% @param  вектор экстремумов
% @return вектор сдвинутой сетки 
function [x_out] = shiftGridAccurate(x, d)
  N = length(x);
  % Скорости движения соседних узлов
  v = zeros(1,N);
  for i = 1:length(d)-1
      % если в списке экстремумов есть соседи одного знака, срыв алгоритма
      if d(i)*d(i+1) > 0
          x_out = [];
          return;
      end
      % добавляем скорости движения узлов
      v_pm = (d(i+1) + d(i))/(d(i+1) - d(i));
      % настроечная константа скорости, рабочий вариант
      MUL_CONST = 1/(2*sqrt(3));
      v(i+1) = (x(i+2) - x(i))*v_pm*MUL_CONST;
  end
  %disp('Векторы скоростей'); disp(v);
  
  % Условно большое tau
  TAU_BIG = 10^10;
  % расчет tau
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
  
  % выборка минимального tau
  tau_min = min(tau);
  % настроечный коэффициент tau, рабочий вариант
  A = 0.2; % Больше 0.2 делать смысла не имеет, меньше 0.1 тоже не имеет, потому что овчинка выделки не стоит
  % домножаем на настроечный коэффициент  
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
