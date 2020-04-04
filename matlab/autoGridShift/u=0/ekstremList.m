% @param - профиль погрешности
% return - 
function [ekstrDelta] = ekstremList(U, N) 
    % Рачет экстремумов    
    disp('Экстремумы U:');
    span = 20;
    mid = span/2;
    if (U(mid) < 0)
        ekstrDelta = min(U(1:span));
    else
        ekstrDelta = max(U(1:span));
    end
    for i=1:N-1
        if (U(i*span + mid) < 0)
            lastEkstr = min((U(i*span:(i+1)*span)));
        else
            lastEkstr = max(U(i*span:(i+1)*span));
        end
        ekstrDelta = [ekstrDelta, lastEkstr];
    end
    disp(ekstrDelta');
    disp('-------------------');
end