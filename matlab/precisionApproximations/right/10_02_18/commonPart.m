function max_max_delta = commonPart(y0, Y, I_base, I_add, k, N)
    baseSize = 2*N;
    C1 = (k+1)*k*(pi^2)/6;
    % disp('---------------------------');
    % disp('y size: '); disp(y0.length());
    % disp('I_base: '); disp(I_base.length());
    % disp('---------------------------');
    z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);
    Z = (((I_add*(k+1)./(Y.^(k+1))).^(2/k) - 1).*(Y.^2))*k/(2*C1);

    % Задаем матрицы A и B
    B = (z(1,:) - ones(1,baseSize));
%     A = setA(baseSize, y0, z, N);
    M = 4;
    A = setExtendedA(baseSize, y0, z, N, M);
%     disp(A); 
%     disp(B');
    E = A\B';
%     disp(E);

    a = E(1:N);
    b = E(N+1:end);

    disp('lg(cond(A)):'); disp(log10(cond(A)));
    disp('--------------------------------');
    disp('Коэффициенты а:'); disp(a);
    disp('----------------------------');
    disp('Коэффициенты b:'); disp(b);
    disp('----------------------------');


%     F_z_base = approximateFz(y0, a, b, N);
    F_z_base = approximateExtendedFz(y0, a, b, N, M);
    delta_z = F_z_base./z-1;
    F_base = getFbyFz(y0, F_z_base, C1, k);
    delta_base = F_base./I_base-1;

    %---------------------------------------
    % Добавим вспомогательную сетку

%     F_z = approximateFz(Y, a, b, N);
    F_z = approximateExtendedFz(Y, a, b, N, M);
    delta_z_additional = F_z./Z-1;
    F = getFbyFz(Y, F_z, C1, k);
    delta_additional = F./I_add-1;

    plotDelta(Y, delta_z_additional, y0, delta_z, 'Относительная погрешность z');
    plotDelta(Y, delta_additional, y0, delta_base, 'Относительная погрешность I');
    % disp((I-F)');

    max_delta = intervalMaximums(baseSize, delta_additional);

    prefix = strcat('extremN', num2str(N));
    filename = strcat(prefix, '.txt');
%     writeFile(filename, '%s\n', num2str(log10(max_delta)));
    f = fopen(filename , 'w');
    fprintf(f, '%e\n', max_delta);
    fclose(f);

%     tau = 0.75;
%     disp(y0);
%     disp('-----eta--------');
%     eta = shiftGrid(y0, max_delta, tau);
%     disp(1./eta);
    
    max_max_delta = max(max_delta);
end

function A = setExtendedA(baseSize, y0, z, N, M)
    for i = 1:baseSize
        for j = 1:baseSize
            if (j > 0 && j < M + 1)
                A(i,j) = y0(i)^(-2*j);
            elseif (j > M && j < N + 1)
                A(i,j) = y0(i)^(-M-j);
            elseif (j > N && j < N + M + 1)
                A(i,j) = -z(i)*y0(i)^(2*N - 2*j);
            else
                A(i,j) = -z(i)*y0(i)^(N - M - j);
            end
        end
    end
end

function F_z = approximateExtendedFz(y, a, b, N, M)
    F_z = [];
    for j = 1:length(y)
        S1 = 0; S2 = 0; 
        for n = 1:N
            if n < M+1
                S1 = S1 + a(n).*y(j).^(-2*n);
            else
                S1 = S1 + a(n).*y(j).^(-M-n);
            end
        end
        for m = 1:N
            if m < M + 1
                S2 = S2 + b(m).*y(j).^(-2*m);
            else
                S2 = S2 + b(m).*y(j).^(-M-m);
            end
        end
        F_z = [F_z, (1 + S1)/(1 + S2)];
    end
end

function A = setA(baseSize, y0, z, N)
    for i = 1:baseSize
        for j = 1:baseSize
            if (j > 0 && j < N + 1)
                A(i,j) = (y0(i))^(-2*j);
            elseif (j > N && j < baseSize+1)
                A(i,j) = -z(i)*y0(i)^(2*(N-j));
            end
        end
    end
end

function F_z = approximateFz(y, a, b, N)
    F_z = [];
    for j = 1:length(y)
        S1 = 0; S2 = 0; 
        for n = 1:N
            S1 = S1 + a(n).*y(j).^(-2*n);
        end
        for m = 1:N
            S2 = S2 + b(m).*y(j).^(-2*m);
        end
        F_z = [F_z, (1 + S1)/(1 + S2)];
    %     delta_z(j) = F_z_base(j)/z(j)-1;
    %%%     z = ((I_base*(k+1)./(y0.^(k+1))).^(2/k) - 1).*(y0.^2)*k/(2*C1);
        
    %     delta_base(j) = F_base(j)/I_base(j)-1;
    end
end

function F = getFbyFz(y, F_z, C1, k)
    F = (((F_z*2*C1./(k*y.^2)) + 1).^(k/2)).*(y.^(k+1))/(k+1);
end

function plotDelta(Y, delta, y, delta_base, desc)
    figure
    grid on, hold on
    xlabel('y'); %ylabel('d*10^1^0');
    title(desc);
    plot(Y, delta, 'k','linewidth', 2.5);
    plot(y, delta_base, 'k*','linewidth',5)
end

function max_delta = intervalMaximums(baseSize, delta_additional)
    disp('lg(dc):');disp(log10(max(abs(delta_additional))));
    disp('-------------------');
    disp('Экстремумы (log10(max_delta)):');
    span = 11;
    max_delta = max(abs(delta_additional(1:span)));
    deltas = [sprintf('i = 1: %f\n', log10(max_delta))];
    for i=1:baseSize-1
%         if (N > 9 && (i == (baseSize-1)))
%             break;
%         end
        if (i == (baseSize-1))
            currentDelta = max(abs(delta_additional((i*span):end)));
            max_delta = [max_delta, currentDelta];
            break;
        end
        currentDelta = max(abs(delta_additional((i*span):(i+1)*span)));
        max_delta = [max_delta, currentDelta];
        deltas = [deltas, sprintf('i = %d: %f\n', i+1, log10(currentDelta))];
    end
    disp(deltas);
end