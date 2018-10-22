import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lin

# Чтение данных
def readData(filename):
    value = []
    with open(filename,'r') as f:
        value = [float(line.strip()) for line in f]

    return value

# Построение графика погрешности
def plot(x1, y1, x2, y2, title):
    plt.plot(x1, y1, 'ko', ms=6, lw=2, mfc='black')
    plt.plot(x2, y2, 'k', lw=2.5)

    plt.grid()
    plt.xlabel('y')
    plt.ylabel('d') #ylabel('d*10^1^0')
    plt.title(title)

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)
    plt.show()

def commonPart(y0, Y, I_base, I_add, k, N):
    baseSize = 2*N # Добавляем 1, чтобы получился размер 2N, особенность питона
    C1 = (k+1)*k*(math.pi**2)/6
    # print('---------------------------')
    # print('y size: ' + y0.length())
    # print('I_base: ' + I_base.length())
    # print('---------------------------')
    z = []
    for i in range(len(y0)):
        z.append(((I_base[i]*(k+1)/(y0[i]**(k+1)))**(2/k) - 1)*(y0[i]**2)*k/(2*C1))
    #print(z[:, np.newaxis])

    Z = []
    for i in range(len(Y)):
        Z.append(((I_add[i]*(k+1)/(Y[i]**(k+1)))**(2/k) - 1)*(Y[i]**2)*k/(2*C1))
    #print(Z)

    # Задаем матрицы A и B
    B = np.array(z, dtype=float) - np.ones(baseSize, dtype=float)
    A = np.zeros((baseSize,baseSize))
    for i in range(baseSize):
        for j in range(baseSize):
            if j < N:
                A[i,j] = (y0[i])**(-2*(j+1))
            elif j > N:
                A[i,j] = -z[i]*(y0[i]**(2*(N-j)))

    #print(A)
    B = B.transpose()
    #print(B[:,np.newaxis])
    E = lin.lstsq(A, B)
    print(E);

    a = np.zeros(N)
    b = np.zeros(N)
    for j in range(len(E)):
        if j < N:
            a[j] = E[j]
        else:
            b[j-N] = E[j]

    #print('lg(cond(A)):'); print(log10(cond(A)));
    #print('--------------------------------');
    #print('Коэффициенты а:'); print(a');
    #print('----------------------------');
    #print('Коэффициенты b:'); print(b');
    #print('----------------------------');

    F_base = zeros(baseSize)
    delta_base = zeros(baseSize)
    for j in range(baseSize):
        S1 = 0; S2 = 0; 
        for n in range(1,N):
            S1 = S1 + a[n]*y0[j]^(-2*n)
        for m in range(1,N):
            S2 = S2 + b[m]*y0[j]^(-2*m)
        F_z_base[j] = (1 + S1)/(1 + S2)
        delta_z[j] = F_z_base[j]/z[j]-1
    #     z = ((I_base*(k+1)/(y0**(k+1)))**(2/k) - 1)*(y0**2)*k/(2*C1);
        F_base[j] = (((F_z_base[j]*2*C1/(k*y0[j]**2)) + 1)**(k/2))*(y0[j]**(k+1))/(k+1)
        delta_base[j] = F_base[j]/I_base[j]-1
    #---------------------------------------
    # Добавим вспомогательную сетку

    F = np.zeros(len(Y));
    delta_additional = np.zeros(len(Y));
    for j in range(1,len(Y)):
        S1 = 0
        S2 = 0
        for n in range(1,N):
            S1 = S1 + a[n]*Y[j]^(-2*n)
        for m in range(1,N):
            S2 = S2 + b[m]*Y[j]^(-2*m)

        F_z[j] = (1 + S1)/(1 + S2)
        delta_z_additional[j] = F_z[j]/Z[j]-1
        F[j] = (((F_z[j]*2*C1/(k*Y[j]**2)) + 1)**(k/2))*(Y[j]**(k+1))/(k+1)
        delta_additional[j] = F[j]/I_add[j]-1

    plot(y0, delta_z, Y, delta_z_additional, 'Относительная погрешность z')
    # print((I-F)')

    #print('lg(dc):');print(log10(max(abs(delta_additional))));
    #print('-------------------');
    #print('Экстремумы:');
    #span = 11;
    #max_delta = max(abs(delta_additional(1:span)));
    #print('log10(max_delta): ')
    #print(log10(max_delta));
    #for i=1:baseSize-1:
        #if (N > 9 && (i == (baseSize-1)))
            #break;
        #end
        #if (i == (baseSize-1))
            #currentDelta = max(abs(delta_additional((i*span):end)));
            #max_delta = [max_delta, currentDelta];
            #break;
        #end
        #currentDelta = max(abs(delta_additional((i*span):(i+1)*span)));
        #max_delta = [max_delta, currentDelta];
        ##disp('log10(max_delta): '); disp(log10(currentDelta));
    #end

    #prefix = 'extremN' + str(N);
    #filename = prefix + '.txt';
    #f = open(filename , 'w');
    #f.write(f, '%e\n', max_delta);
    #f.close(); 

    plot(y0, delta_base, Y, delta_additional, 'Относительная погрешность I')

   # max_max_delta = max(max_delta);


if __name__ == '__main__':
    k = 1/2;

    # N = [3, 5, 7, 9] #, 10, 11, 13]
    N = 9;
    #for item in N:
        # prefix = 'k12/n' + str(N)
        # prefix = 'linearGrid/k12/n' + str(N)
        # prefix = 'linTrigAlpha07/k12/n' + str(N)
        # prefix = 'changedWeights/k12/n' +  str(N)
    prefix = 'testShift/k12/n' + str(N)
    #prefix = 'testShift/tau1/k12/n' + str(N)
        # prefix = 'x5linTrig/k12/n' + str(N)

    filename = prefix + '/y0.txt'
    y0 = readData(filename)
    filename = prefix + '/Y.txt'
    Y = readData(filename)
    filename = prefix + '/I_base.txt'
    I_base = readData(filename)
    filename = prefix + '/I_add.txt'
    I_add = readData(filename)

    commonPart(y0, Y, I_base, I_add, k, N)
