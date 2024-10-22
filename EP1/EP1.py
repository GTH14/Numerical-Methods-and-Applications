import numpy as np
from time import time
def QR(A0):
    t1 = time()
    Ak1 = A0
    n = len(A0)
    I = np.identity(n)
    Vk1 = I
    k = 0   
    erro = 1
    erro_max = 10**(-6)
    # Autovalores = []
    # muk = 0

    k=0
    while erro>erro_max:
        Ak = Ak1
        Vk = Vk1
        Qk = I
        Rk = Ak 
        for j in range(n-1):
            Rj,Qj = Givens(j,Rk)
            Qj = np.transpose(Qj)
            Qk = np.dot(Qk,Qj)
            Rk = Rj
            # print(Qj)
        Ak1 = np.dot(Rk,Qk)
        Vk1 = np.dot(Vk,Qk)
        erro = np.abs(Ak1[1][0])
        for i in range(len(Ak1)-2):
            if erro < np.abs(Ak1[i+2][i+1]):
                erro = np.abs(Ak1[i+2][i+1])
        k+=1
    Autovalores = Ak1
    Autovetores = Vk1
    t2 = time()

    return Autovetores, Autovalores,k,t2-t1

def cut(A):
    n = len(A)
    B = np.identity(n-1)
    for i in range(n-1):
        for j in range(n-1):
            B[i][j] = A[i][j]
    return B    
    
def Givens(i,A):
    alpha = A[i][i]
    beta = A[i+1][i]
    n = len(A)
    if abs(alpha)>abs(beta):
        tau = -beta/alpha
        c = 1/np.sqrt(1+tau**2)
        s = c*tau
    else:
        tau = -alpha/beta
        s = 1/np.sqrt(1+tau**2)
        c = s*tau
    Q = np.identity(n)
    Q[i][i] = c
    Q[i][i+1] = -s
    Q[i+1][i] = s
    Q[i+1][i+1] = c
    R = np.dot(Q,A)
    return R,Q

def Calc_muk(Ak,m):
    alphak_n1 = Ak[m-1][m-1]
    alphak_n = Ak[m][m]
    betak_n1 = Ak[m][m-1]
    dk = (alphak_n1-alphak_n)/2
    muk = alphak_n + dk-np.sign(dk)*(dk**2 + (betak_n1)**2)**(1/2)
    return muk

def matrizA(n):
    A = np.identity(n)
    for i in range(n):
        A[i][i] = 2
        if i < n-1:
            A[i][i+1] = -1
            A[i+1][i] = -1
    return A
A = matrizA(4)
# print(A)

V, Lambda, k, t = QR(A)
print(k)
Lambda_analitico = []
V_analitico = []
n = len(A)
for i in range(n):
    v = []
    Lambda_analitico.append(2*(1-np.cos((i+1)*np.pi/(n+1))))
    for j in range(n):
        v.append(np.sin((j+1)*(i+1)*np.pi/(n+1)))
    
    V_analitico.append(v)

Lambda_analitico = np.array(Lambda_analitico)
V_analitico = np.transpose(np.array(V_analitico))

print(V_analitico)
print(np.dot(A,np.transpose(np.transpose(V_analitico)[0]))/Lambda_analitico[0])
# # print(V_analitico)

# # print(Lambda)

# print(t)
# print(Lambda_analitico)
