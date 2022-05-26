import numpy as np
import numpy as np
from math import cos, sin, atan, pi
from utils import *


def getMatrizA():
    A = []

    with open("task02/src/mat_A.dat") as file:
        for line in file.readlines():
            data = line.split()
            A.append(data)
    return A


def getMatrizB():
    B = []

    with open("task02/src/vet_B.dat") as file:
        for line in file:
            B.append(line.rstrip('\n'));
    return B


def main():

    ICOD = input('Entre com o ICOD do método desejado: ')
    DET = input('Digite 1 para calcular determinante e 0 para não calcular ')

    Amatrix = getMatrizA()
    Bmatrix = getMatrizB()

    A = []
    B = []

    for i in range(len(Amatrix)):
        A.append(list(np.float_(Amatrix[i])))

    B = list(np.float_(Bmatrix))

    if (ICOD == '1'):
        metodo_potencia(A, B, 1, 10**-3, DET)
    
    elif(ICOD == '2'):
        jacobi(A, 10**-3, DET)

    else:
        print("Insira um código válido")


def metodo_potencia(A, x, lamb, tolerancia, DET):
    n = len(A)
    res = 1
    iteracao = 0
    lamb = x[0]  # primeiro elemento do vetor deve ser 1
    while res > tolerancia:
        y = multiplicacao_matriz_vetor(A, x)
        lamb_aux = y[0]

        for i in range(n):
            y[i] = y[i] / lamb_aux

        res = modulo(lamb - lamb_aux) / modulo(lamb)
        lamb = lamb_aux
        x = y

        iteracao += 1
    print("lambda , x , nº iteração")
    if (DET == '1'):
        print(np.linalg.det(A))
    print(round(lamb, 2), arredondar(x,2), iteracao)
    return lamb, x, iteracao


def jacobi(A, tolerancia, DET):
    Aux = np.copy(A)
    n = np.shape(A)[0]
    X = np.eye(n)
    P = np.eye(n)
    iteracao = 0

    At, boo = simetrica(A)
    if boo == False:
        return (A, X, iteracao, boo)

    ok, linha, coluna = teste_convergencia(Aux, n, tolerancia)

    while (ok == False):
        ok, linha, coluna = teste_convergencia(Aux, n, tolerancia)

        phi = angulo_phi(Aux, linha, coluna)
        P = matrizRotacao(n, phi, linha, coluna)
        P_t = matrizTransposta(P)

        Aux = np.dot(P_t, Aux)
        Aux = np.dot(Aux, P)

        X = np.dot(X, P)
        iteracao += 1

        print('\n Iteracao nº ', iteracao, '\nangulo phi: ', phi,
              '\nmatriz autovalores A=\n', arredondar(Aux, 2), '\nmatriz Autovetores X=\n', arredondar(X, 2))
        if (DET == '1'):
            print(np.linalg.det(A))
    return (Aux, X, iteracao, ok)

    
if __name__ == '__main__':
    main()
