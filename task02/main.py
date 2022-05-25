import numpy as np


def getMatrizA():
    A = []

    with open("mat_A.dat") as file:
        for line in file.readlines():
            data = line.split()
            A.append(data)
    return A


def getMatrizB():
    B = []

    with open("vet_B.dat") as file:
        for line in file:
            B.append(line.rstrip('\n'));
    return B


def main(ICOD):
    Amatrix = getMatrizA()
    Bmatrix = getMatrizB()

    A = []
    B = []

    for i in range(len(Amatrix)):
        A.append(list(np.float_(Amatrix[i])))

    B = list(np.float_(Bmatrix))

    A = [[1, 2, 2], [4, 4, 2], [4, 6, 4]]
    B = [1, 6, 10]

    if (ICOD == 1):
        metodo_potencia(A, B)


    else:
        print("Insira um código válido")


# Funções Auxiliares

def multiplicacao_matriz_vetor(A, x):
  n = len(A)
  y = np.zeros(n)

  for i in range(n):
      soma = 0
      for j in range(n):
          soma += A[i][j] * x[j]
      y[i] = soma

  return y


def modulo(x):
    if x < 0:
      return -x
    return x


# Task 2

def metodo_potencia(A, x, lamb=1, tolerancia=10**-3):
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

    print(lamb, x, iteracao)
    return lamb, x, iteracao


if __name__ == '__main__':
    main(1)
