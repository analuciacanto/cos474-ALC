import numpy as np
from math import cos, sin, atan, pi


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


def matrizTransposta(A):
    n = len(A) # tamanho da matriz A
    L = np.zeros((n, n)) # Cria uma matriz vazia do tamanho de A

    for i in range(len(A)): # percorre linhas da matriz 
        for j in range(len(A)): # percorre colunas da matriz 
            L[j][i] = A[i][j] # preenche a matriz transposta
    return L


def simetrica(A):
    n = np.shape(A)[0]
    m = np.shape(A)[1]
    if (n != m):
        return (A, False)
    At = matrizTransposta(A)

    for i in np.arange(n):
        for j in np.arange(n):
            if (At[i][j] != A[i][j]):
                return (At, False)

    return (At, True)


def matrizRotacao(n, phi, i, j):
    I = np.eye(n)
    I[i, i] = cos(phi)
    I[j, j] = cos(phi)
    I[i, j] = -sin(phi)
    I[j, i] = sin(phi)
    return I

def arredondar(A, x):
    n = np.shape(A)
    aux = len(n)

    n = np.shape(A)[0]

    if aux == 1:
        for i in np.arange(n):
            A[i] = round(A[i], x)
    else:
        m = np.shape(A)[1]
        for i in np.arange(n):
            for j in np.arange(m):
                A[i, j] = round(A[i, j], x)

    return A

def angulo_phi(A, i, j):
    if (A[i, i] == A[j, j]):
        phi = (pi / 4)
        return phi

    phi = atan((2*A[j, i])/(A[i, i] - A[j, j]))
    phi /= 2
    return phi



def teste_convergencia(A, n, tolerancia):
  
    linha, coluna = 0, 1

    boolean = True
    for i in np.arange(n):
        for j in np.arange(i+1, n):
            a = modulo(A[i, j])
            if (a > tolerancia):
                boolean = False
                if (a > modulo(A[linha, coluna])):
                    linha = i
                    coluna = j

    return (boolean, linha, coluna)

