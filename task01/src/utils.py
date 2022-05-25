# utils.py>
import numpy as np



#Funções auxiliares

def matrizQuadrada(A): # verifica se uma matriz é quadrada
    n = len(A)  # numero de linhas 
    m = len(A[0])  # numero de colunas
    if n != m:
        return False
    return True


def matrizPositiva(A): # verifica se os números de uma matriz quadrada são positivos
  n = len(A) # tamanho da matriz quadrada
  for i in np.arange(n): # percorre linhas da matriz quadrada
    for j in np.arange(n): # percorre colunas da matriz quadrada
        if A[i][j] <= 0:
          return False
  return True


def determinanteLU(A): # multiplicar diagonal principal
    n = np.shape(A)[0]  # numero de linhas
    det = 1
    for i in np.arange(n): # percorre o tamanho da matriz
        det *= A[i][i] # multiplica diagonal principal
    return det



def identidade(N):
    'Cria uma matriz identidade'
    identidade = []
    for i in range(0,N):
        linha = [0] * N
        linha[i] = 1
        identidade.append(linha)
    return identidade


def modulo(x): # calcula o módulo de um número 
  if (x < 0):
      return (-x)
  return x


def matrizTransposta(A):
    n = len(A) # tamanho da matriz A
    L = np.zeros((n, n)) # Cria uma matriz vazia do tamanho de A

    for i in range(n): # percorre linhas da matriz 
        for j in range(n): # percorre colunas da matriz 
            L[j][i] = A[i][j] # preenche a matriz transposta
    return L


def condicao_convergencia(A): # para convergir, a matriz A deve ser diagonal dominante
    converge = True
    n = len(A) # tamanho da matriz
    for i in range(n): # percorre as linhas da matriz
            somai, somaj = 0, 0
            for j in range(n): # percorre as colunas da matriz
                if (i != j): #verifica se não é a diagonal principal
                  somai += modulo(A[i][j]) # soma linha
                  somaj += modulo(A[j][i]) # soma coluna
            a = modulo(A[i][i]) # pega o valor da diagonal principal
            if ((a < somai) or (a < somaj)): # o valor da diagonal principal precisa ser maior ou igual que a soma das linhas e maior ou igual que a soma das colunas 
              converge = False
    return converge