import numpy as np
from utils import *

def getMatrizA():
    A = []

    with open("task01/src/mat_A.dat") as file:
        for line in file.readlines():
            data = line.split()
            A.append(data)
    return A

def getMatrizB():
    B = []

    with open("task01/src/vet_B.dat") as file:
        for line in file:
            B.append(line.rstrip('\n'));
    return B


def main(): 

    ICOD = input('Entre com o ICOD do método desejado: ')
    DET = input("Digite 1 para calcular determinante e 0 para não: ")

    Amatrix = getMatrizA()
    Bmatrix = getMatrizB()

    A = []
    B = []

    for i in range(len(Amatrix)):
        A.append(list(np.float_(Amatrix[i])))
  
    B = list(np.float_(Bmatrix))    
 
    if (ICOD == '1'):
        decomposicaoLU(A, B, len(A), DET)        
    elif (ICOD == '2'):
        cholesky(A, B, len(A), DET)
    elif (ICOD == '3'):
        jacobi(A, B, 10**-3, DET)
    elif (ICOD == '4'):
        gauss_seidel(A, B, 10**-3, DET)      
    else:
        print("Insira um código válido")



def residuo_vetor(A, x, x_aux): # calcula o resíduo para o método iterativo de jacobi e gauss_seidel usando a norma euclidiana
  n = len(A) # tamanho da matriz
  num, den = 0.0, 0.0 # inicializa numerador e denominador com zero 
  for i in range(n): # percorre os vetores somando cada elemento elevado ao quadrado 
      num += (x_aux[i] - x[i])**2 # faz a diferença entre x_aux menos x ao quadrado
      den += (x_aux[i])**2 # eleva o valor de x_aux ao quadrado
  num = num**0.5 # tira a raiz da soma para o numerador 
  den = den**0.5 # tira a raiz da soma para o denominador 
  res = num/den # faz a divisão final do resíduo
  return res

def substituicoes_sucessivas(A, B, N):    

    'Resolve o sistema linear triangular inferior Ax = b'
    'A deve ser portanto uma matriz triangular inferior'

    y = N * [0]
    y[0] = B[0]/A[0][0]

    for i in range(1, N):
        sum = 0
        for j in range(0, i):
            sum = sum + A[i][j]*y[j]
        y[i] = (B[i]-sum)/A[i][i]
    return y

def substituicoes_retroativas(A, Y, N):
    'Resolve o sistema linear triangular superior Ax = b'
    'A deve ser portanto uma matriz triangular superior'

    #Inicializar vetor x com tamanho N
    x = N* [0]
    x[N-1] = Y[N-1] / A[N-1][N-1]

    for i in range(N-1, -1, -1):
        sum = 0
        for j in range(i+1, N):
            sum = sum + A[i][j] * x[j]
        x[i] = (Y[i] - sum)/A[i][i]

    return x

def decomposicaoLU(A, B, N, DET):
    quadrada = matrizQuadrada(A)

    if quadrada:
        n = N  # tamanho da matriz A
        L = identidade(n)
        
        for k in range(0, n-1):
            #Para cada linha i
            for i in range(k+1,n):
                #Calcula o fator m de gauss
                m = - A[i][k]/ A[k][k]
                L[i][k] = -m

                #Atualiza a linha i da matriz percorrendo todas as colunas j
                for j in range(k+1, n):
                    A[i][j] = m* A[k][j] + A[i][j]                
                A[i][k] = 0

     
        # Calculando Ly = B
        y = substituicoes_sucessivas(L, B, N)

        #Calculando Ux = B
        x = substituicoes_retroativas(A, y, N)
        print(arredondar(x,2)) 
        if (DET == '1'):
            print(np.linalg.det(A))
        return x
   
    else: 
        print("Insira uma matriz quadrada")   


def cholesky(A, B, N, DET):

    n = len(A)
    L = np.zeros((n, n))
    # Conferindo se é simétrica e positiva
    quadrada = matrizQuadrada(A)
    positiva = matrizPositiva(A)

    if (quadrada & positiva) == False:
        print("Insira uma matriz válida")
        return (L, False) # retorna a matriz vazia e false indicando que a operação não foi realizada

    for i in range(n):
      for j in range(i+1):
        if (i == j):
          soma = 0
          for k in range(i):
            soma += L[i][k]**2
          L[i][i] = (A[i][i] - soma) ** 0.5
        else:
          soma = 0
          for k in range(i):
            soma += L[i][k] * L[j][k]
          L[i][j] = ((A[i][j] - soma) / L[j][j]) 

    Lt = matrizTransposta(L)
     # Calculando Ly = B
    y = substituicoes_sucessivas(L, B, N)

    #Calculando Ux = B
    x = substituicoes_retroativas(A, y, N)
    print(arredondar(x, 2)) 
    if (DET == '1'):
        print(np.linalg.det(A))
    return (x, True) # retorna a matriz L triangular inferior e True indicando que a operação foi realizada


def jacobi(A, b, tolerancia, DET):
    n = np.shape(A)[0] # tamanho da matriz
    x = len(A) * [1.0] # vetor solução, por enquanto definido como um vetor de um
    x_aux = len(A) * [1.0] # vetor auxiliar, para alocar as respostas das operações
    iteracao = 0 # contador de iterações
    res = 1.0 # resíduo inicial
    if condicao_convergencia(A):
        while res > tolerancia:
            for i in range(n): # percorre as linhas da matriz
                soma = 0.0
                for j in range(n): # percorre as colunas da matriz
                    if i != j: #verifica se não é a diagonal principal
                        soma += A[i][j]*x[j] # usa o vetor que não sofre alteração para completar a soma
                x_aux[i] = (b[i] - soma) / A[i][i]
            iteracao += 1
            res = residuo_vetor(A, x, x_aux) # calcula o residuo a cada iteração
            print(res)
            x = np.copy(x_aux) # atualiza o vetor x da solução mais aproximada
        print(arredondar(x,2), iteracao)
        if (DET == '1'):
            print(np.linalg.det(A))
        return x, iteracao
    else:
        print("Não converge.")


def gauss_seidel(A, b, tolerancia, DET): # matriz A; vetor resultado b; tolerancia de aproximação, por exemplo 10**-3
    n = len(A) # tamanho da matriz
    x = len(A) * [1.0] # vetor solução, por enquanto definido como um vetor de um
    x_aux = len(A) * [1.0] # vetor auxiliar, para alocar as respostas das operações
    iteracao = 0 # contador de iterações
    res = 1 # resíduo inicial
    if condicao_convergencia:
        while (res > tolerancia):
            for i in range(n): # percorre as linhas da matriz
                soma = 0
                for j in range(n): # percorre as colunas da matriz
                   if (i != j): #verifica se não é a diagonal principal
                        soma += A[i][j]*x_aux[j] # usa sempre o vetor com os valores mais atualizados para a soma
                x_aux[i] = (b[i] - soma) / A[i][i]
            iteracao += 1
            res = residuo_vetor(A, x, x_aux) # calcula o residuo a cada iteração
            x = np.copy(x_aux) # atualiza o vetor x da solução mais aproximada
        print(arredondar(x,2), iteracao, res)
        if (DET == '1'):
            print(np.linalg.det(A))
        return (arredondar(x,2), iteracao)
           
    else:
        print("A matriz não converge.")


if __name__ == '__main__':
    main()
