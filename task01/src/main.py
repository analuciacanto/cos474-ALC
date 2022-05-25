import numpy as np

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



def main(ICOD): 

    Amatrix = getMatrizA()
    Bmatrix = getMatrizB()

    A = []
    B = []

    for i in range(len(Amatrix)):
        A.append(list(np.float_(Amatrix[i])))
  
    B = list(np.float_(Bmatrix))    


    A = [[1,2,2],[4,4,2],[4,6,4]]
    B = [3,6,10]
    
    if (ICOD == 1):
        decomposicaoLU(A, B, len(A))
       

    else:
        print("Insira um código válido")


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



def decomposicaoLU(A, B, N):
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

        print(L)
        print(A)
        
        # Calculando Ly = B
        y = substituicoes_sucessivas(L, B, N)
        print(y)

        #Calculando Ux = B
        x = substituicoes_retroativas(A, y, N)
        print(x)
        return x
   
    else: 
        print("Insira uma matriz quadrada")   


main(1)