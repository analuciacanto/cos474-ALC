import numpy as np
from utils import *

def getMatrizN():
    pontos = []

    with open("task03/src/coord.dat") as file:
        for line in file.readlines():
            data = line.split()
            pontos.append(data)
    return pontos

def main(): 

    ICOD = input('Entre com o ICOD do método desejado: ')
    
    matrix = getMatrizN()
  
    coord = []

    for i in range(len(matrix)):
        coord.append(list(np.float_(matrix[i])))

    if (ICOD == '1'):
        interpolation(len(coord),coord, 1)
    elif (ICOD == '2'):
       regression( "2.0 + x", [[1.0,2.0], [2.0, 3.5], [3.0,6.5]], 1)
  
    else:
        print("Insira um código válido")

def interpolation( n,coord, x):
    
    coefL = []  
    y = 0
    for i in range(n):
      L = 1 
      num = 1
      deno = 1
      for j in range(n):
            if i != j:
                num *=  (x - coord[j][0])
                deno *= (coord[i][0] - coord[j][0])
      y += (num/deno) * coord[i][1]       

   
    print("y(" + str(x) +") = " + str(y))

def regression(f, coord, x):
    # Os coeficientes da função deverão ser omitidos, presume-se estarem em ordem: b1, b2, b3 ... 

    P = []
    Y = []

    fx = f.split(" + ")

    for i in range(len(coord)): #Número de pontos
        Y.append(coord[i][1]) 
        pLinha = []
        
        for j in range(len(fx)):   #Quantidade de termos (b)   
            fx = f.split(" + ")
            if ('x' in fx[j]):        
                fx[j] = fx[j].replace('x', str(coord[i][0]))          
            
            pLinha.append(float(fx[j]))   

        P.append(pLinha) #Montamos a matriz P

    print("P")
    print(P)

    print("Y")
    print(Y)
    #Calculando A(transposta) * A
    
    A = np.dot(np.transpose(P), P) 

    print(A)

    #Calculando C

    C = np.dot(np.transpose(P), Y)
    print("C")
    print(C)
    
    B = np.dot(np.linalg.inv(A), C)

    print("B")
    print(B)

if __name__ == '__main__':
    main()
