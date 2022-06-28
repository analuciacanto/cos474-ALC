
from valoresIntegracaoPolinomial import *
from valoresQuadraturaGaussiana import *
import math
import sympy as sym

def main():

    ICOD = input('Entre com o ICOD do método desejado: ')

    c1 = 1
    c2 = 2
    c3 = -3
    c4 = 3
    a = 0
    b = 1
    epsilon = 5*10**-4
    maxIter = 20
    x0 = 1
    numeroPontos = 5
    deltax1 = 0.5
    deltax2 = 0.25
    p = 1


    if (ICOD == '1'):
        METODO = input(
            'Digite 1 para o método da bisseção e 2 para o método de Newton: ')
        if (METODO == "1"):
            metodoBissecao(f, c1, c2, c3, c4, a, b, epsilon, maxIter)
        elif(METODO == "2"):
            metodoNewton(f, fderivada, c1, c2, c3, c4, x0, epsilon, maxIter)

    elif (ICOD == "2"):
        METODO = input(
            'Digite 1 para o método da Integração Polinomial e 2 para o método da Quadratura Gaussiana: ')
        if (METODO == "1"):
             algoritmoIntegracaoPolinomial(f, c1, c2, c3, c4, a, b, numeroPontos)
        elif (METODO == "2"):
            algoritmoQuadraturaGaussiana(f, c1, c2, c3, c4, a, b, numeroPontos)

    elif (ICOD == "3"):
        METODO = input(
            'Digite 1 para o a derivada passo a frente, 2 para derivada passo para trás e 3 para derivada central: ')
        if (METODO == "1"):
             derivadaPassoParaFrente(f, c1, c2, c3, c4, deltax1, x0)
        elif (METODO == "2"):
            derivadaPassoParaTras(f, c1, c2, c3, c4, deltax1, x0)
        elif (METODO == "3"):
            derivadaCentral(f, c1, c2, c3, c4, deltax1, x0)

    elif (ICOD == "4"):
        extrapolacaoRichard(f, c1, c2, c3, c4, deltax1, deltax2, x0, p)
       
    else:
        print("Insira um código válido")


def f(x, c1, c2, c3, c4):
    #return c1**(c2*x) + c3*x**c4
    return math.sin(x)
    

def fderivada(x0, c1, c2, c3, c4):
    x = sym.Symbol('x')
    df = sym.diff(c1**(c2*x) + c3*x**c4)
    return df.subs(x, x0)


def algoritmoBissecao(f, c1, c2, c3, c4, a, b, epsilon, maxIter):
    # Inicializa as variáveis Fa e Fb com os valores de f(a) e f(b), respectivamente
    Fa = f(a, c1, c2, c3, c4)
    Fb = f(b, c1, c2, c3, c4)

    # Teste para saber se a função muda de sinal. Se não mudar, mostrar
    # mensagem de erro
    if (Fa * Fb > 0):
        print("Erro! A função não muda de sinal")
        return (False, None)

    # Inicializa tamanho do intervalo intervX usando a função abs, x e Fx
    intervX = abs(b - a)
    x = (a + b)/2
    Fx = f(x, c1, c2, c3, c4)

    # Mostra dados de inicialização
    print("-\t%e\t%e\t%e\t%e\t%e\t%e\t%e" % (a, Fa, b, Fb, x, Fx, intervX))

    # Testa se o intervalo já é do tamanho da precisão e retorna a raiz sem erros
    if (intervX <= epsilon):
        return (True, x)

    # Inicializa k
    k = 1

    while k <= maxIter:
        # Se a função não mudar de sinal entre a e x, então atualiza o a e Fa.
        # Senão, atualiza o b e Fb

        if Fa * Fx > 0:
            a = x
            Fa = Fx
        else:
            b = x
            Fb = Fx

        # Atualiza intervX, x, e Fx
        intervX = abs(b - a)
        x = (a + b)/2
        Fx = f(x, c1, c2, c3, c4)

        print("-\t%e\t%e\t%e\t%e\t%e\t%e\t%e" % (a, Fa, b, Fb, x, Fx, intervX))

        # Teste de critério de parada
        if (intervX <= epsilon):
            return (True, x)

        # Incrementa k
        k += 1

    # Se chegar aqui é porque o número máximo de iterações foi atingido
    # Mostrar uma mensagem de erro e retorna que houve erro e a última raiz encontrada
    print("ERRO! número máximo de iterações atingido.")
    return (False, x)


def metodoBissecao(f, c1, c2, c3, c4, a, b, epsilon, maxIter):
    (executou, raiz) = algoritmoBissecao(
        f, c1, c2, c3, c4, a, b, epsilon, maxIter)
    if executou == False:
        print("O Método da Bisseção retornou um erro.")
    if raiz is not None:
        print("Raiz encontrada: %s" % round(raiz, 3))


def algoritmoNewton(f, fderivada, c1, c2, c3, c4, x0, epsilon, maxIter):
    # Teste se x0 já é logo a raiz

    if (abs(f(x0, c1, c2, c3, c4)) <= epsilon):
        return(True, x0)

    # Escreva o valor da aproximação inicial
    print("k\t x0\t f(x0)\t f'(x0)\t f'(x0)\t")

    # Inicie as iterações (pode ser um for)

    for i in range(1, maxIter + 1):

        # Em cada iteração:
        # Calcule x1 a partir de x0
        x1 = x0 - (f(x0, c1, c2, c3, c4) / fderivada(x0, c1, c2, c3, c4))

    # Escreva os valores de k, x1, f(x1)
        print("-\t%e\t%e\t%e\t%e\t%e" %
              (i, x0, f(x0, c1, c2, c3, c4), fderivada(x0, c1, c2, c3, c4), x1))

    # Teste para o critério de parada usando módulo da função
        if (abs(f(x0, c1, c2, c3, c4)) <= epsilon):
            return(True, x1)

    # Atualize o valor de x0

        x0 = x1

    # Se atingir o número máximo de iterações mostra mensagem de erro e retorna
    # a última raiz encontrada
    print("Erro! Número máximo de iterações atingido")
    return(False, x1)


def metodoNewton(f, fderivada, c1, c2, c3, c4, x0, epsilon, maxIter):
    (executou, raiz) = algoritmoNewton(
        f, fderivada, c1, c2, c3, c4, x0, epsilon, maxIter)
    if executou == False:
        print("O Método de Newton retornou um erro.")
    if raiz is not None:
        print("Raiz encontrada: %s" % round(raiz, 3))


def algoritmoQuadraturaGaussiana(f, c1, c2, c3, c4, a, b, numeroPontos): 
    
    I = 0

    for i in range(numeroPontos):
        A = wiValues()        
        X = xiValues()   

        x = ((b - a)/2)*X[numeroPontos][i+1] + (a+b)/2 
                      
        I += A[numeroPontos][i+1]*f(x, c1, c2, c3, c4)
            
        
    print(I)
    return I



def algoritmoIntegracaoPolinomial(f, c1, c2, c3, c4, a, b, numeroPontos): 
    
    I = 0
    L = b-a

    for i in range(numeroPontos):
        W = wiValuesIP(L)        
        X = xiValuesIP(a,b, (L)/(numeroPontos-1))   
       
        I +=  W[numeroPontos][i+1] * f(X[numeroPontos][i+1], c1, c2, c3, c4)            
           
    print(I)
    return I

def derivadaPassoParaFrente(f, c1, c2, c3, c4, deltax, xa):
    derivada = (f(xa + deltax,  c1, c2, c3, c4) - f(xa,  c1, c2, c3, c4)) / deltax
    print("Derivada passo para frente: "  + str(derivada))
    return derivada

    
def derivadaPassoParaTras(f, c1, c2, c3, c4, deltax, xa):
    derivada = (f(xa,  c1, c2, c3, c4) - f(xa - deltax,  c1, c2, c3, c4)) / deltax
    print(derivada)
    return derivada

def derivadaCentral(f, c1, c2, c3, c4, deltax, xa):
    derivada = (f(xa + deltax,  c1, c2, c3, c4) - f(xa - deltax,  c1, c2, c3, c4)) / 2*deltax
    print(derivada)
    return derivada


def extrapolacaoRichard(f, c1, c2, c3, c4, deltax1, deltax2, xa, p):
    q = deltax1/deltax2

    d1 =  derivadaPassoParaFrente(f, c1, c2, c3, c4, deltax1, xa)
    d2 = derivadaPassoParaFrente(f, c1, c2, c3, c4, deltax2, xa)

    derivada = d1 + (d1 - d2)/(q**(-p) - 1)

    print("Extrapolação de Richard: " + str(derivada))
    return derivada


main()
