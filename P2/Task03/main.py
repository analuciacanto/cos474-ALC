import math

def main():
    m=1 
    c=0.1
    k=2
    a1 = 1
    a2 = 2
    a3 = 1.5
    w1 = 0.05
    w2 = 1
    w3 = 2
    x = 0
    dx = 0
    h = 0.1 #passo
    T = 3.0 #tempo total de integração

    RungeKutta(m,c,k,a1,a2,a3,w1,w2,w3,x, dx, h, T)


def F(t, a1, a2, a3, w1, w2, w3):
    return a1*math.sin(w1*t) + a2*math.sin(w2*t)+ a3*math.cos(w3*t)      

def Y_derivate(t, c, dx, k, x, m, a1, a2, a3, w1, w2, w3 ):
        return (F(t, a1, a2, a3, w1, w2, w3) - c*dx - k*x)/m


def  RungeKutta(m,c,k, a1,a2,a3,w1,w2,w3, x, dx, h, T):
    #Definindo tamanho do passo
    t = 0
    Y = 0

    while(t <= T):
        k1 = (h/2)*Y_derivate(t, c, dx, k, x, m, a1, a2, a3, w1, w2, w3)
        k2 = (h/2)*Y_derivate(t+(h/2), c, dx + k1, k, x + (h/2)*(dx+0.5*k1), m , a1, a2, a3, w1, w2, w3)
        k3 = (h/2)*Y_derivate(t+(h/2), c, dx + k2, k, x + (h/2)*(dx+0.5*k1), m , a1, a2, a3, w1, w2, w3)
        k4 = (h/2)*Y_derivate(t+(h/2), c, dx + 2*k3, k, x + (h/2)*(dx+k3), m , a1, a2, a3, w1, w2, w3)      

        Y += (h/2)*(dx + (1/3)*(k1+k2+k3))
        dx += 1/3*(k1+2*k2+2*k3+k4)
        t += h

        print('tempo:' + str(round(t,2)), 'deslocamento:' + str(round(Y,4)), 'velocidade:' + str(round(dx,4)), 'aceleracao:' + str(round(Y_derivate(t, c, dx, k, x, m , a1, a2, a3, w1, w2, w3), 4)))

  

main()