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
    y0 = 0
    RungeKutta()


def f(t, a1, a2, a3, w1, w2, w3):
    return a1*math.sin(w1*t) + a2*math.sin(w2*t)+ a3*math.cos(w3*t)      

def RungeKutta(m, c, k, a1, a2, a3, w1, w2, w3, y0):
    t0 = 0.0
    h = 0.1
    tk = k*h

    for i in range(1,10):
        k1 = f()
        Y = y0 + 1/6(k1 + 4*k2 + k3)*h


