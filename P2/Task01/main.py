
import numpy as np

def main():
    c2 = 1.0
    c3 = 0.0
    c4 = 0.0
    th1 = 0.75
    th2 = 6.5 
    TOL = 10**(-4)
    N = 10

    algoritmoNewton(c2,c3,c4, th1, th2, TOL, N)


def y(c2,c3,c4, th1, th2):
    Y = []
    y1 = 2*(c3**2) + c2**2 + 6*c4**2  - 1.0
    y2 = 8*c3**3 + 6*c3*(c2**2) + 36*c3*c2*c4 + 108*c3*c4**2 - th1
    y3 = 60*c3**4 + 60*(c3**2)*(c2**2) + 576*c3**2*c2*c4 + 2232*c3**2*c4**2 + 252*c4**2*c2**2 + 1296*c4**3*c2 + 3348*c4**4 + 24*c2**3*c4 + 3*c2 -th2
    Y.append(y1)
    Y.append(y2)
    Y.append(y3)
    return Y


def jacobianaY(c2,c3,c4, th1, th2):
    M = []
    a = c2
    b = c3
    c = c4
    d = th1
    y1 = [2*a, 4*b, 12*c]
    y2 = [(12*b)*(a + 3*c), 6*(a**2 + 4*b**2 + 6*a*c + 18*c**2),  36*b*(a + 6*c)]
    y3 = [3 + 72*a**2*c + 576*b**2*c + 24*a * (5* b**2 + 21*c**2), 
    24*(5*a**2 * b + 10 * b**3 + 48* a* b* c + 186* b *c**2 + 54* c**3),
    24*(a**3 + 24* a* b**2 + 21* a**2* c + 6* c* (31* b**2 + 27* b* c + 93 *c**2))
    ]
    M.append(y1)
    M.append(y2)
    M.append(y3)
    return M


def algoritmoNewton(c2,c3,c4, th1, th2, TOL, N):
    c = [c2,c3,c4] 

    k = 1

    while(k <= N):
        
        fx = y(c[0],c[1],c[2], th1, th2)
        jx = jacobianaY(c[0],c[1],c[2], th1, th2)
        
        delta = - np.linalg.inv(jx).dot(fx)
        c = [c[0] + delta[0], c[1] + delta[1], c[2] + delta[2]]

        print({"Iteração":  k , "c1,c2, c3": c })
     
        if (np.linalg.norm(delta) / np.linalg.norm(c) < TOL):
            print("Tolerância máxima atingida")
            print(c)
            break

        k += 1




main()
