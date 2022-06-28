# o delta é delta = L/(numP-1)
# e o L tem no slide, que é b-a
# Site com tabela: https://gcpeixoto.github.io/ipynb-lab-metodos-numericos/aula-19-quadratura-gaussiana.html#formulas-com-mais-pontos

def wiValues():
    return {
        1: {
        1: 2.0,
        },
        2: {
        1: 1.0,
        2: 1.0,
        },
        3: {
        1: 0.55555556,
        2: 0.88888889,
        3: 0.55555556,
        },
        4: {
        1: 0.34785485,
        2: 0.65214515,
        3: 0.65214515,
        4: 0.34785485,
        },
        5: {
        1: 0.23692689,
        2: 0.47862867,
        3: 0.56888889,
        4: 0.47862867,
        5: 0.23692689,
        },
        6: {
        1: 0.17132449,
        2: 0.36076157,
        3: 0.46791393,
        4: 0.46791393,
        5: 0.36076157,
        6: 0.17132449,
        },
        7: {
        1: 0.12948497,
        2: 0.27970539,
        3: 0.38183005,
        4: 0.41795918,
        5: 0.38183005,
        6: 0.27970539,
        7: 0.12948497,
        },
        8: {
        1: 0.10122854,
        2: 0.22238103,
        3: 0.31370665,
        4: 0.36268378,
        5: 0.36268378,
        6: 0.31370665,
        7: 0.22238103,
        8: 0.10122854,
        },
    9: {
        1: 0.08127439,
        2: 0.18064816,
        3: 0.2606107 ,
        4 :0.31234708,
        5: 0.33023936,
        6: 0.31234708,
        7: 0.2606107 ,
        8: 0.18064816,
        9: 0.08127439,
        },
        10: {
        1: 0.06667134,
        2: 0.14945135,
        3: 0.21908636,
        4: 0.26926672,
        5: 0.29552422,
        6: 0.29552422,
        7: 0.26926672,
        8: 0.21908636,
        9: 0.14945135,
        10: 0.06667134,
        },
    }


def xiValues():
  return  {
    1: {
      1: 0.0,
    },
    2: {
      1: -0.5773502691,
      2: +0.5773502691,
    },
    3: {
      1: -0.77459667,
      2: 0.0,
      3: 0.77459667,
    },
    4: {
      1: -0.86113631,
      2: -0.33998104,
      3: 0.33998104,
      4: 0.86113631,
    },
    5: {
      1: -0.90617985,
      2: -0.53846931,
      3: 0.0,
      4: 0.53846931,
      5: 0.90617985,
    },
    6:{
      1: -0.93246951,
      2: -0.66120939,
      3: -0.23861919,
      4:  0.23861919 ,
      5:  0.66120939,
      6:  0.93246951,
    },
    7: {
      1: -0.94910791 ,
      2: -0.74153119,
      3: -0.40584515,
      4: 0,
      5: 0.40584515,
      6: 0.74153119,
      7: 0.94910791 ,
    },
    8: {
      1: -0.96028986,
      2: -0.79666648,
      3: -0.52553241,
      4: -0.18343464 ,
      5: 0.18343464 ,
      6: 0.52553241,
      7: 0.79666648,
      8: 0.96028986,
    },
    9: {
      1: -0.96816024,
      2: -0.83603111,
      3: -0.61337143,
      4: -0.32425342,
      5: 0.0,
      6: 0.32425342,
      7: 0.61337143,
      8: 0.83603111,
      9: 0.96816024,
    },
    10: {
      1: -0.97390653,
      2: -0.86506337,
      3: -0.67940957,
      4: -0.43339539,
      5: -0.14887434,
      6: 0.14887434,
      7: 0.43339539,
      8: 0.67940957,
      9: 0.86506337,
      10: 0.97390653,
    },
  }
