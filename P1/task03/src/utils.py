def multiplica_matriz(mat1, mat2, n):  
   mat3 = [[0 for j in range(n) ] for i in range(n)]
   for i in range(n):
     for j in range(n):
       for k in range(n): #calcula prod. interno da linha i por coluna j
         mat3[i][j] = mat3[i][j] + (mat1[i][k] * mat2[k][j])
   return mat3


