#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int count_pm=0;
int count_ud = 0;
#pragma warning(disable:4996)
void matrixFree(long double **matrix, long double **Lmatrix, long double *Dmatrix,int  size) {//Clear Matrix
  for (int i = 0; i < size; i++) {
    free(matrix[i]);
    free(Lmatrix[i]);
  }
  free(matrix);
  free(Lmatrix);
  free(Dmatrix);
}

void MatrixRead(int size, long double **matrix, long double * vector) { //Read Matrix and vector
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      scanf("%lf", &matrix[i][j]);
    }
    
  }
  for (int i = 0; i < size; i++) {
    scanf("%lf", &vector[i]);
  }
}

long double sumLD(int j, long double* d, long double** l) { //Step summing
  long double result=0;
  for (int k = 0; k < j; k++) {
    result += l[j][k] * l[j][k] * d[k];
    count_pm++;
    count_ud += 2;
  }
  return result;
}

long double sumLLD(int i,int j, long double*d, long double** l) {//Step Summing
  long double result=0;
  for (int k = 0; k < j; k++) {
    result += l[i][k] * l[j][k] * d[k];
    count_pm++;
    count_ud += 2;
  }
  return result;
}
void LDLt(long double **matrix, long double *d, long double **l, long double **lt, int size) {//Hole Method
  for (int j = 0; j < size; j++) {
    d[j] = matrix[j][j] - sumLD(j, d, l);
    count_pm++;
    l[j][j] = 1;
    for (int i = j+1; i < size; i++) {
      l[i][j] = (1 / d[j])*(matrix[i][j] - sumLLD(i, j, d, l));
      count_pm += 1;
      count_ud += 2;
    }
  }
}

void PrintMatrix(long double **matrix, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      printf("%g ", matrix[i][j]);
    }
    printf("\n");
  }
}

void PrintTranMatrix(long double** matrix, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      printf("%g ", matrix[j][i]);
    }
    printf("\n");
  }
}
void PrintdoubleMatrix(long double **matrix, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      printf("%g ", matrix[i][j]);
    }
    printf("\n");
  }
}
void PrintDiagonalMatrix(long double *matrix, int size) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      if (i == j) {
        printf("%g ", matrix[j]);
      }
      else {
        printf("0 ");
      }
    }
    printf("\n");
  }
}

void PrintAll(long double** matrix, long double * Dmatrix, long double** Lmatrix, int matrixSize) {//Print different Matrix
  printf("matrix i read\n");
  PrintMatrix(matrix, matrixSize);
  printf("\n");
  printf("L matrix\n");
  PrintdoubleMatrix(Lmatrix, matrixSize);
  printf("\n");
  printf("D matrix\n");
  PrintDiagonalMatrix(Dmatrix, matrixSize);
  printf("\n");
  printf("Lt matrix\n");
  PrintTranMatrix(Lmatrix, matrixSize);
  printf("\n");
  printf("+- %i\n */ %i\n", count_pm, count_ud);
}
void SolveSystem(int size, long double**Lmatrix, long double*Dmatrix, long double* vectorB, long double** matrix) {//Use LDLt to solve SLAE
  long double *y = malloc(size * sizeof(long double));
  long double *z = malloc(size * sizeof(long double));
  long double *answers = malloc(size * sizeof(long double));
  long double **LDL = malloc(size * sizeof(long double*));
  for (int i = 0; i < size; i++) {
    LDL[i] = calloc(size, sizeof(long double));
  }
  for (int i = 0; i < size; i++) {
    scanf("%lf", &answers[i]);
  }
  long double know;
  z[0] = vectorB[0];
  for (int i = 1; i < size; i++) {
    know = 0;
    for (int j = 0; j < i; j++) {
      know += Lmatrix[i][j] * z[j];
      count_pm++;
      count_ud++;
    }
    z[i] = vectorB[i] - know;
    count_pm++;
  }
  for (int i = 0; i < size; i++) {
    y[i] = z[i] / Dmatrix[i];
    count_ud++;
  }

  for (int i = size-1; i >=0; i--) {
    know = 0;
    for (int j = i+1; j < size; j++) {
      know += Lmatrix[j][i] * z[j];
      count_pm++;
      count_ud++;
    }
    z[i] = y[i] - know;
    count_pm++;
  }
  PrintAll(matrix, Dmatrix, Lmatrix, size);
  printf("result and difference with matlab:\n");
  for (int i = 0; i < size; i++) {
    printf("%.15f\tdif  %.15f\n", z[i],fabs(z[i]-answers[i]) );
  }
  //LDLtmatrix
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      know = 0;
      know += Lmatrix[j][i] * Dmatrix[i];
      count_pm++;
      count_ud++;
      LDL[i][j] = know;
      printf("%.12f ", know);
    }
   // printf("\n");
  }
 printf("a-LdDt\n");//print different info

   for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      know = 0;
      for (int z = 0; z < size; z++) {
        know += Lmatrix[i][z] * LDL[z][j];
        count_pm++;
        count_ud++;
      }
      printf("%g ", matrix[i][j] - know);
    }
    printf("\n");
  }
  printf("Ax-b\n");
  for (int i = 0; i < size; i++) {
    know = 0;
    for (int j = 0; j < size; j++) {
      know += matrix[i][j] * z[j];
    }
    printf("%.15f\n", fabs(vectorB[i] - know));
  }
  free(z);
  free(y);

}

int main(void) {
  int matrixSize;
  //int elem;
  FILE *fp = freopen("matrix.txt", "r", stdin);  
  scanf("%i", &matrixSize);
  long double **matrix = malloc(matrixSize * sizeof(long double *));
  long double **Lmatrix = malloc(matrixSize * sizeof(long double*));
  long double *Dmatrix = malloc(matrixSize * sizeof(long double));
  long double *Bvector = malloc(matrixSize * sizeof(long double));
  for (int i = 0; i < matrixSize; i++) {
    matrix[i] = calloc(matrixSize , sizeof(long double));
    Lmatrix[i] = calloc(matrixSize , sizeof(long double));
  }
  MatrixRead(matrixSize, matrix,Bvector);
  LDLt(matrix, Dmatrix, Lmatrix, Lmatrix, matrixSize);
  //PrintAll(matrix, Dmatrix,Lmatrix, matrixSize);
  SolveSystem(matrixSize, Lmatrix, Dmatrix,Bvector, matrix);
  matrixFree(matrix, Lmatrix,Dmatrix, matrixSize);
  free(Bvector);
  getch();  
  return 0;
}