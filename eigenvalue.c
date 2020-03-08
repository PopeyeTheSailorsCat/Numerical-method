#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable:4996)
#define delta 1e-15
void ReadFromInput(int size,long double**matrix, long double *h) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      scanf("%lf", &matrix[i][j]);
    }
  }
  //for (int i = 0; i < size; i++) {
  scanf("%lf", h);
  //printf("%lf", *h);
  //}
}
long double InfiniteNorm(int size, long double *vector) {//Find an infinite norm
  int max = 0;
  for (int i = 1; i < size; i++) {
    if (vector[max] < vector[i]) {
      max = i;
    }
  }
  return vector[max];
}
void InitXY(int size, long double** matrix, long double *x, long double* y) {//Create the initial vector from the available data
  for (int i = 0; i < size; i++) {
    y[i] = matrix[i][i];
    }
  long double norm;
  norm = InfiniteNorm(size, y);
  for (int i = 0; i < size; i++) {
    x[i] = y[i]/norm;
  }
}

void ReinitX(int size, long double* y, long double* x) {//Recalculation of the vector
  long double norm = InfiniteNorm(size, y);
  for (int i = 0; i < size; i++) {
    x[i] = y[i] / norm;
  }
}

void Ax(int size, long double** matrix, long double *x, long double *y) {//Multiplying a matrix by a vector
  //long double* result = malloc(size * sizeof(long double));
  long double stepResult;
  for (int i = 0; i < size; i++) {
    stepResult = 0;
    for (int j = 0; j < size; j++) {
      stepResult += matrix[i][j] * x[j];
    }
    y[i] = stepResult;
  }
  //return result;
}

void InitHOld(int size, long double* HOld) {//Creating an initial approximation
  //for (int i = 0; i < size; i++) {
    *HOld = -10 + rand() % 20;
  //}
}
long double MaxDifference(int size, long double* x, long double* y) {//calculating the difference
  long double max;
 // max = x[0] - y[0];
  //printf("point");
  //for (int i = 0; i < size; i++) {
 //   if (fabs(x[i] - y[i]) > max); {
//      max = fabs(x[i] - y[i]); 
 //   }

  //}
  //printf("mpoint");
  return fabs(*x-*y);
}

int CalculateH(int size, long double* hold, long double* hnew, long double * x, long double *y, long double eps) {//calculating the coefficient
  //for (int i = 0; i < size; i++) {
  int k = 1;
    *hnew = y[k] / x[k];
    if (fabs(x[k]) < delta) {
      printf("happen xi<delta\n");
    }
  //}
  //printf("mpoint");
  if (MaxDifference(size, hold, hnew) < eps) {
    return 0;
  }
  else {
    return 1;
  }

}

long double ResultH(int size, long double * vector) {
  long double sum = 0;
  for (int i = 0; i < size; i++) {
    sum += vector[i];
  }
  return(sum / size);
}
void ChangeH(int size, long double* New, long double* Old) {
  *Old = *New;
}

void DiffAxHx(int size, long double**matrix, long double h, long double* x) {//The test data
  long double* result = malloc(size * sizeof(long double));
  Ax(size, matrix, x, result);
  for (int i = 0; i < size; i++) {
    x[i] = h*x[i];
  }
  printf("norm(Ax-hx) %g\n", fabs(InfiniteNorm(size, result) - InfiniteNorm(size, x)));
  free(result);
}


void PMmax(int size, long double** matrix, long double answer) {//using the method with different accuracy
  long double* vectorX = malloc(size * sizeof(long double));
  long double* vectorY = malloc(size * sizeof(long double));
  //long double* hold = malloc(size * sizeof(long double));
  //long double* hnew = malloc(size * sizeof(long double));
  long double hold;
  long double hnew;
  long double mineResult;
  int iter;
  for (long double eps = 1e-1; eps >= 1e-12; eps /= 10) {
    iter = 0;
    InitXY(size, matrix, vectorX, vectorY);
    Ax(size, matrix, vectorX, vectorY);
    InitHOld(size, &hold);
    while (CalculateH(size, &hold, &hnew, vectorX, vectorY,eps)) {
      ReinitX(size, vectorY, vectorX);
      Ax(size, matrix, vectorX, vectorY);
      ChangeH(size, &hnew, &hold);
      iter++;
    }
    mineResult = hnew;// ResultH(size, hnew);
    printf("eps = %g\n", eps);
    printf("iter %i\n", iter);
    printf("%.15f difmlb %g  ", mineResult, fabs(mineResult-answer));
    DiffAxHx(size, matrix, mineResult, vectorX);
  }
  free(vectorX);
  free(vectorY);
 // free(hold);
 // free(hnew);
}

int main(void) {
  srand(time(NULL));
  FILE * fp;
  FILE * fe;
  fp = freopen("input.txt", "r", stdin);
  int size;
  fe = freopen("output.txt", "w", stdout);
  scanf("%i", &size);
  long double **matrix = malloc(size * sizeof(long double*));
  long double answer;
  for (int i = 0; i < size; i++) {
    matrix[i] = malloc(size * sizeof(long double));
  }
  ReadFromInput(size, matrix, &answer);
  PMmax(size, matrix, answer);
  for (int i = 0; i < size; i++) {
    free(matrix[i]);
  }
  free(matrix);
  getch();
  return 0;
}