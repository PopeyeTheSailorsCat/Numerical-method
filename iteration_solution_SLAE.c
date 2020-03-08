#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#pragma warning(disable:4996)
int count_pm = 0;
int count_ud = 0;
int operation = 647;
void matrixFree(long double **matrix, long double **Lmatrix, long double *Dmatrix, int  size) {
  for (int i = 0; i < size; i++) {
    free(matrix[i]);
    free(Lmatrix[i]);
  }
  free(matrix);
  free(Lmatrix);
  free(Dmatrix);
}

long double Scalar(int size,long double* Vector, long double* secondVector) {//scalar multiplication of two vectors
  long double result=0;
  for (int i = 0; i < size; i++) {
    result += Vector[i] * secondVector[i];
    count_pm++;
    count_ud++;
  }
 // printf("sc %lf\n", result);
  return result;
}
void MatrixRead(int size, long double **matrix, long double * vector, long double * solution) {
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      scanf("%lf", &matrix[i][j]);
      //printf("%lf", matrix[i][j]);
    }
    //printf("\n");
  }
  for (int i = 0; i < size; i++) {
    scanf("%lf", &vector[i]);
    //printf("%lf", vector[i]);
  }
  for (int i = 0; i < size; i++) {
    scanf("%LF", &solution[i]);
    //printf("%lf", vector[i]);
  }
}
void FirstInitX(int size, long double** matrixA, long double* vectorB, long double*vectorX) {//creating a vector from initial data
  for (int i = 0; i < size; i++) {
    if (matrixA[i][i] != 0) {
      vectorX[i] = fabs(vectorB[i] / matrixA[i][i]);
      count_ud++;
    }
    else {
      vectorX[i] = fabs(vectorB[i]);
    }
    //printf("%lf\n", vectorX[i]);
  }

}
void InitVectorGrad(int size, long double** matrixA, long double* vectorB, long double* vectorX, long double* vectorG) {//creating a gradient vector
  long double Aix;
  for (int i = 0; i < size; i++) {
    Aix = 0;
    for (int j = 0; j < size; j++) {
      Aix += matrixA[i][j] * vectorX[j];
      count_pm++;
      count_ud++;
    }

    vectorG[i] = (Aix - vectorB[i]);
    count_pm++;
  }
 /* for (int i = 0; i < size; i++) {
    printf("%g ", vectorDiff[i]);
  }
  getch();*/
}

_Bool CheckDiff(int size,long double* vectorGrad, long double epsi, long double coef) {//check out
  long double step;
  for (int i = 0; i < size; i++) {
    //step = 0;
    //for (int j = 0; j < size; j++) {
    //  step += matrix[i][j] * vectorX[j];
    //}
    if (fabs(vectorGrad[i]*coef) > epsi) {
      return 1;
    }
  }
  return 0;
}


long double InitCoef(int size, long double** matrixA, long double* vectorGrad) {//calculating the coefficient
  long double* vectorAg = calloc(size, sizeof(long double));
  long double step;
  for (int i = 0; i < size; i++) {
    step = 0;
    for (int j = 0; j < size; j++) {
      step += matrixA[i][j] * vectorGrad[j];
      count_pm++;
      count_ud++;
    }
    vectorAg[i] = step;
  }
  step = 2 * Scalar(size, vectorAg, vectorGrad);
  count_ud++;
  free(vectorAg);
  count_ud++;
  return (Scalar(size, vectorGrad, vectorGrad) / step );
}

void ReInitX(int size, long double* vectorX,  long double Coef, long double* vectorGrad) {//recalculation of the vector
  for (int i = 0; i < size; i++) {
    vectorX[i] -= Coef*vectorGrad[i];
    count_pm++;
    count_ud++;
    //printf("%lf\n", vectorX[i]);

  }
  //printf("\n");
  //getch();
}
void ReInitGrad(int size, long double** matrixA, long double* vectorG, long double* vectorB, long double* vectorX) {//recalculating the gradient
  long double Ax;
  for (int i = 0; i < size; i++) {
    Ax = 0;
    for (int j = 0; j < size; j++) {
      Ax += matrixA[i][j] * vectorX[j];
      count_pm++;
      count_ud++;
    }
    vectorG[i] = (Ax - vectorB[i]);
    count_pm++;
    //printf("g %lf\n", vectorG[i]);
  }
}
void GradientMin(int size, long double** matrixA, long double* vectorB, long double* solutionMatlab) {//using the gradient method with different accuracy
  long double* vectorGrad = malloc(size * sizeof(long double));                                       //and counting the number of operations
  long double* vectorX = calloc(size, sizeof(long double));
  int iter;
  long double Coef;
  _Bool flag_1 = 0;
  _Bool flag_2 = 0;
  _Bool flag_3 = 0;
  for (long double epsi = 1e-1; epsi >= 1e-12; epsi /= 10) {
    iter = 0;
   
    FirstInitX(size, matrixA, vectorB, vectorX);//создание вектора х
    ReInitGrad(size, matrixA, vectorGrad, vectorB, vectorX);
    Coef = InitCoef(size, matrixA, vectorGrad);
    while(CheckDiff(size, vectorGrad, epsi,Coef)){
      //printf("point\n");
      if (((float)operation / 2)<(count_pm + count_ud) && flag_1 == 0) {
        printf("n/2  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_1 = 1;
      }
      if (operation<(count_pm + count_ud) && flag_2 == 0) {
        printf("n  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_2 = 1;
      }
      if (2 * operation<(count_pm + count_ud) && flag_3 == 0) {
        printf("2n   %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_3 = 1;
      }
      Coef=InitCoef(size, matrixA, vectorGrad);
      if (((float)operation / 2)<(count_pm + count_ud) && flag_1 == 0) {
        printf("n/2  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_1 = 1;
      }
      if (operation<(count_pm + count_ud) && flag_2 == 0) {
        printf("n  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_2 = 1;
      }
      if (2 * operation<(count_pm + count_ud) && flag_3 == 0) {
        printf("2n   %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_3 = 1;
      }
      //printf(" coef %lf\n", Coef);
      ReInitX(size, vectorX, Coef, vectorGrad);//Пересчет x
      if (((float)operation / 2)<(count_pm + count_ud) && flag_1 == 0) {
        printf("n/2  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_1 = 1;
      }
      if (operation<(count_pm + count_ud) && flag_2 == 0) {
        printf("n  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_2 = 1;
      }
      if (2 * operation<(count_pm + count_ud) && flag_3 == 0) {
        printf("2n   %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_3 = 1;
      }
      ReInitGrad(size, matrixA, vectorGrad, vectorB, vectorX);
    // printf("%lf\n",vectorGrad[0]);
    // getch();
      if (((float)operation / 2)<(count_pm + count_ud) && flag_1 == 0) {
        printf("n/2  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_1 = 1;
      }
      if (operation<(count_pm + count_ud) && flag_2 == 0) {
        printf("n  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_2 = 1;
      }
      if (2 * operation<(count_pm + count_ud) && flag_3 == 0) {
        printf("2n   %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_3 = 1;
      }
      iter++;
      //printf("%i", iter);
      //printf("%i", count_pm + count_ud);
      if (((float)operation/2)<(count_pm+count_ud)&& flag_1==0){
        printf("n/2  %i   %.12f\n", iter ,fabs(vectorX[5] - solutionMatlab[5]));
        flag_1 = 1;
      }
      if (operation<(count_pm + count_ud) && flag_2 == 0) {
        printf("n  %i   %.12f\n", iter, fabs(vectorX[5] - solutionMatlab[5]));
        flag_2 = 1;
      }
      if (2*operation<(count_pm + count_ud) && flag_3 == 0) {
        printf("2n   %i   %.12f\n", iter ,fabs(vectorX[5] - solutionMatlab[5]));
        flag_3 = 1;
      }
    }
    //printf("%i", iter);
   // printf("%i", count_pm + count_ud);
   printf("iter %i\n", iter);
    printf("eps %g\n", epsi);
    printf("result\n :");
    for (int i = 5; i < 6; i++) {
      Coef = 0;
      for (int j = 0; j < size; j++) {  
        Coef += matrixA[i][j] * vectorX[j];
      }
      printf("%.12f   difmatlab:%g  Ax-b %g\n ", vectorX[i],fabs( vectorX[i] - solutionMatlab[i]), Coef-vectorB[i] );
    }
    printf("\t\t+- %i\n", count_pm);
    printf("\t\t/ %i\n", count_ud); 
    //printf("\t\t+- %i\n", count_pm);
  //  printf("\t\t+- %i\n", count_pm);
  //  printf("\t\t/ %i\n", count_ud);
    count_pm = 0;
    count_ud = 0;
  }
  free(vectorGrad);
  free(vectorX);
}
int main(void) {
  int  size;
  //scanf("%lf", &xStart);
  
  FILE* fp;
  FILE * Fp;
  fp = freopen("input.txt", "r", stdin);
  Fp = freopen("output.txt", "w", stdout);
  scanf("%i", &size);
 // printf("%i\n", size);
  long double **matrixA = malloc(size*sizeof(long double*));
  long double * vectorB = malloc(size*sizeof(long double));
  long double * solutiomMatlab = malloc(size * sizeof(long double));
  for (int i = 0; i < size; i++) {
    matrixA[i] = malloc(size * sizeof(long double));
  }
  MatrixRead(size, matrixA, vectorB, solutiomMatlab);
  //close(fp);
  GradientMin(size, matrixA, vectorB, solutiomMatlab);
  getch();
  return 0;
}