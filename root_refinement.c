#include <stdio.h>
#include <math.h>
#pragma warning(disable:4996)

typedef long double(*function)
(long double x);

typedef long double(*dirfunction)
(long double x);


long double funct(long double x) {//the function
  return (pow(x, 4) - 18 * pow(x, 3) - 10);
}

long double dirfunct(long double x) {//the function dir
   return(4 * pow(x, 3) - 54 * pow(x, 2));
}

long double seconddirfucnt(long double x) {//the function second dir
  return(12 * pow(x, 2) - 108 * x);
}

long double seconddirfunct2(long double x) {//the function2 second dir
  return -(2*x*(3*pow(x,7)+3*pow(x,4)-x*x*x-6*x+1))/pow((pow(x,6)-2*x*x*x+x*x+1),2);
}

long double dirfunct2(long double x) {//the function2 dir
  return -(pow(x, 2)*(pow(x, 4) - 4 * x + 1) / (pow(x, 6) - 2 * pow(x, 3) + pow(x, 2) + 1));
}

long double funct2(long double x) {////the function2
  return (atan(pow(x,2)+pow(x,-1)) - x);
}

long double minf(function functinside, long double left, long double right) {//minimum of function
  double iRun, min;
  iRun = left;
  min = functinside(iRun);
  while(iRun <= right) {
    if (min > functinside(iRun)) {
      min = functinside(iRun);
    }
    iRun += 0.2;
  }
  return min;
}

long double maxf(function functinside, long double left, long double right) {//maximum of function
  double iRun, max;
  iRun = left;
  max = functinside(iRun);
  while (iRun <= right) {
    if (max < functinside(iRun)) {
      max = functinside(iRun);
    }
    iRun += 0.2;
  }
  return max;
}

/*int printMass(int x[20]) {
  int iRun;
  printf("[");
  for (iRun = 1; iRun < 20; iRun++) {
    printf("10^%i; ", x[iRun]);
  }
  printf("]\n");
}*/

halfelem(long double mainLeft, long double mainRight, function functinside) {//The method of half division
  int iRun, iRuntoo;
  long double eps, mid, right, left;
  for (iRun = -1; iRun > -20; iRun--) {//different accuracy
    eps = pow(10, iRun);
    iRuntoo = 0;
    left = mainLeft;
    right = mainRight;
    printf("eps=e%i    ", iRun);

    while ((right - left) / pow(2, iRuntoo) >= eps) {
      iRuntoo++;
      mid = (left + right) / 2;
      if (functinside(mid) <eps) {
        break;
      }
      if (functinside(left)*functinside(mid) < 0) {
        right = mid;
      }
      else {
        left = mid;
      }
    }
    mid = (left + right) / 2;
    //printfz("+");

    printf("%e  it %i \n", mid, iRuntoo);
  }
    printf("//////////////////////////\n");
}

modifiedNewton(long double mainLeft, long double mainRight, function functinside, dirfunction dirfun,function secondir, long double left, long double right) {
  long double dir,firstpoint, eps, xn,xnext, coef; //modified Newton method
  int whatTime, iRun, iRuntoo;
 // int iter[20];
//  int epsi[20];
  printf("Give first point:");
  scanf("%Lf", &firstpoint);
  coef = maxf(secondir, left, right / (2 * minf(dirfun, left, right)));
  
  for (iRun = -1; iRun > -20; iRun--) {//different accuracy
    eps = pow(10, iRun);
    iRuntoo = 0;
    printf("eps=e%i    ", iRun);
    dir = dirfun(firstpoint);
    xn = firstpoint;
    xnext = firstpoint - funct(firstpoint) / dir;

    while (fabs(xnext - xn) > coef*eps) {
      iRuntoo++;
      xn = xnext;
      xnext = xnext - functinside(xnext) / dir;
      //printfz("+");
    }
 //  epsi[fabs(iRun)] = eps;
  //  iter[fabs(iRun)] = iRuntoo;
    printf("%e  it %i \n", xnext, iRuntoo);

  }
    printf("//////////////////////////\n");
//    printMass(epsi);

}

int main(void) {
  long double mainLeft, mainRight;
  scanf("%Lf%Lf", &mainLeft, &mainRight);
  //printf("%e", funct(1));
  halfelem(mainLeft, mainRight, funct);
  halfelem(mainLeft, mainRight, funct2);
  modifiedNewton(mainLeft, mainRight, funct, dirfunct, seconddirfucnt, mainLeft, mainRight);
  modifiedNewton(mainLeft, mainRight, funct2, dirfunct2, seconddirfunct2, mainLeft, mainRight);
  scanf("%Lf", &mainLeft);
  return 0;
}