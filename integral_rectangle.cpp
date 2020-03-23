#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime> 
using namespace std;

double MatlabResult = 0.329487402603441;
double MyFunction(double x) {
  return(log(1 + x*x));
}

double RectangeIntegral(double a, double b, int dote) {//Calculating the integral using the rectangle formula
  double h = (b - a) / dote;
  double sumOfFunction = 0;
  for (int i = 1; i <= dote; i++) {
    sumOfFunction += MyFunction(a + i*h - h / 2);
  }
  return h*sumOfFunction;
}

void IntegrateWithDifDotes(double a, double b) {//Calculating the integral with different specified accuracy
  ofstream outfiler;
  outfiler.open("out.txt", ios::out);
  /*for (int i = 4; i <= 1000; i++) {
    double res = RectangeIntegral(a, b, i);
    outfiler << i << ' ' << fabs(res - MatlabResult) << endl;
  }*/
  for (double eps = 1e-1; eps > 1e-16; eps /= 10) {//Different accuracy
    int n = 2;
    int iter = 0;
    double halfWay = -10000;
    double currentWay = 1;
    while ((fabs(currentWay - halfWay) / 3) > eps) {
      halfWay = currentWay;
      currentWay = RectangeIntegral(a, b, n);
      n *= 2;
      iter++;
    }
    //outfiler << fabs(currentWay - MatlabResult) << " eps " << eps<< endl;
    outfiler << iter << endl;
  }
  return;
}


void main(void) {
  double a = 0.3;
  double b = 1.1;
  IntegrateWithDifDotes(a, b);
  cin >> a;
}