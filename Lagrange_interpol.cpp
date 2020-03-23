#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime> 
#define Pi 3.14159265           
using namespace std;


void ChebKnot(double a, double b, int dote, double* resultDotes) { //Borders+ number of points+ pointer to an array for points
  double t;
  for (int i = 0; i < dote; i++) {
    t = cos(3.14159265* (2*(i + 1) +1) /(2* (dote + 1)));
    resultDotes[i] = (a + b) / 2 + t*(b - a) / 2;
  }
  return;
}

void Ravn(double a, double b, int dote, double* resultDotes) {//Uniform grid
  for (int i = 0; i < dote; i++) {
    resultDotes[i] = a + (b - a) / dote*i;
  }
}

void HalfStep(double a, double b, int dote, double* resultDotes) {
  double step = (b - a) / 4;//Half way of half segment
  cout << "get into "<< dote  << endl;
  resultDotes[0] = a;//put the beginning of the segment as the first point
  resultDotes[dote - 1] = b;
  for (int i = 1; i <= dote - 2 / 2; i++) {
    resultDotes[i] = resultDotes[i - 1] + step;
    cout << resultDotes[i] << " ";
    resultDotes[dote - 1 - i] = resultDotes[dote - i] - step;
    cout << resultDotes[dote - 1 - i] << endl;
    step /= 2;
  }
 
  if (dote % 2 == 1) {
    resultDotes[dote / 2 + 1];
  }
  return;
}

void RandomDotes(double a, double b, int dote, double* resultDotes) {
  for (int i = 0; i < dote+1; i++) {
    resultDotes[i]= 0.3 + 0.0001 * (rand() % 8001);// 0.3+0.8
  }
  return;
}

double MyFunction(double x) {
    return(log(1 + x*x));
  //return sin(x)+cos(x);
}

void MultiplyPolynpom(int dote, double* dotes, double* ResultPolynom) {//Multiplication of polynomials is not being used but works
  double *workPolynom = new double[dote+1];//to save intermediate results
  double *XworkPolynom = new double[dote+1];//For a polynomial*on x
  double *CworkPolynom = new double[dote+1];//For a polynomial to a constant
  double x= 1; //to store the denominator
  for (int i = 0; i < dote; i++) {
    int realj = 1; //Due to the problem of jumping the step where i=j, because of it in the final results there is an offset, I use in polydromes
    workPolynom[1]=1; 
    XworkPolynom[0] = 0;    //In order to get the first step and below too
    if (i != 0) {
      workPolynom[0] = -dotes[0];
    }
    else {
      workPolynom[0] = -dotes[1];
    }
    //In order to get the first step
    x = 1;
    for (int j = 0; j < dote; j++) {
      if (i != j ) {

        x *= (dotes[i] - dotes[j]);
        if ((i == 0 && j != 1) || (i!=0 && j!=0)) {
          for (int n = 0; n < realj; n++) {//Multiplication of polynomials
            XworkPolynom[n + 1] = workPolynom[n];//Elements that were multiplied by x increased the power of
            CworkPolynom[n] = (-dotes[j]) * workPolynom[n]; //The rest were multiplied by subtracting xj
          }
          for (int n = 0; n < realj; n++) {
            workPolynom[n] = XworkPolynom[n] + CworkPolynom[n];
          }
          workPolynom[realj] = XworkPolynom[realj];
          realj++;
        }
        else {
          realj++;
        }
        
      }
    }
    for (int n = realj-1; n >=0; n--) {
      
      ResultPolynom[n] += MyFunction(dotes[i]) * workPolynom[n] / x;
    }
  }
  delete[] workPolynom;
  delete[] XworkPolynom;
  delete[] CworkPolynom;
}

double CountDotes(int dote, double* dotes) {
  double result=0;
  double dif=0;
  for (int i = 1; i < dote; i++) {
    result = 0;
    double currentDote = (dotes[i - 1] + dotes[i]) / 2;
    for (int j = 0; j < dote; j++) {
      double stepresult = 1;
      for (int k = 0; k < dote; k++) {
        if (j != k) {
          stepresult = stepresult* (currentDote - dotes[k]) / (dotes[j] - dotes[k]);
         
        }
      }
      result += stepresult*MyFunction(dotes[j]);
    }
    if (fabs(result - MyFunction(currentDote)) > dif) {
      dif = fabs(result - MyFunction(currentDote));
    }
  }
  return dif;
}

void LagrangeInterpol(double a, double b) {
  ofstream outfiler;
  outfiler.open("out.txt", ios::out);
  for (int dote = 3; dote < 30; dote++) {
    double *dotes = new double[dote];
    double *resultPolynom = new double[dote+1]();
    ChebKnot(a, b, dote, dotes);
    outfiler << CountDotes(dote, dotes) << endl;//<<' ' <<dote<< endl;
    delete[] dotes;
    delete[] resultPolynom;
  }
  outfiler << endl << endl;

  outfiler << "random " << endl;
  for (int dote = 3; dote < 30; dote++) { //another grid
    double *dotes = new double[dote];
    double srednDif=0;
    HalfStep(a, b, dote, dotes);

    for (int i = 0; i < 10; i++) {
      srednDif += CountDotes(dote, dotes);
    }
    outfiler << srednDif / 10 << endl;
    delete[] dotes;
  }

}
void main(void) {
  double a = 0.3;
  double b = 1.1;
  LagrangeInterpol(a, b);
  int response;
  cin >> response;
}