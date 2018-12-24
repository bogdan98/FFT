#include <iostream>
#include <fstream>
#include <string>
#include <windows.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
using namespace std;

double pi = 3.14159265359;

//Implementation of Complex Number class with relevant operations for complex numbers
class Cmplx
{
public:

  double re;
  double im;

  double mod(){
    return sqrt(pow(re, 2.0) + pow(im, 2.0));
  }

  double arg(){
    double res;
    if ((re == 0.0) && (im == 0)){
      res = 0.0;
    }
    else if ((re > 0.0) && (im > 0.0)){
      res = atan(im/re);
    }
    else if ((re > 0.0) && (im < 0.0)){
      res = 2*pi + atan(im/re);
    }
    else if ((re < 0.0) && (im > 0.0)){
      res = pi + atan(im/re);
    }
    else if ((re < 0.0) && (im < 0.0)){
      res = pi + atan(im/re);
    }
    else if ((re == 0) && (im > 0)){
      res = pi/2.0;
    }
    else if ((re == 0) && (im < 0)){
      res = 1.5*pi;
    }
    else if ((re > 0) && (im == 0)){
      res = 0.0;
    }
    else if((re < 0) && (im == 0)){
      res = pi;
    }
    return res;
  }

  Cmplx mult(Cmplx c2){
    Cmplx res;
    res.re = re*c2.re - im*c2.im;
    res.im = re*c2.im + im*c2.re;
    return res;
  }

  Cmplx div(Cmplx c2){
    Cmplx res;
    res.re = (re*c2.re + im*c2.im)/pow(c2.mod(), 2.0);
    res.im = ( - re*c2.im + im*c2.re)/pow(c2.mod(), 2.0);
    return res;
  }

  Cmplx add(Cmplx c2){
    Cmplx res;
    res.re = re + c2.re;
    res.im = im + c2.im;
    return res;
  }

  Cmplx scale(double c){
    Cmplx res;
    res.re = re*c;
    res.im = im*c;
    return res;
  }

  void display(){
    if (im >= 0.0){
      printf("%f + %fi\n", re, im);
    }
    else {
      printf("%f - %fi\n", re, abs(im));
    }
  }

};

//Function that returns a Cmplx object provided modulus and argument
Cmplx polar(double mod, double arg){
  Cmplx res;
  res.re = mod*cos(arg);
  res.im = mod*sin(arg);
  return res;
}

//create an array with REAL values of signal(input as a real function of time) for later use in Fourier transform
double* signal (double dt, int N, double(*f)(double)){
  double* result = new double[N];
  for (int i = 0; i < N; i++){
    result[i] = f(dt*i);
  }
  return result;
}

// generalization of the previous function
Cmplx* signalg (double dt, int N, Cmplx(*f)(double)){
  Cmplx* result = new Cmplx[N];
  for (int i = 0; i < N; i++){
    result[i] = f(dt*i);
  }
  return result;
}

//some function of the input signal
double f(double x){
  return cos(5000*x);
}

//add two complex numbers
Cmplx add(Cmplx c1, Cmplx c2){
  Cmplx res;
  res.re = c1.re + c2.re;
  res.im = c1.im + c2.im;
  return res;
}

//sum of elements in an array of complex numbers
Cmplx sum(Cmplx *arr){
  Cmplx res;
  res.im = 0.0;
  res.re = 0.0;
  int N = sizeof(arr);
  for (int i = 0; i < N; i++){
    res.im = res.im + (arr[i]).im;
    res.re = res.re + (arr[i]).re;
  }
  return res;
}

// Fast Fourier transform: F(w) = int(-inf to inf) [f(t)*exp(-iwt)] - Cambridge convention.
// this becomes sum[f(t_j)*exp(-i*w*t_j)]; t_j =t_0 + dt*j; assume t_0 is zero, if necessary,
// the result can then be multiplied by t_0;
Cmplx Fourier(Cmplx *f, int N, double w, double dt){
  Cmplx result;
  int half = N/2;
  int rem = N%2;
  Cmplx *evenf = new Cmplx[half + N%2]; //starting from 0
  Cmplx *oddf = new Cmplx[half];
  if (half == 0){
    result = (f[0]).mult(polar(1.0, 0.0));
  }
  else {
    for (int i = 0; i < half + rem; i++){
      evenf[i] = f[2*i];
    }
    for (int i = 0; i < half; i++){
      oddf[i] = f[2*i + 1];
    }
    result = add(Fourier(evenf, half + rem, 2.0*w, dt), (polar(1.0, -w*dt)).mult(Fourier(oddf, half, 2.0*w, dt)));
  }

  delete[] evenf;
  delete[] oddf;

  return result;
}

//now inverse FT, but no division by 2pi (yet)
Cmplx iFourier(Cmplx *f, int N, double t, double dw){
  Cmplx result;
  int half = N/2;
  int rem = N%2;
  Cmplx *evenf = new Cmplx[half + N%2]; //starting from 0
  Cmplx *oddf = new Cmplx[half];
  if (half == 0){
    result = (f[0]).mult(polar(1.0, 0.0));
  }
  else {
    for (int i = 0; i < half + rem; i++){
      evenf[i] = f[2*i];
    }
    for (int i = 0; i < half; i++){
      oddf[i] = f[2*i + 1];
    }
    result = add(iFourier(evenf, half + rem, 2.0*t, dw), (polar(1.0, t*dw)).mult(iFourier(oddf, half, 2.0*t, dw)));
  }
  delete[] evenf;
  delete[] oddf;
  return result;
}

//Fast Fourier transform using pts number of points
Cmplx *FFT(Cmplx *f, int N, int pts, double dt){
  Cmplx *res = new Cmplx[pts];
  double fmax = 0.5/dt;
  for (int i = 0; i < pts; i++){
    res[i] = Fourier(f, N, -2*pi*fmax + 4*pi*fmax*i/(pts-1), dt);
  }
  return res;
}

//given a REAL input sampled signal, uses FFT to calculate the power spectrum
double* PowerSpectrum(double *f, int N, int pts, double dt){
  double* res = new double[pts];
  Cmplx* f1 = new Cmplx[N];
  Cmplx* f2 = new Cmplx[pts];
  for (int i = 0; i < N; i++){
    (f1[i]).re = f[i];
    (f1[i]).im = 0.0;
  }
  f2 = FFT(f1, N, pts, dt);
  for (int i = 0; i < pts; i++){
    res[i] = pow((f2[i]).mod(), 2.0);
  }
  return res;
}

//Fast inverse Fourier transform using pts number of points
Cmplx *iFFT(Cmplx *f, int N, int pts, double dt){
  Cmplx *res = new Cmplx[pts];
  double fmax = 0.5/dt;
  double *w = new double[pts];
  for (int i = 0; i < pts; i++){
    w[i] = -2*pi*fmax + 4*pi*fmax*i/(pts-1);
    res[i] = (iFourier(f, N, w[i], dt)).div(polar(2*pi, 0.0));
  }
  return res;
}

//gives the value of true signal at time s(rectangular window), given
//sampling frequency and sampled signal at the input
//this is convolution of the sampled signal and sinc function
double TrueSignalval (double *vals, int N, double dt, double s){
  //assume x starts from zero
  double res = 0.0;
  for (int i = 0; i < N; i++){
    if (s == dt*i){
      res = res + vals[i];
    }
    else{
      res = res + vals[i]*sin(pi*(s - dt*i)/dt)/(pi*(s - dt*i)/dt);
    }
  }
  return res;
}

//calculates the true signal on the grid provided
double* TrueSignal(double *vals, double *grid, int N, int Ngrid, double dt){
  double *res = new double[Ngrid];
  for (int i = 0; i < Ngrid; i++){
    res[i] = TrueSignalval(vals, N, dt, grid[i]);
  }
  return res;
}

int main(){
  double* c1 = new double[1000];
  c1 = signal(0.001, 1000, f);
  double *c3 = new double[1000];
  c3 = PowerSpectrum(c1, 1000, 1000, 0.001);
  ofstream myfile;
	myfile.open("FFT2.txt");
	for (int i = 0; i<1000; i++){
    myfile << -2*pi*500.0 + 4*pi*500.0*i/999.0 <<" "<< c3[i] << endl;
	}
  //sample code that writes the power spectrum of signal given by f
  //into a text file

}
