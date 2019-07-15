#ifndef _HEADER_FILE_NAME_H
#define _HEADER_FILE_NAME_H

#include <string>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TMath.h"
#include <algorithm>
#include <math.h>
#include "Math/IFunction.h"
#include "InitialCondition.cpp"
#include "ExactSolution.cpp"

#endif

TFile *fileoutput;
TGraph *gr;


const double lvc = 1.;
const int GridS = 1000;

const int numberofoutput = 5;
//const double timest = 0.15;
const double Spstep = 0.15;
const double timest = 0.1;
bool exactsolution = true;



double lim = GridS*Spstep;

// Initial conditions
//initial conditions can be set in InitialCondition.cpp   DoEval()

//old way but still work
TF1 init1 = TF1("f1", "cos(1.0*(x-5))*exp(-(x-5)*(x-5)/(2*5.0))",0,lim+0.01);
TF1 init2 = TF1("f1", "1.*sin(1.0*(x-5))*exp(-(x-5)*(x-5)/2.) + (x-5)* cos(1.0*(x-5))*exp(-(x-5)*(x-5)/2.)",0,lim+0.01);


// Exact solutions

TF2 exactsol = TF2 ("f2","cos(1.0*(x-50-y))*exp(-(x-50-y)*(x-50-y)/(2*5.0))",0,lim+0.01,0,lim+0.01);//wave
TF2 dexactsol = TF2 ("df2","-(sin(1.0*(x-5-y)) + (x-5-y)*cos(1.0*(x-5-y)))*exp(-(x-5-y)*(x-5-y)/2.)",0,lim+0.01,0,lim+0.01);

struct Sol{
  vector<double> Solvector;
  double time;
};



/* Classes for mathematical functions and methods, i. e. functions for initial conditions, boundary conditions, probably some known exact solutions, form of nonliear potentialetc*/

// class for the functions that set the from of the potential, i. e. phi3 or sin-gordon etc
class Potential{
public:
  Potential(string s, double a, double g);
  double GetValue();
private:
  double A;
};

class HyperbolicEquation{
public:
  double c;
  double CourantFactor();
  double Timestep;
  double Spatialstep;
  int Gridsize;
  vector<Sol> Adjective(string method, InitialCondition initial1);
  vector<Sol> Wave(string method, string pot, double g, string boundary, InitialCondition initial1, InitialCondition initial2, double time);
private:
  vector<Sol> LaxWave3(string pot, double g, string boundary,InitialCondition initial1, InitialCondition initial2, double time);
  vector<Sol> LaxWavePeriodic3(string pot, double g, InitialCondition initial1, InitialCondition initial2, double time);
  vector<Sol> LaxWendloffWave3(string pot, double g, string boundary, InitialCondition initial1, InitialCondition initial2, double time);
  vector<Sol> LaxWendloffWavePeriodic3(string pot, double g, InitialCondition initial1, InitialCondition initial2, double time);
  vector<Sol> LeapfrogWave(string pot, double g, string boundary,InitialCondition initial1, InitialCondition initial2, double time);
  vector<Sol> LeapfrogWavePeriodic(string pot, double g, InitialCondition initial1, InitialCondition initial2, double time);
  vector<double> transform(InitialCondition f);
  vector<double> transformLeapfrog(InitialCondition f);
  vector<double> VectOfDer(InitialCondition f);
  vector<double> VectOfDerLeapfrog(InitialCondition f);
};



Potential::Potential(string s, double a, double g) {

  if(s == "0")
    {
      A=0;
    }
  else if (s == "phi3")
    {
      A =g*a*a*a;
    }
  // default
  else
    {
      cout << "Wrong potential specification";
      A=0;
    };
};
double Potential::GetValue(){
  return A;
};



vector<Sol> HyperbolicEquation::Wave(string method, string pot, double g, string boundary, InitialCondition initial1, InitialCondition initial2, double time){
  if (method == "Lax3")
    {
     return LaxWave3(pot, g, boundary, initial1, initial2, time);
    }
  else if (method == "LW3")
    {
      return LaxWendloffWave3(pot, g, boundary, initial1, initial2, time);
    }
  else if (method == "Leapfrog")
    {
            return LeapfrogWave(pot, g, boundary, initial1, initial2, time);
    }
  // default gives empty vetor
  else
    {
      cout << "Wrong method";
      return {};
    };
};
vector<double> HyperbolicEquation::transform(InitialCondition f){
  vector<double> temp= {};
  for (int i=0;i<Gridsize;i++){
    temp.push_back(f.DoEval(i*Spatialstep));
  };
  return temp;
}

vector<double> HyperbolicEquation::transformLeapfrog(InitialCondition f){ //double Gridsize for halpsteps
  vector<double> temp= {};
  for (int i=0;i<2 * Gridsize;i++){
    temp.push_back(f.DoEval(i*Spatialstep));
  };
  return temp;
}

vector<double> HyperbolicEquation::VectOfDer(InitialCondition f){
  vector<double> temp = {};
  double a;
  for (int i = 0; i<Gridsize;i++){
    if (i == 0) a=(f.DoEval(0.01*Spatialstep)-f.DoEval(0))/(0.01*Spatialstep);
    a = (f.DoEval(i*Spatialstep + 0.01*Spatialstep) - f.DoEval(i*Spatialstep - 0.01*Spatialstep))/(0.02*Spatialstep);
    if (i == Spatialstep) a=(f.DoEval(i*Spatialstep)-f.DoEval(i*Spatialstep - 0.01*Spatialstep))/(0.01*Spatialstep);
    temp.push_back(a);
  };
  return temp;
};

vector<double> HyperbolicEquation::VectOfDerLeapfrog(InitialCondition f){ //double Gridsize for halpsteps
  vector<double> temp = {};
  double a;
  for (int i = 0; i<2 * Gridsize;i++){
    if (i == 0) a=(f.DoEval(0.01*Spatialstep)-f.DoEval(0))/(0.01*Spatialstep);
    a = (f.DoEval(i*Spatialstep + 0.01*Spatialstep) - f.DoEval(i*Spatialstep - 0.01*Spatialstep))/(0.02*Spatialstep);
    if (i == Spatialstep) a=(f.DoEval(i*Spatialstep)-f.DoEval(i*Spatialstep - 0.01*Spatialstep))/(0.01*Spatialstep);
    temp.push_back(a);
  };
  return temp;
};

vector<Sol> HyperbolicEquation::LaxWave3(string pot, double g, string boundary, InitialCondition initial1, InitialCondition initial2, double time){
  if (boundary == "periodic")
    {
      return LaxWavePeriodic3(pot, g, initial1, initial2, time);
    }
  // default -- periodic
  else
    {
      return LaxWavePeriodic3(pot, g, initial1, initial2, time);
    };
};

vector<Sol> HyperbolicEquation::LaxWendloffWave3(string pot, double g, string boundary, InitialCondition initial1, InitialCondition initial2, double time){
  if (boundary == "periodic")
    {
      return LaxWendloffWavePeriodic3(pot, g, initial1, initial2, time);
    }
  // default -- periodic
  else
    {
      return LaxWendloffWavePeriodic3(pot, g, initial1, initial2, time);
    };
};

vector<Sol> HyperbolicEquation::LeapfrogWave(string pot, double g, string boundary, InitialCondition initial1, InitialCondition initial2, double time){
  if (boundary == "periodic")
    {
      return LeapfrogWavePeriodic(pot, g, initial1, initial2, time);
    }
  // default -- periodic
  else
    {
      return LeapfrogWavePeriodic(pot, g, initial1, initial2, time);
    };
};


// three-variable system
vector<Sol> HyperbolicEquation::LaxWavePeriodic3(string pot, double g, InitialCondition initial1, InitialCondition initial2, double time){
  vector<double> winit = VectOfDer(initial1);
  vector<double> uinit = transform(initial1);
  vector<double> vinit = transform(initial2);
  vector<Sol> uOuttemp = {};
  Sol ustemp;
  int it = int(time/Timestep + 0.5);
  double ut, wt, vt;
  vector<double> u, w, v;
  vector<double> utemp, vtemp, wtemp;
  for (int i=0;i<it;i++){
   utemp = {};
   vtemp = {};
   wtemp = {};
    if (i==0){
      u = uinit;
      v = vinit;
      w = winit;
    }
    // periodic boundary conditions
    else {
      for (int j=0;j<Gridsize;j++){
       	Potential Potent(pot,g,v.at(j));
	if (j>0){
	  ut = 0.5 * (u.at((j+1)%Gridsize) + u.at((j-1)%Gridsize)) + Timestep * v.at((j)%Gridsize);
	  vt = 0.5 *  (v.at((j+1)%Gridsize) + v.at((j-1)%Gridsize)) + (c*c*Timestep/(2*Spatialstep)) * (w.at((j+1)%Gridsize) - w.at((j-1)%Gridsize)) - Potent.GetValue();
	  wt = 0.5 * (w.at((j+1)%Gridsize) + w.at((j-1)%Gridsize)) + (Timestep/(2* Spatialstep)) * (v.at((j+1)%Gridsize) - v.at((j-1)%Gridsize));
	}
	else {
	  ut = 0.5 * (u.at((j+1)%Gridsize) + u.at((j-1)%Gridsize + Gridsize)) + Timestep * v.at((j)%Gridsize);
	  vt = 0.5 * (v.at((j+1)%Gridsize) + v.at((j-1)%Gridsize + Gridsize)) + (c*c*Timestep/(2* Spatialstep)) * (w.at((j+1)%Gridsize) - w.at((j-1)%Gridsize + Gridsize)) + Potent.GetValue();
	  wt = 0.5 * (w.at((j+1)%Gridsize) + w.at((j-1)%Gridsize + Gridsize)) + (Timestep/(2*Spatialstep)) * (v.at((j+1)%Gridsize) - v.at((j-1)%Gridsize + Gridsize));
	}
	utemp.push_back(ut);
	vtemp.push_back(vt);
	wtemp.push_back(wt);
      };
      u = utemp;
      v = vtemp;
      w = wtemp;
    };
    ustemp.Solvector = u;
    ustemp.time = i*Timestep;
    uOuttemp.push_back(ustemp);

  };
  cout << u.size() << "  " << u.at(7);
  return uOuttemp; //u;
};


vector<Sol> HyperbolicEquation::LaxWendloffWavePeriodic3(string pot, double g, InitialCondition initial1, InitialCondition initial2, double time){
  vector<double> winit = VectOfDer(initial1);
  vector<double> uinit = transform(initial1);
  vector<double> vinit = transform(initial2);
  vector<Sol> uOuttemp = {};
  Sol ustemp;
  int it = int(time/Timestep + 0.5);
  double ut, wt, vt;
  vector<double> u, w, v;
  vector<double> utemp, vtemp, wtemp;
  for (int i=0;i<it;i++){
   utemp = {};
   vtemp = {};
   wtemp = {};
    if (i==0){
      u = uinit;
      v = vinit;
      w = winit;
    }
    // periodic boundary conditions
    else {
      for (int j=0;j<Gridsize;j++){
       	Potential Potent(pot,g,v.at(j));
	if (j>0){
	  ut = u.at((j)%Gridsize) + Timestep * v.at((j)%Gridsize)+ ((Timestep*Timestep*c*c)/(2*Spatialstep*Spatialstep))*(u.at((j+1)%Gridsize) + u.at((j-1)%Gridsize) - 2 * u.at((j)%Gridsize));
	  vt = v.at((j)%Gridsize) + (c*c*Timestep/(2*Spatialstep)) * (w.at((j+1)%Gridsize) - w.at((j-1)%Gridsize)) - Potent.GetValue() + ((Timestep*Timestep*c*c)/(2*Spatialstep*Spatialstep))*(v.at((j+1)%Gridsize) + v.at((j-1)%Gridsize) - 2 * v.at((j)%Gridsize));
	  wt = w.at((j)%Gridsize) + (Timestep/(2* Spatialstep)) * (v.at((j+1)%Gridsize) - v.at((j-1)%Gridsize)) + ((Timestep*Timestep*c*c)/(2*Spatialstep*Spatialstep))*(w.at((j+1)%Gridsize) + w.at((j-1)%Gridsize) - 2 * w.at((j)%Gridsize));
	}
	else {
	  ut = u.at((j)%Gridsize) + Timestep * v.at((j)%Gridsize) + ((Timestep*Timestep*c*c)/(2*Spatialstep*Spatialstep))*(u.at((j+1)%Gridsize) + u.at((j-1)%Gridsize + Gridsize) - 2 * u.at((j)%Gridsize));;
	  vt = v.at((j)%Gridsize) + (c*c*Timestep/(2* Spatialstep)) * (w.at((j+1)%Gridsize) - w.at((j-1)%Gridsize + Gridsize)) + Potent.GetValue() + ((Timestep*Timestep*c*c)/(2*Spatialstep*Spatialstep))*(v.at((j+1)%Gridsize) + v.at((j-1)%Gridsize + Gridsize) - 2 * v.at((j)%Gridsize));;
	  wt = w.at((j)%Gridsize) + (Timestep/(2*Spatialstep)) * (v.at((j+1)%Gridsize) - v.at((j-1)%Gridsize + Gridsize)) + ((Timestep*Timestep*c*c)/(2*Spatialstep*Spatialstep))*(w.at((j+1)%Gridsize) + w.at((j-1)%Gridsize + Gridsize) - 2 * w.at((j)%Gridsize));
	};
	utemp.push_back(ut);
	vtemp.push_back(vt);
	wtemp.push_back(wt);
      };
      u = utemp;
      v = vtemp;
      w = wtemp;
    };
    ustemp.Solvector = u;
    ustemp.time = i*Timestep;
    uOuttemp.push_back(ustemp);

  };
  cout << u.size() << "  " << u.at(7);
  return uOuttemp; //u;
};


//Leapfrog begin

vector<Sol> HyperbolicEquation::LeapfrogWavePeriodic(string pot, double g, InitialCondition initial1, InitialCondition initial2, double time){
  //vector<double> winit = VectOfDer(initial1);
  vector<double> uinit = transform(initial1);
  vector<double> vinit = transform(initial2);
  vector<Sol> uOuttemp = {};
  Sol ustemp;
  int it = int(time/Timestep + 0.5);
  double ut, wt, vt;
  vector<double> u, w, v;
  vector<double> utemp, vtemp, wtemp;
  for (int i=0;i<it;i++){
   utemp = {};
   vtemp = {};
   wtemp = {};
    if (i==0){  //t==t0
      u = uinit;
      v = vinit;
      //w = winit;
      // initial half step for w_{n+1/2}



      for (int j=0;j<Gridsize;j++){
    wt = 1/Spatialstep * (u.at((j+1)%Gridsize)-u.at((j)%Gridsize));
    wtemp.push_back(wt);
      };
      w = wtemp;

  // initial half step for v^{n+1/2}
      for (int j=0;j<Gridsize;j++){
        Potential Potent(pot,g,v.at(j));
  if (j>0){
    vt =  v.at((j)%Gridsize) + 0.5*(c*c*Timestep/(Spatialstep)) * (w.at((j)%Gridsize) - w.at((j-1)%Gridsize)) - Potent.GetValue() ;
  }
  else {
    vt =   v.at((j)%Gridsize) + 0.5*(c*c*Timestep / Spatialstep) * (w.at((j)%Gridsize) - w.at((j-1)%Gridsize + Gridsize)) - Potent.GetValue() ;
    }
  vtemp.push_back(vt);
      };
     v = vtemp;
      //initial condition for w space halfpstep
}
    // periodic boundary conditions
    else { // t>t0  i>0
      for (int j=0;j<Gridsize;j++){
        //w^{n+1/2} and u^{n+1}

    wt = w.at((j)%Gridsize) + (Timestep/(Spatialstep)) * (v.at((j+1)%Gridsize) - v.at((j)%Gridsize));
    ut = u.at((j)%Gridsize) + Timestep * v.at((j)%Gridsize);

  wtemp.push_back(wt);
  utemp.push_back(ut);
      };
        w = wtemp;
        u = utemp;

    //v^{n+1/2}
for (int j=0;j<Gridsize;j++){
   	Potential Potent(pot,g,v.at(j));
	if (j>0){
    vt = v.at((j)%Gridsize) + (c*c*Timestep/(Spatialstep)) * (w.at((j)%Gridsize) - w.at((j-1)%Gridsize)) - Potent.GetValue();
	}
	else {
  vt =  v.at((j)%Gridsize) +  (c*c*Timestep / Spatialstep) * (w.at((j)%Gridsize) - w.at((j-1)%Gridsize + Gridsize)) - Potent.GetValue();
	}
	vtemp.push_back(vt);
      };
  v = vtemp;

    };
    ustemp.Solvector = u;
    ustemp.time = i*Timestep;
    uOuttemp.push_back(ustemp);

  };
    return uOuttemp; //u;
};

//Leapfrog end





double HyperbolicEquation::CourantFactor(){
  return c * Timestep/Spatialstep;
};



void beginJob(){
  char name[50];
  sprintf(name,"Klein-Gordon.root");
  fileoutput = TFile::Open(name,"RECREATE");
};

void MakeGraph(int N,vector<Sol> vect, double spstep, int timeouti){
  TGraph *gr = new TGraph();
  char name[50];
  sprintf(name,"Solution, f(x), t = %f", vect.at(timeouti).time);
  gr -> SetName(name);
  sprintf(name,"time = %f", vect.at(timeouti).time);
  gr->SetTitle(name);
  if (vect.at(timeouti).Solvector.size() == N){
    for (int i = 0; i<N;i++){
      gr->SetPoint(i,i*spstep,vect.at(timeouti).Solvector.at(i));
    };
  };
  if (gr) {
    gr -> Draw("AC");
    gr -> Write();
   };
};



void MakeGraph(int N,vector<Sol> vect, double spstep, int timeouti, ExactSolution f){
  TGraph *gr = new TGraph();
  TGraph *gr2 = new TGraph();
  TMultiGraph *mg = new TMultiGraph();
  char name[50];
  double temp[2];
  sprintf(name,"Solution, f(x), t = %f", vect.at(timeouti).time);
  mg -> SetName(name);
  sprintf(name,"time = %f", vect.at(timeouti).time);
  mg->SetTitle(name);
  if (vect.at(timeouti).Solvector.size() == N){
    for (int i = 0; i<N;i++){
      temp[0] = i*spstep;
      temp[1] = vect.at(timeouti).time;
      gr->SetPoint(i,i*spstep,vect.at(timeouti).Solvector.at(i));
      gr2->SetPoint(i,i*spstep,f.DoEval(temp));
    };
  };
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerSize(0.5);
  gr->SetMarkerStyle(7);

  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerSize(0.5);
  gr2->SetMarkerStyle(8);

  if ((gr) && (gr2)) {
    mg->Add(gr);
    mg->Add(gr2);
    mg -> Draw("AC");
    mg -> Write();
   };
};

void endJob(){
  fileoutput->Close();
};



void MainKG(string method, string pot, double g, string boundary, double time){
  InitialCondition initial1 = InitialCondition("phi0");
  InitialCondition initial2 = InitialCondition("dphi0");
  beginJob();
  HyperbolicEquation b;
  ExactSolution solution(exactsol);
  b.c = lvc;
  b.Gridsize = GridS;
  b.Spatialstep = Spstep;
  b.Timestep = timest;
  double lim = b.Gridsize*b.Spatialstep;
  vector<Sol> vect = b.Wave(method, pot, g, boundary, initial1, initial2, time);
  if (exactsolution == true){
    for (int i = 0;i<numberofoutput;i++){
      double outtime = int(time*i/(numberofoutput*timest)+0.5);
      MakeGraph(b.Gridsize, vect, b.Spatialstep,outtime,solution);
    }
   }
  else {
    for (int i = 0;i<numberofoutput;i++){
      double outtime = int(time*i/(numberofoutput*timest)+0.5);
      MakeGraph(b.Gridsize, vect, b.Spatialstep,outtime);
  };
  }
  endJob();
};
