#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Element4_2D.h"
#include "Jacobian.h"

using namespace std;
struct Jacobian;

struct Node {
    double x,y;
    int BC = 0;     // boundary condition
    int t0;
};

struct Element {
    int ID[4];
    double H[4][4]; //H matrix = Hpc1+Hpc2+...
    double HBC[4][4];
    double C[4][4]; //C matrix = Cpc1+Cpc2+...
    double P[4];

    Element();
    void showH(int elementIndex);
    void showHBC(int elementIndex);
    void showC(int elementIndex);
    void showP(int elementIndex);
};

struct Grid {
    double H;   //height
    double B;   //width
    int nH;     //nodes in height
    int nB;     //nodes in width
    int nN;
    int nE;

    double alpha;  // [W/m^2*K]
    double tAmb;   // ambient temperature [C]
    double t0; //initial temperature
    double k;  // conductivity [W/(m*C(temp))]
    double ro; // density [kg/m^3]
    double c;  // specific heat [J/(kg*C(temp))] (cieplo wlasciwe)

    Node* nodes;
    Element* elements;
    Element4_2D* universalElement;
    GaussQuadrature* gaussQuadrature;
    double** H_global;
    double ** C;
    double *P;
    Grid(double h, double b, int nh, int nb, int integrationPoints, double alpha, double tAmb, double t0, double K, double ro, double c);
    void showNodes();
    void showElements();
    ~Grid();
    void calculate_H_C_in_elements(double K, double ro, double c);
    double distance_nodes(Node node1, Node node2);
    void calculate_Hbc_P_in_elements(double alpha, double tAmb);
    void show_HBC(int elementIndex);
    void aggregation_H_Hbc_P_C_global();
    void show_H_global();
    void show_P();
    void show_C();
};

double** calculate_Hpc(int pktCalk, Jacobian *jac, Element4_2D *el, double K);
double** calculate_Cpc(int pktCalk, Jacobian *jac, Element4_2D *el, double c, double ro);