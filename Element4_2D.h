#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

struct GaussQuadrature{
    int points;         //liczba punktów
    double* nodes;      //węzły
    double* scales;     //wagi
    GaussQuadrature(int points);    //2 or 3 points
};

struct ShapeFunction{
    double N1, N2, N3, N4;

	void setValues(double valueN1, double valueN2, double valueN3, double valueN4);
    void showValues();
};

struct Wall {
    double** N1234; //shape function value table with dim [int point][nodes] ([2][4] or [3][4])
    double* scale;  //waga
    double* ksi;    //coordinates in local system
    double* eta;
};

struct IntegrationPoint{
    double ksi;
    double eta;
    double NValues[4];  //shape function values in integration point
    double scale[2];
};

struct Element4_2D {
    int points;     // 4 or 9
    ShapeFunction* d_ksi;   // dN/dksi
    ShapeFunction* d_eta;   // dN/deta
    IntegrationPoint* intPoints;
    Wall walls[4];

    Element4_2D(GaussQuadrature gq);

    ~Element4_2D();

    void print();
    void showIntPointsInfo();
};