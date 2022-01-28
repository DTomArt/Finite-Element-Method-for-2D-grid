#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include "Element4_2D.h"
#include "Grid.h"

using namespace std;
struct Grid;

struct Jacobian {
    double J[2][2];
    double J_inv[2][2]; //wiersz, kolumna
    double det;
    Jacobian(int el_sk, int pktCalk, Element4_2D* el4_2d, Grid* g);

    void showJ();

    void showJ_inv();
};