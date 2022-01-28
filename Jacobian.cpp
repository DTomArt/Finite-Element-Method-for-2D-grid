#include "Jacobian.h"

using namespace std;

Jacobian::Jacobian(int el_sk, int pktCalk, Element4_2D* el4_2d, Grid* g){ 
    // J[0][0]
    const double dx_dKsi =
    el4_2d->d_ksi[pktCalk].N1 * g->nodes[g->elements[el_sk].ID[0] - 1].x +
    el4_2d->d_ksi[pktCalk].N2 * g->nodes[g->elements[el_sk].ID[1] - 1].x +
    el4_2d->d_ksi[pktCalk].N3 * g->nodes[g->elements[el_sk].ID[2] - 1].x +
    el4_2d->d_ksi[pktCalk].N4 * g->nodes[g->elements[el_sk].ID[3] - 1].x;
    
    // J[0][1]
    const double dy_dKsi =
    el4_2d->d_ksi[pktCalk].N1 * g->nodes[g->elements[el_sk].ID[0] - 1].y +
    el4_2d->d_ksi[pktCalk].N2 * g->nodes[g->elements[el_sk].ID[1] - 1].y +
    el4_2d->d_ksi[pktCalk].N3 * g->nodes[g->elements[el_sk].ID[2] - 1].y +
    el4_2d->d_ksi[pktCalk].N4 * g->nodes[g->elements[el_sk].ID[3] - 1].y;

    // J[1][0]
    const double dx_dEta =
    el4_2d->d_eta[pktCalk].N1 * g->nodes[g->elements[el_sk].ID[0] - 1].x +
    el4_2d->d_eta[pktCalk].N2 * g->nodes[g->elements[el_sk].ID[1] - 1].x +
    el4_2d->d_eta[pktCalk].N3 * g->nodes[g->elements[el_sk].ID[2] - 1].x +
    el4_2d->d_eta[pktCalk].N4 * g->nodes[g->elements[el_sk].ID[3] - 1].x;

    // J[1][1]
    const double dy_dEta =
    el4_2d->d_eta[pktCalk].N1 * g->nodes[g->elements[el_sk].ID[0] - 1].y +
    el4_2d->d_eta[pktCalk].N2 * g->nodes[g->elements[el_sk].ID[1] - 1].y +
    el4_2d->d_eta[pktCalk].N3 * g->nodes[g->elements[el_sk].ID[2] - 1].y +
    el4_2d->d_eta[pktCalk].N4 * g->nodes[g->elements[el_sk].ID[3] - 1].y;

    J[0][0] = dx_dKsi;
    J[0][1] = dy_dKsi;
    J[1][0] = dx_dEta;
    J[1][1] = dy_dEta;

    //test
    // J[0][0] = 0.0125;
    // J[0][1] = 0;
    // J[1][0] = 0;
    // J[1][1] = 0.0125;
    //
    
    det = (J[0][0] * J[1][1]) - (J[0][1] * J[1][0]);

    J_inv[0][0] = J[1][1];
    J_inv[0][1] = - J[0][1];
    J_inv[1][0] = - J[1][0];
    J_inv[1][1] = J[0][0];

    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++){
            J_inv[i][j] /= det;
        }
    }
}

void Jacobian::showJ(){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            cout << J[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Jacobian::showJ_inv(){
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            cout << J_inv[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}