#include "Grid.h"

// double** calculate_Hpc(int pktCalk, Jacobian *jac, Element4_2D *el, double K);

Element::Element() {
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            H[i][j] = 0.0;
            HBC[i][j] = 0.0;
            C[i][j] = 0.0;
        }
        P[i] = 0.0;
    }
}

void Element::showH(int elementIndex = -1) {
    if(elementIndex != -1)
        printf("\n\n%d element HBC matrix:\n", elementIndex + 1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++) {
            printf("%f  ", H[i][j]);
        }
        printf("\n");
    }
}

void Element::showHBC(int elementIndex = -1) {
    if(elementIndex != -1)
        printf("\n\n%d element HBC matrix:\n", elementIndex + 1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++) {
            printf("%f  ", HBC[i][j]);
        }
        printf("\n");
    }
}

void Element::showC(int elementIndex = -1) {
    if(elementIndex != -1)
        printf("\n\n%d element C matrix:\n", elementIndex + 1);
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++) {
            printf("%f  ", C[i][j]);
        }
        printf("\n");
    }
}

void Element::showP(int elementIndex = -1) {
    if(elementIndex != -1)
        printf("\n\n%d element P matrix:\n", elementIndex + 1);
    for(int i = 0; i < 4; i++){
            printf("%f  ", P[i]);
        printf("\n");
    }
}

Grid::Grid(double h, double b, int nh, int nb, int integrationPoints, double alpha, 
    double tAmb, double t0, double K, double ro, double c){
    H = h;
    B = b;
    nH = nh;
    nB = nb;
    nN = nH * nB;
    nE = (nH - 1) * (nB - 1);
    gaussQuadrature = new GaussQuadrature(integrationPoints);
    universalElement = new Element4_2D(*gaussQuadrature);

    this->alpha = alpha;
    this->tAmb = tAmb;
    this->t0 = t0;
    this->k = K;
    this->ro = ro;
    this->c = c;
    
    nodes = new Node[nN];
    const double dx = B/(nB-1);
    const double dy = H/(nH-1);
    double x = 0;

    for(int i = 0; i < nN; i++) {
        nodes[i].y = (i % nH) * dy;
        nodes[i].x = x;
        nodes[i].t0 = this->t0;

        //set boundary condition
        if(nodes[i].y == 0. || nodes[i].y == H || nodes[i].x == 0. || nodes[i].x == B)
            nodes[i].BC = 1;

        if((i+1) % nH == 0 && i != 0){
            x += dx;
        }
    }

    elements = new Element[nE];
    int k = 1;

    for(int i = 0; i < nE; i++) {
        elements[i].ID[0] = k;
        elements[i].ID[1] = k + nH;
        elements[i].ID[2] = elements[i].ID[1] + 1;
        elements[i].ID[3] = elements[i].ID[0] + 1;
        k++;
        if(k % nH == 0){ k++; }
    }

    calculate_H_C_in_elements(this->k, this->ro, this->c);
    calculate_Hbc_P_in_elements(this->alpha, this->tAmb);
    aggregation_H_Hbc_P_C_global();
}

void Grid::showNodes(){
    for(int i = 0; i < nN; i++){
        printf("N%d: (%f,%f) BC=%d\n", i+1, nodes[i].x, nodes[i].y, nodes[i].BC);
    }
}

void Grid::showElements(){
    for(int i = 0; i < nE; i++){
        printf("E%d: ",i+1);
        for(int j = 0; j < 4; j++)
            printf("ID%d: %d  ", j+1, elements[i].ID[j]);
        printf("\n");
    }
}

Grid::~Grid(){
    delete [] nodes;
    delete [] elements;
}

void Grid::calculate_H_C_in_elements(double K, double ro, double c){
    for(int i = 0; i < this->nE; i++) {
        for(int j = 0; j < universalElement->points; j++) {
            // cout << "element: " << i + 1 << " pkt: " << j + 1 << endl;
            Jacobian jac(i, j, universalElement, this);
            double** Hmat = calculate_Hpc(j, &jac, universalElement, K);
            double** Cmat = calculate_Cpc(j, &jac, universalElement, c, ro);

            for(int k = 0; k < 4; k++)
                for(int l = 0; l < 4; l++) {
                    this->elements[i].H[k][l] += Hmat[k][l];
                    this->elements[i].C[k][l] += Cmat[k][l];
                }

            for(int k = 0; k < 4; k++){
                delete[] Hmat[k];
                delete[] Cmat[k];
            }
            delete[] Hmat;
            delete[] Cmat;
        }
        // this->elements[i].showH(i);
        // this->elements[i].showC(i);
    }
}

double Grid::distance_nodes(Node node1, Node node2) {
        return sqrt(pow(node1.x - node2.x, 2) + pow(node1.y - node2.y, 2));
}

void Grid::calculate_Hbc_P_in_elements(double alpha, double tAmb) {
    for (int index = 0; index < nE; index++) {
        //for each wall
        for(int w = 0; w < 4; w++) {
            int wN = w;    //wall node index
            int wNN = (w+1) % 4;    //wall next node index
            if (nodes[elements[index].ID[wN] - 1].BC == 1 && nodes[elements[index].ID[wNN] - 1].BC == 1) {
                // // different temperature on left wall
                // if(nodes[elements[index].ID[wN] - 1].x == 0 && nodes[elements[index].ID[wNN] - 1].x == 0)   tAmb = 600;
                // else tAmb = 1200;

                for (int i = 0; i < gaussQuadrature->points; i++) { //for each integration point
                    for (int x = 0; x < 4; x++) {   // for each HBC matrix column and P
                        double detJ = distance_nodes(nodes[elements[index].ID[wN]-1], nodes[elements[index].ID[wNN]-1]) / 2.0;
                        for (int y = 0; y < 4; y++) {   // for each HBC matrix line 
                            elements[index].HBC[x][y] += alpha * gaussQuadrature->scales[i] 
                            * (universalElement->walls[w].N1234[i][x]) 
                            * (universalElement->walls[w].N1234[i][y]) * detJ;
                        }
                        elements[index].P[x] += alpha * gaussQuadrature->scales[i] 
                        * (universalElement->walls[w].N1234[i][x]) * tAmb * detJ;
                    }
                }
            }
        }
        // show HBC or P for each element
        // elements[index].showHBC(index);
        // elements[index].showP(index);
    }
}

void Grid::show_HBC(int elementIndex) {
    printf("\n\n* * * HBC for %d element * * *\n", elementIndex+1);
    for (int x = 0; x < universalElement->points; x++) {
        for (int y = 0; y < universalElement->points; y++) {
            cout << elements[elementIndex].HBC[x][y] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void Grid::aggregation_H_Hbc_P_C_global() {
    // init global matrix H, P, C
    H_global = new double*[nN];
    C = new double* [nN];
    P = new double [nN];
    for(int i = 0; i < nN; i++){
        P[i] = 0.0;
        H_global[i] = new double[nN];
        C[i] = new double[nN];
        for(int j = 0; j < nN; j++){
            H_global[i][j] = 0.0;
            C[i][j] = 0.0;
        }
    }

    //calculate H_global and P vector
    for(int i = 0; i < nE; i++) {   //for every element
        for(int x = 0; x < 4; x++) {    //for every node
            for(int y = 0; y < 4; y++){
                H_global[elements[i].ID[x] - 1][elements[i].ID[y] - 1] += elements[i].H[x][y] + elements[i].HBC[x][y];
                C[elements[i].ID[x] - 1][elements[i].ID[y] - 1] += elements[i].C[x][y];
            }
            P[elements[i].ID[x] - 1] += elements[i].P[x];
        }
    }
}

void Grid::show_H_global() {
    printf("\n\n* * * H_global * * *\n");
    for (int x = 0; x < nN; x++) {
        for (int y = 0; y < nN; y++) {
            cout << this->H_global[x][y] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Grid::show_C() {
    printf("\n\n* * * C * * *\n");
    for (int x = 0; x < nN; x++) {
        for (int y = 0; y < nN; y++) {
            cout << this->C[x][y] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Grid::show_P() {
    printf("\n\n* * * P vector * * *\n");
    for (int x = 0; x < nN; x++)
        cout << this->P[x] << " ";
    cout << endl;
}


double** calculate_Hpc(int pktCalk, Jacobian *jac, Element4_2D *el, double K) {
    double pc_dN_dx[el->points][4];
    double pc_dN_dy[el->points][4];
    pc_dN_dx[pktCalk][0] = jac->J_inv[0][0] * el->d_ksi[pktCalk].N1 + jac->J_inv[0][1] * el->d_eta[pktCalk].N1;
    pc_dN_dx[pktCalk][1] = jac->J_inv[0][0] * el->d_ksi[pktCalk].N2 + jac->J_inv[0][1] * el->d_eta[pktCalk].N2;
    pc_dN_dx[pktCalk][2] = jac->J_inv[0][0] * el->d_ksi[pktCalk].N3 + jac->J_inv[0][1] * el->d_eta[pktCalk].N3;
    pc_dN_dx[pktCalk][3] = jac->J_inv[0][0] * el->d_ksi[pktCalk].N4 + jac->J_inv[0][1] * el->d_eta[pktCalk].N4;

    pc_dN_dy[pktCalk][0] = jac->J_inv[1][0] * el->d_ksi[pktCalk].N1 + jac->J_inv[1][1] * el->d_eta[pktCalk].N1;
    pc_dN_dy[pktCalk][1] = jac->J_inv[1][0] * el->d_ksi[pktCalk].N2 + jac->J_inv[1][1] * el->d_eta[pktCalk].N2;
    pc_dN_dy[pktCalk][2] = jac->J_inv[1][0] * el->d_ksi[pktCalk].N3 + jac->J_inv[1][1] * el->d_eta[pktCalk].N3;
    pc_dN_dy[pktCalk][3] = jac->J_inv[1][0] * el->d_ksi[pktCalk].N4 + jac->J_inv[1][1] * el->d_eta[pktCalk].N4;
    
    // FOR TESTING
    // for(int i = 0; i < 4;i++) {
    //     for(int j = 0; j < 4; j++)
    //         cout << pc_dN_dx[i][j] << " ";
    //     cout << endl;
    // }
    // cout << "\n\n\n";

    // for(int i = 0; i < 4;i++) {
    //     for(int j = 0; j < 4; j++)
    //         cout << pc_dN_dy[i][j] << " ";
    //     cout << endl;
    // }
    // cout << "\n\n\n";
    
    //calculating H formula - calculate dN/dx and dN/dy
    double dn_dx[4][4], dn_dy[4][4];
    for(int j = 0; j < 4; j++){
        for(int k = 0; k < 4; k++){
            dn_dx[j][k] = pc_dN_dx[pktCalk][j] * pc_dN_dx[pktCalk][k];
            dn_dy[j][k] = pc_dN_dy[pktCalk][j] * pc_dN_dy[pktCalk][k];
        }
    }
    // FOR TESTING
    // for(int i = 0; i < 4;i++) {
    //     for(int j = 0; j < 4; j++)
    //         cout << dn_dx[i][j] << " ";
    //     cout << endl;
    // }
    // cout << "\n\n\n";

    // for(int i = 0; i < 4;i++) {
    //     for(int j = 0; j < 4; j++)
    //         cout << dn_dy[i][j] << " ";
    //     cout << endl;
    // }
    // cout << "\n\n\n";

    //add matrixes
    double sum[4][4];
    for(int j = 0; j < 4; j++){
        for(int k = 0; k < 4; k++){
            sum[j][k] = dn_dx[j][k] + dn_dy[j][k];
        }
    }
    
    // init Hpc
    double** Hpc = 0;
    Hpc = new double*[4];
    for(int i = 0;i < 4; i++){
        Hpc[i] = new double[4];
        // for(int j = 0; j < 4; j++){
        //     Hpc[i][j] = 0;
        // }
    }

    //Hpc = K ( thermal conductivity ) * sum * dV( Jacobian det ) * SCALES
    for(int j = 0; j < 4; j++){
        for(int k = 0; k < 4; k++){
            Hpc[j][k] = K * sum[j][k] * jac->det * el->intPoints[pktCalk].scale[0] * el->intPoints[pktCalk].scale[1];
        }
    }
    return Hpc;
}

double** calculate_Cpc(int pktCalk, Jacobian *jac, Element4_2D *el, double c, double ro) {
    //calculate N matrix
    double sum[4][4];
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            sum[i][j] = el->intPoints[pktCalk].NValues[i] * el->intPoints[pktCalk].NValues[j];

    // printf("\n\n");
    // for(int i = 0; i < 4; i++){
    //     for(int j = 0; j < 4; j++)
    //         printf("%lf ", sum[i][j]);
    //     printf("\n");
    // }

    // init Cpc
    double** Cpc;
    Cpc = new double*[4];
    for(int i = 0; i < 4; i++)
        Cpc[i] = new double[4];


    //Cpc = ro * c * N matrix * dV( Jacobian det ) * SCALES
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            Cpc[i][j] = ro * c * sum[i][j] * jac->det * el->intPoints[pktCalk].scale[0] * el->intPoints[pktCalk].scale[1];
    

    // show matrix for specific integration point
    // printf("\n\n* * Matrix C for %d integration point * *\n", pktCalk+1);
    // for(int i = 0; i < 4; i++){
    //     for(int j = 0; j < 4; j++)
    //         printf("%lf ", Cpc[i][j]);
    //     printf("\n");
    // }

    return Cpc;
}