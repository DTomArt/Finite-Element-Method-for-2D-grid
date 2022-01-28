#include "Element4_2D.h"

using namespace std;

GaussQuadrature::GaussQuadrature(int points){    //2 or 3 points
    this->points = points;
    this->nodes = new double[points];
    this->scales = new double[points];
    if(points == 2){
        nodes[0] = -1/sqrt(3);
        nodes[1] = 1/sqrt(3);
        scales[0] = 1;
        scales[1] = 1;
    } else if(points == 3){
        nodes[0] = -sqrt(3.0 / 5.0);
        nodes[1] = 0;
        nodes[2] = sqrt(3.0 / 5.0);
        scales[0] = 5.0 / 9.0;
        scales[1] = 8.0 / 9.0;
        scales[2] = 5.0 / 9.0;
    }
}

void ShapeFunction::setValues(double valueN1, double valueN2, double valueN3, double valueN4) {
    this->N1 = valueN1;
    this->N2 = valueN2;
    this->N3 = valueN3;
    this->N4 = valueN4;
}

void ShapeFunction::showValues(){
    cout << setprecision(5) << N1 << " "; 
    cout << setprecision(5) << N2 << " "; 
    cout << setprecision(5) << N3 << " "; 
    cout << setprecision(5) << N4 << " ";
    cout << endl;
}

Element4_2D::Element4_2D(GaussQuadrature gq){
    this->points = gq.points*gq.points;
    double n1,n2,n3,n4;
    d_ksi = new ShapeFunction[points];
    d_eta = new ShapeFunction[points];
    intPoints = new IntegrationPoint[points];

    int temp1 = 0;
    int temp2 = 0;
    for(int i = 0; i < points; i++) {

        if (i % gq.points == 0 && i != 0)
        {
            temp1++;
            temp2 = 0;
        }

        n1 = -(1.0 / 4.0) * (1.0 - gq.nodes[temp1]);
        n2 = (1.0 / 4.0) * (1.0 - gq.nodes[temp1]);
        n3 = (1.0 / 4.0) * (1.0 + gq.nodes[temp1]);
        n4 = -(1.0 / 4.0) * (1.0 + gq.nodes[temp1]);
        d_ksi[i].setValues(n1, n2, n3, n4);

        n1 = -(1.0 / 4.0) * (1.0 - gq.nodes[temp2]);
        n2 = -(1.0 / 4.0) * (1.0 + gq.nodes[temp2]);
        n3 = (1.0 / 4.0) * (1.0 + gq.nodes[temp2]);
        n4 = (1.0 / 4.0) * (1.0 - gq.nodes[temp2]);
        d_eta[i].setValues(n1, n2, n3, n4);

        intPoints[i].ksi = gq.nodes[temp2];
        intPoints[i].eta = gq.nodes[temp1];

        intPoints[i].NValues[0] = 0.25 * (1.0 - intPoints[i].ksi) * (1.0 - intPoints[i].eta);
        intPoints[i].NValues[1] = 0.25 * (1.0 + intPoints[i].ksi) * (1.0 - intPoints[i].eta);
        intPoints[i].NValues[2] = 0.25 * (1.0 + intPoints[i].ksi) * (1.0 + intPoints[i].eta);
        intPoints[i].NValues[3] = 0.25 * (1.0 - intPoints[i].ksi) * (1.0 + intPoints[i].eta);

        intPoints[i].scale[0] = gq.scales[temp1];
        intPoints[i].scale[1] = gq.scales[temp2];
        // cout << e << " " << k << endl;
        // cout << intPoints[i].scale[0] << " " << intPoints[i].scale[1] << endl;
        
        temp2++;
    }

    //walls: gq.points - quantity of integration points on walls: 2 or 3
    for(int i = 0; i < 4; i++){
        walls[i].N1234 = new double*[gq.points];
        walls[i].ksi = new double[gq.points];
        walls[i].eta = new double[gq.points];
        walls[i].scale = new double[gq.points];
        for(int j=0;j<gq.points;j++){
            walls[i].N1234[j] = new double[4];
            // for(int k=0;k<4;k++)
            //     walls[i].N1234[j][k] = 0;
        }
    }

    //wall 1
    for(int i = 0; i < gq.points; i++){
        walls[0].ksi[i] = gq.nodes[i];
        walls[0].eta[i] = -1;
        walls[0].scale[i] = gq.scales[i];
    }

    //wall 2
    for(int i = 0; i < gq.points; i++){
        walls[1].ksi[i] = 1;
        walls[1].eta[i] = gq.nodes[i];
        walls[1].scale[i] = gq.scales[i];
    }

    //wall 3
    for(int i = 0; i < gq.points; i++){
        walls[2].ksi[i] = gq.nodes[i];
        walls[2].eta[i] = 1;
        walls[2].scale[i] = gq.scales[i];
    }

    //wall 4
    for(int i = 0; i < gq.points; i++){
        walls[3].ksi[i] = -1;
        walls[3].eta[i] = gq.nodes[i];
        walls[3].scale[i] = gq.scales[i];
    }

    //N (shape functions) on walls in integration points
    for (int j = 0; j < 4; j++) {
        for (int i = 0; i < gq.points; i++) {
            walls[j].N1234[i][0] = 0.25 * (1.0 - walls[j].ksi[i]) * (1.0 - walls[j].eta[i]);
            walls[j].N1234[i][1] = 0.25 * (1.0 + walls[j].ksi[i]) * (1.0 - walls[j].eta[i]);
            walls[j].N1234[i][2] = 0.25 * (1.0 + walls[j].ksi[i]) * (1.0 + walls[j].eta[i]);
            walls[j].N1234[i][3] = 0.25 * (1.0 - walls[j].ksi[i]) * (1.0 + walls[j].eta[i]);
        }
    }

}

Element4_2D::~Element4_2D(){
    delete [] d_ksi;
    delete [] d_eta;
    delete [] intPoints;
}

void Element4_2D::print() {
    printf("\n\n * * Universal element 2D * *\n");
    cout << "d_N/d_ksi:" << endl;
    for(int i = 0; i < points; i++){
        d_ksi[i].showValues();
    }
    cout << endl;

    cout << "d_N/d_eta:" << endl;
    for(int i = 0; i < points; i++){
        d_eta[i].showValues();
    }
    cout << endl;
}

void Element4_2D::showIntPointsInfo() {
    printf("\n\n * * Universal element 2D * *\n");
    printf(" Integration points:\n");
    printf(" Pc\tksi\t\teta\t\t N1\tN2\tN3\tN4\n");
    for(int i = 0; i < points; i++){
        printf(" %d\t%lf\t%lf\t", i+1, intPoints[i].ksi, intPoints[i].eta);
        for(int j = 0; j < 4; j++ )
            printf("%.5g ",intPoints[i].NValues[j]);
        printf("\n");
    }
}