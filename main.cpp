#include <iostream>
#include "Grid.h"
#include "GaussEliminationMethod.h"

int main() {
    //grid initial parameters for edit
        double B = 0.1, H = 0.1;    //height and width
        int nB = 31, nH = 31;   //number of nodes in width and height
        int points = 3;   // Gauss integration scheme number of points (in this program can be 2 or 3)

        double alpha = 300.;  // [W/m^2*K]
        double tAmb = 1200.;   // ambient temperature [C]
        double t0 = 100.; //initial temperature
        double K = 25.;  // conductivity [W/(m*C(temp))]
        double ro = 7800.; // density [kg/m^3]
        double c = 700.;  // specific heat [J/(kg*C(temp))] (cieplo wlasciwe)
        
        int tau = 100;    // simulation time
        int dTau = 1;     //simulation step

    // initializing grid
    Grid grid(H, B, nH, nB, points, alpha, tAmb, t0, K, ro, c);
    
    int Nn = grid.nN;   //number of nodes in grid
    double** AB = new double* [Nn];
    double* t1 = new double[Nn];
    for(int i = 0; i < Nn; i++) AB[i] = new double[Nn+1];

    double t0Vec[Nn];
    fill_n(t0Vec, Nn, t0);    //fill table with tInit

    // simulation
    printf("Time[s]\tMinTemp[s]\tMaxTemp[s]\n");
    for(int time = dTau; time <= tau; time+=dTau) {
        double Ct0[Nn] = {0.};
        grid.calculate_H_C_in_elements(K,ro,c); //simulation is performed with constant initial parameters
        grid.calculate_Hbc_P_in_elements(alpha,tAmb);

        // right side of equation ([C]/dTau*{t0} + P)
        for(int i = 0; i < Nn; i++)
            for(int j = 0; j < Nn; j++)
                Ct0[i] += grid.C[i][j]/dTau * t0Vec[j];


        for(int i = 0; i < Nn; i++)
            for(int j = 0; j <= Nn; j++)
                if( j == Nn)
                    AB[i][j] = Ct0[i] + grid.P[i];
                else AB[i][j] = grid.H_global[i][j] + grid.C[i][j]/dTau;
        
        // printf("{P} = {P}+{[C]/dT}*{T0}\n");
        // for(int i = 0; i < Nn; i++){
        //     cout << AB[i][Nn] << " ";
        // }
        // cout << endl << endl;

        // printf(" MATRIX [H] = [H] + [C]/dT and {P} = {P}+{[C]/dT}*{T0}\n");
        // for(int i = 0; i < Nn; i++){
        //     for(int j = 0; j <= Nn; j++)
        //         cout << AB[i][j] << " ";
        //     cout << endl;
        // }

        if (!gaussEliminationMethod(Nn, AB, t1))
        {
            cout << "DIVIDER EQUALS 0, ERROR GAUSS\n";
            return -1;
        }

        // calculate min, max
        double min = t1[0], max = t1[0];
        for(int i = 1; i < Nn; i++){
            if(t1[i] < min)  min = t1[i];
            else if (t1[i] > max) max = t1[i];
        }

        //final result
        printf("%d\t%.3f   \t%.3f\n", time, min, max);

        // // show output vector
        // printf("\n");
        // for(int i = 0; i< Nn; i++) printf("x%d = %lf\n", i+1, t[i]);
        // printf("\n");

        for(int i = 0; i < Nn; i++) t0Vec[i] = t1[i];
    }
    
    for(int i = 0; i < Nn ; i++) delete [] AB[i];
    delete [] AB;
    delete [] t1;
    return 0;
}

//g++ .cpp