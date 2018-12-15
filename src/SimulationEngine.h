//
// Created by lad on 12/15/18.
//

#ifndef SIMULATIONS_SIMULATIONENGINE_H
#define SIMULATIONS_SIMULATIONENGINE_H


class SimulationEngine {
public:
    SimulationEngine();
};

double ***zeroDouble3DArray(int size1, int size2, int size3);
double *getColumnFromTensor(double ***tensor, int i, int size, int k);
double *getGamma(double *px, double *py, int size);
double *numberMinusColumn(double t, double *x, int size);
double *computeE0(double *phi, double chirp, double t, double* x, double shift, double *y, double R0, int size);
double *computeB0(double *phi, double chirp, double t, double* x, double shift, double *y, double R0, int size);
double *computePy(double *py, double *E0, double dt, double *B0, double *px, double *gamma, int size);
double *computePx(double *py, double *E0, double dt, double *B0, double *px, double *gamma, int size);
double *computeY(double *y, double *py, double dt, double *gamma, int size);
double *computeX(double *x, double *px, double dt, double *gamma, int size);
void putColumnIntoTensor(double ***PP, int i, int size, int j, double *column);
double *get5(double *E0, double *B0, double *px, double *gamma, int size);
double *get6(double *BO, double *py, double *gamma, int size);
double *get7(double t, int size);
void printArrayToHDF5(double ****array);


#endif //SIMULATIONS_SIMULATIONENGINE_H
