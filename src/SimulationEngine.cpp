//
// Created by lad on 12/15/18.
//

#include <cstdio>
#include <cmath>
#include "SimulationEngine.h"

/**
 * Размерность пространства
 */
const int DIMENSIONS = 9;

/**
 * Число шагов
 */
const int N = 20;

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

SimulationEngine::SimulationEngine() {
    printf("\nНачало симуляции\n");

    int nmax = 100;
    int kmax = nmax * 100;

    double ***PP = zeroDouble3DArray(DIMENSIONS, kmax, N+1);
    //double ***EY = zeroDouble3DArray(DIMENSIONS, kmax, N+1);
    //double ***HZ = zeroDouble3DArray(DIMENSIONS, kmax, N+1);

    for (int i = 0; i < nmax; i++) {
        for (int tempI = i*nmax; tempI < i*nmax + nmax; tempI++) {
            PP[0][tempI][0] = 10.0 + (-1.0 + 2.0 / nmax * (tempI - i*nmax));
        }
        for (int tempI = i*nmax; tempI < i*nmax + nmax; tempI++) {
            PP[1][tempI][0] = 10.0;
        }
        for (int tempI = i*nmax; tempI < i*nmax + nmax; tempI++) {
            PP[2][tempI][0] = 10.0 * i / nmax;
        }
        for (int tempI = i*nmax; tempI < i*nmax + nmax; tempI++) {
            PP[3][tempI][0] = 0.0 + -1.0 + 2.0 / nmax * (tempI - i*nmax);
        }
    }

    double *x = getColumnFromTensor(PP, 0, kmax, 0);
    double *y = getColumnFromTensor(PP, 1, kmax, 0);
    double *px = getColumnFromTensor(PP, 2, kmax, 0);
    double *py = getColumnFromTensor(PP, 3, kmax, 0);

    double *gamma = getGamma(px, py, kmax);

    int ii = 10;
    double R0 = 2.5;
    double dt=0.01;
    double shift = +10.0;

    for (int i = 0; i < 10*N; i++) {
        double t = ((double) i)/100.0;
        double chirp = -0.0;
        double *phi = numberMinusColumn(t, x, kmax);
        double *E0 = computeE0(phi, chirp, t, x, shift, y, R0, kmax);
        double *B0 = computeB0(phi, chirp, t, x , shift, y, R0, kmax);
        py = computePy(py, E0, dt, B0, px, gamma, kmax);
        px = computePx(py, E0, dt, B0, px, gamma, kmax);
        y = computeY(y, py, dt, gamma, kmax);
        x = computeX(x, px, dt, gamma, kmax);
        gamma = getGamma(px, py, kmax);
        if (ii == 10) {
            ii=0;
            //print?
            int iii = (int) i/10;
            putColumnIntoTensor(PP, 0, kmax, iii + 1, x); //TODO check this out
            putColumnIntoTensor(PP, 1, kmax, iii + 1, y);
            putColumnIntoTensor(PP, 2, kmax, iii + 1, px);
            putColumnIntoTensor(PP, 3, kmax, iii + 1, py);
            //EY
            //HZ
            putColumnIntoTensor(PP, 4, kmax, iii + 1, get5(E0, B0, px, gamma, kmax));
            putColumnIntoTensor(PP, 5, kmax, iii + 1, get6(B0, py, gamma, kmax));
            putColumnIntoTensor(PP, 6, kmax, iii + 1, get7(t, kmax));

        }
        ii++;
    }

    printf("\nКонец симуляции\n");
}

double ***zeroDouble3DArray(int size1, int size2, int size3) {

    double ***result;
    result = new double**[size1];
    /*try {*/
    for (int index1 = 0; index1 < size1; index1++) {
        result[index1] = new double *[size2];
        for (int index2 = 0; index2 < size2; index2++) {
            result[index1][index2] = new double[size3];
            for (int index3 = 0; index3 < size3; index3++) {
                result[index1][index2][index3] = 0.0;
            }
        }
    }
    /*  } catch (const std::exception exception) {
          printf( exception.what() );
          result = nullptr;
      }*/

    return result;
}

double *getColumnFromTensor(double ***tensor, int i, int size, int k) {

    auto *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index]=tensor[i][index][k];
    }

    return result;
}

double *getGamma(double *px, double *py, int size) {

    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = sqrt(1.0 + px[index]*px[index] + py[index]*py[index]);
    }

    return result;
}

double *numberMinusColumn(double t, double *x, int size) {

    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = t - x[index];
    }

    return result;
}

double *computeE0(double *phi, double chirp, double t, double* x, double shift, double *y, double R0, int size) {

    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = 2.0*
                        cos(phi[index]*6.28*(1+phi[index]*chirp))
                        *exp(-pow(t-x[index]+shift, 2.0)/50.0)
                        *exp(-pow(y[index]-10.0, 2)/5.0)
                        +
                        0.05*(y[index]-10.0)/R0
                        *exp(-pow(y[index]-10.0, 2)/(R0*R0));
    }

    return result;
}

double *computeB0(double *phi, double chirp, double t, double* x, double shift, double *y, double R0, int size) {

    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = 2.0*
                        cos(phi[index]*6.28*(1+phi[index]*chirp))
                        *exp(-pow(t-x[index]+shift, 2.0)/50.0)
                        *exp(-pow(y[index]-10.0, 2)/5.0)
                        +
                        0.05*(y[index]-10.0)/R0
                        *exp(-pow(y[index]-10.0, 2)/(R0*R0));
    }

    return result;
}

double *computePy(double *py, double *E0, double dt, double *B0, double *px, double *gamma, int size) {
    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = py[index] - 6.28 * E0[index] * dt + 6.28 * B0[index] * dt * px[index] / gamma[index];
    }

    return result;
}


double *computePx(double *py, double *E0, double dt, double *B0, double *px, double *gamma, int size) {
    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = px[index] - 6.28*B0[index]*dt*py[index]/gamma[index];
    }

    return result;
}

double *computeY(double *y, double *py, double dt, double *gamma, int size) {
    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = y[index] + py[index] * dt / gamma[index];
    }

    return result;
}

double *computeX(double *x, double *px, double dt, double *gamma, int size) {
    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = x[index] + px[index] * dt / gamma[index];
    }

    return result;
}

void putColumnIntoTensor(double ***PP, int i, int size, int k, double *column) {

    for (int index = 0; index < size; index++) {
        PP[i][index][k] = column[index];
    }
}

double *get5(double *E0, double *B0, double *px, double *gamma, int size) {
    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = -0*6.28*E0[index]+6.28*B0[index]*px[index]/gamma[index];
    }

    return result;
}

double *get6(double *B0, double *py, double *gamma, int size) {
    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = -6.28*B0[index]*py[index]/gamma[index];
    }

    return result;
}

double *get7(double t, int size) {
    double *result = new double[size];

    for (int index = 0; index < size; index++) {
        result[index] = t;
    }

    return result;
}