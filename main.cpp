#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

using namespace std;

double gamma(double T);

double gamma_0(double T);

double heatCapacity(int i, double T);

double heatCapacity_0(int i, double T);

double func1(double p, double nu, double T);

double func2(double p, double nu, double T);

double func3(double p, double nu, double T);

double T_sum();

void copyv(double *a, const int n, double *b);

void l1_func(double p, double nu, double T, double *l1);

void l2_func(double p, double nu, double T, double *l2);

void l3_func(double p, double nu, double T, double *l3);

double heatCapacity_1(int i, double T);

double derivate_gamma(double T);

double getInversedJacobian(const double *l1, const double *l2, const double *l3);

void setParams(const double *l1, const double *l2, const double *l3, double *L, double det);

const double R = 8.31;
const double w_bef[3] = {0.055, 0.220, 0.725};
const double w_aft[3] = {0.151, 0.124, 0.725};
const double mu_bef[3] = {0.016, 0.032, 0.028};
const double mu_aft[3] = {0.044, 0.018, 0.028};

const int p_0 = (const int) 1e5;
const int T_0 = 298;
const int Q = (const int) 2.76e6;
const double ro_0 = 1.1;
const double nu_0 = 1 / ro_0;
const double eps = 1e-5;

int main() {
    //initializing in 3 1 go!
    double p, nu, T, v, D;
    double f1, f2, f3;
    double gm, ro;
    double d = eps + 1;
    int count = 0;
    double *P1 = (double *) malloc(3 * sizeof(double));
    double *P2 = (double *) malloc(3 * sizeof(double));
    double *buf = (double *) malloc(3 * sizeof(double));
    double *l1 = (double *) malloc(3 * sizeof(double));
    double *l2 = (double *) malloc(3 * sizeof(double));
    double *l3 = (double *) malloc(3 * sizeof(double));
    double *L = (double *) malloc(9 * sizeof(double));
    //initial approximation//
    gm = gamma(T_0);
    ro = ((gm + 1) / gm) * ro_0;
    p = ro_0 * 2 * (gm - 1) * Q;
    nu = 1 / ro;
    T = p / (ro * R * T_sum());
    //the vector of initial conditions//
    P1[0] = 50000;//p;
    P1[1] = 20;//nu;
    P1[2] = 5000;//T;

    cout << "Initial approximation" << endl;
    printf("p = %f   nu = %f   T = %f\n", P1[0], P1[1], P1[2]);
    printf("***********************************************************\n");
    //Jacobian//
    l1_func(p, nu, T, l1);
    l2_func(p, nu, T, l2);
    l3_func(p, nu, T, l3);
    //Inversed Jacobian//
    double det = getInversedJacobian(l1, l2, l3);
    setParams(l1, l2, l3, L, det);
    while (d > eps) {
        f1 = func1(P1[0], P1[1], P1[2]);
        f2 = func2(P1[0], P1[1], P1[2]);
        f3 = func3(P1[0], P1[1], P1[2]);
        P2[0] = P1[0] - (L[0] * f1 + L[1] * f2 + L[2] * f3);
        //f2 = func2(P2[0], P1[1], P1[2]);
        P2[1] = P1[1] - (L[3] * f1 + L[4] * f2 + L[5] * f3);
        //f3 = func3(P2[0], P2[1], P1[2]);
        P2[2] = P1[2] - (L[6] * f1 + L[7] * f2 + L[8] * f3);
        d = sqrt(pow(f1, 2) + pow(f2, 2) + pow(f3, 2));  //discrepancy//
        count++;
        printf("iteration: %d  error %f\n", count, d);
        printf("f1 = %f  f2 = %f  f3 = %f\n", f1, f2, f3);
        printf("p = %f   nu = %f   T = %f\n", P2[0], P2[1], P2[2]);
        printf("-----------------------------------------------------------\n");
        //Jacobian//
        l1_func(P2[0], P2[1], P2[2], l1);
        l2_func(P2[0], P2[1], P2[2], l2);
        l3_func(P2[0], P2[1], P2[2], l3);
        //Inversed Jacobian//
        det = getInversedJacobian(l1, l2, l3);
        setParams(l1, l2, l3, L, det);
        //getchar();

        if (d <= eps) {
            p = P2[0];
            nu = P2[1];
            T = P2[2];
        }

        copyv(P1, 3, buf);
        copyv(P2, 3, P1);
        copyv(buf, 3, P2);
    }

    D = nu_0 * sqrt((p - p_0) / (nu_0 - nu));
    v = (nu_0 - nu) * sqrt((p - p_0) / (nu_0 - nu));
    gm = gamma(T);

    cout << "***********************************************************" << endl;
    printf("p = %f  ro = %f  T = %f\n", p, 1 / nu, T);
    printf("v = %f  D = %f  gamma = %f\n", v, D, gm);
    cout << count << endl;

//todo: use smart pointers, hotya rabotaet i ladno
    free(P1);
    free(P2);
    free(buf);
    getchar();
    return 0;
}

void setParams(const double *l1, const double *l2, const double *l3, double *L, double det) {
    L[0] = (l2[1] * l3[2] - l3[1] * l2[2]) / det;
    L[1] = -(l1[1] * l3[2] - l3[1] * l1[2]) / det;
    L[2] = (l1[1] * l2[2] - l2[1] * l1[2]) / det;
    L[3] = -(l2[0] * l3[2] - l3[0] * l2[2]) / det;
    L[4] = (l1[0] * l3[2] - l3[0] * l1[2]) / det;
    L[5] = -(l1[0] * l2[2] - l2[0] * l1[2]) / det;
    L[6] = (l2[0] * l3[1] - l3[0] * l2[1]) / det;
    L[7] = -(l1[0] * l3[1] - l3[0] * l1[1]) / det;
    L[8] = (l1[0] * l2[1] - l2[0] * l1[1]) / det;
}

double getInversedJacobian(const double *l1, const double *l2, const double *l3) {
    double det = l1[0] * (l2[1] * l3[2] - l3[1] * l2[2]) - l1[1] * (l2[0] * l3[2] - l3[0] * l2[2]) +
                 l1[2] * (l2[0] * l3[1] - l3[0] * l2[1]);
    return det;
}

//TODO: Replace with lambda
double gamma(double T) {
    double gm, up = 0, down = 0;
    int i;
    for (i = 0; i < 3; i++) {
        up += w_aft[i] / mu_aft[i];
        down += (w_aft[i] * heatCapacity(i, T)) / mu_aft[i];
    }

    gm = 1 + R * (up / down);
    return gm;
}

//TODO: Replace with lambda
double gamma_0(double T) {
    double gm, up = 0, down = 0;
    int i;
    for (i = 0; i < 3; i++) {
        up += w_bef[i] / mu_bef[i];
        down += (w_bef[i] * heatCapacity_0(i, T)) / mu_bef[i];
    }
    gm = 1 + R * (up / down);
    return gm;
}

//TODO:REFACTOR!!!! MAKE LOADING FROM XML,CSV,TXT,ETC
double heatCapacity(int i, double T) {
    double c_v = 0;
    if (i == 0) {
        if (T <= 1000) {
            c_v = R * (2.35677352 + 8.98459677e-3 * T - 7.12356269e-6 * pow(T, 2) + 2.45919022e-9 * pow(T, 3) -
                       1.43699548e-13 * pow(T, 4) - 1);
        }
        else
            c_v = R * (3.85746029 + 4.41437026e-3 * T - 2.21481404e-6 * pow(T, 2) + 5.23490188e-10 * pow(T, 3) -
                       4.72084164e-14 * pow(T, 4) - 1);
    }
    if (i == 1) {
        if (T <= 1000) {
            c_v = R * (4.19864056 - 2.03643410e-3 * T + 6.52040211e-6 * pow(T, 2) - 5.48797062e-9 * pow(T, 3) +
                       1.77197817e-12 * pow(T, 4) - 1);
        }
        else
            c_v = R * (3.03399249 + 2.17691804e-3 * T - 1.64072518e-7 * pow(T, 2) - 9.70419870e-11 * pow(T, 3) +
                       1.68200992e-14 * pow(T, 4) - 1);
    }
    if (i == 2) {
        if (T <= 1000) {
            c_v = R * (0.03298677e2 + 0.14082404e-2 * T - 0.03963222e-4 * pow(T, 2) + 0.05641515e-7 * pow(T, 3) -
                       0.02444854e-10 * pow(T, 4) - 1);
        }
        else
            c_v = R * (0.02926640e2 + 0.14879768e-2 * T - 0.05684760e-5 * pow(T, 2) + 0.10097038e-9 * pow(T, 3) -
                       0.06753351e-13 * pow(T, 4) - 1);
    }
    return c_v;
}
//TODO:REFACTOR!!!! MAKE LOADING FROM XML,CSV,TXT,ETC

double heatCapacity_0(int i, double T) {

    double c_v = 0;
    if (i == 0) {
        if (T <= 1000) {
            c_v = R * (5.14987613 - 1.36709788e-2 * T + 4.91800599e-5 * pow(T, 2) - 4.84743026e-8 * pow(T, 3) +
                       1.66693956e-11 * pow(T, 4) - 1);
        }
        else
            c_v = R * (7.48514950e-2 + 1.33909467e-2 * T - 5.73285809e-6 * pow(T, 2) + 1.22292535e-9 * pow(T, 3) -
                       1.01815230e-13 * pow(T, 4) - 1);
    }
    if (i == 1) {
        if (T <= 1000) {
            c_v = R * (3.78245636 - 2.99673416e-3 * T + 9.84730201e-6 * pow(T, 2) - 9.68129509e-9 * pow(T, 3) +
                       3.24372837e-12 * pow(T, 4) - 1);
        }
        else
            c_v = R * (3.28253784 + 1.48308754e-3 * T - 7.57966669e-7 * pow(T, 2) + 2.09470555e-10 * pow(T, 3) -
                       2.16717794e-14 * pow(T, 4) - 1);
    }
    if (i == 2) {
        if (T <= 1000) {
            c_v = R * (0.03298677e2 + 0.14082404e-2 * T - 0.03963222e-4 * pow(T, 2) + 0.05641515e-7 * pow(T, 3) -
                       0.02444854e-10 * pow(T, 4) - 1);
        }
        else
            c_v = R * (0.02926640e2 + 0.14879768e-2 * T - 0.05684760e-5 * pow(T, 2) + 0.10097038e-9 * pow(T, 3) -
                       0.06753351e-13 * pow(T, 4) - 1);
    }
    return c_v;
}
//TODO:REFACTOR!!!! MAKE LOADING FROM XML,CSV,TXT,ETC

double heatCapacity_1(int i, double T) {

    double c_v = 0;
    if (i == 1) {
        if (T <= 1000) {
            c_v = R * (-2.99673416e-3 + 2 * 9.84730201e-6 * pow(T, 1) - 3 * 9.68129509e-9 * pow(T, 2) +
                       4 * 3.24372837e-12 * pow(T, 3));
        }
        else
            c_v = R * (1.48308754e-3 - 2 * 7.57966669e-7 * pow(T, 1) + 3 * 2.09470555e-10 * pow(T, 2) -
                       4 * 2.16717794e-14 * pow(T, 3));
    }
    if (i == 2) {
        if (T <= 1000) {
            c_v = R * (0.14082404e-2 - 2 * 0.03963222e-4 * pow(T, 1) + 3 * 0.05641515e-7 * pow(T, 2) -
                       4 * 0.02444854e-10 * pow(T, 30));
        }
        else
            c_v = R * (0.14879768e-2 - 2 * 0.05684760e-5 * pow(T, 1) + 3 * 0.10097038e-9 * pow(T, 2) -
                       4 * 0.06753351e-13 * pow(T, 3));
    }
    return c_v;
}

double func1(double p, double nu, double T) {
    double res;
    res = 0.5 * (p_0 + p) * (nu_0 - nu) + Q - ((p * nu) / (gamma(T) - 1)) + ((p_0 * nu_0) / (gamma_0(T_0) - 1));
    return res;
}

double func2(double p, double nu, double T) {
    double res;
    res = ((gamma(T) * p) / nu) - ((p - p_0) / (nu_0 - nu));
    return res;
}

double func3(double p, double nu, double T) {
    double res;
    res = (R * T * T_sum()) / nu - p;
    return res;
}

double T_sum() {
    double res = 0;
    for (int i = 0; i < 3; i++) {
        res += w_aft[i] / mu_aft[i];
    }
    return res;
}

void copyv(double *a, const int n, double *b) {
    memcpy(b, a, sizeof(double) * n);
}

double derivate_gamma(double T) {
    double gm, up = 0, down = 0;
    int i;
    for (i = 0; i < 3; i++) {
        up += w_aft[i] / mu_aft[i];
        down += -heatCapacity_1(i, T) / ((w_aft[i] * pow(heatCapacity(i, T), 2)) / mu_aft[i]);
    }

    gm = R * (up * down);
    return gm;
}

void l1_func(double p, double nu, double T, double *l1) {
    l1[0] = 0.5 * (1 / 1.1 - nu) - nu / (gamma(T) - 1);
    l1[1] = -0.5 * (100000 + p) - p / (gamma(T) - 1);
    l1[2] = p * nu * derivate_gamma(T) / pow(gamma(T) - 1, 2);
}

void l2_func(double p, double nu, double T, double *l2) {
    l2[0] = gamma(T) / nu - 1 / (1 / 1.1 - nu);
    l2[1] = -gamma(T) * p / pow(nu, 2) - (p - 100000) / pow(1 / 1.1 - nu, 2);
    l2[2] = derivate_gamma(T) * p / nu;
}

//p- nahuj ne nughna, no pust budet
void l3_func(double p, double nu, double T, double *l3) {
    l3[0] = -1;
    l3[1] = -R * T * T_sum() / pow(nu, 2);
    l3[2] = R * T_sum() / nu;
}
