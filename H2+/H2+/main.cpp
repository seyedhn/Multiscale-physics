using namespace std;

#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;
//using Eigen::MatrixXd;
//using Eigen::EigenSolver;


const double PI  =3.141592653589793238463;

//For Hydrogen
double alpha(int i){
    double x;
    
    switch(i){
        case 0 : x = 13.00773;
            break;
        case 1 : x = 1.962079;
            break;
        case 2 : x = 0.444529;
            break;
        case 3 : x = 0.1219492;
            break;
        case 4 : x = 13.00773;
            break;
        case 5 : x = 1.962079;
            break;
        case 6 : x = 0.444529;
            break;
        case 7 : x = 0.1219492;
            break;
        default : x = 100;
    }
    return x;
}

int R(int i){
    int x;
    
    switch(i){
        case 0 ... 3 : x = 4;
            break;
        case 4 ... 7 : x = 3;
            break;
        default : x = 100;
    }
    return x;
}

double F_0(double t){
    return 0.5*sqrt(PI/t)*erf(sqrt(t));
}

double S(int p, int q){
    return pow(PI/(alpha(p)+alpha(q)), 1.5)*
    exp(-alpha(p)*alpha(q)/(alpha(p)+alpha(q))*(R(p) - R(q))*(R(p) - R(q)));
}

double T(int p, int q){
    return 0.5*pow(PI/(alpha(p)+alpha(q)), 1.5)*
    (alpha(p)*alpha(q)/(alpha(p)+alpha(q)))*
    (6-4*(alpha(p)*alpha(q)/(alpha(p)+alpha(q)))*(R(p) - R(q))*(R(p) - R(q)))*
    exp(-alpha(p)*alpha(q)/(alpha(p)+alpha(q))*(R(p) - R(q))*(R(p) - R(q)));
}

double A(int p, int q, int R_C){
    double R_p = (alpha(p)*R(p) + alpha(q)*R(q))/(alpha(p)+alpha(q));
    
    if (R(p) == R(q) && R(p) == R_C) {
        return -2*PI/(alpha(p)+alpha(q));
    }
    else {
    return -2*PI/(alpha(p)+alpha(q))*
    exp(-alpha(p)*alpha(q)/(alpha(p)+alpha(q))*(R(p) - R(q))*(R(p) - R(q)))*
        F_0((alpha(p)+alpha(q))*(R_p - R_C)*(R_p - R_C)); }
}

double H_Core(int p, int q, int Z){
    return T(p,q) + Z*A(p,q, R(0)) + Z*A(p,q, R(4));
}

MatrixXd Calculate_H(int Z){
    MatrixXd h(8,8);
    for(int p = 0; p<8; p++){
        for(int q = 0; q<8; q++){
            h(p,q) = H_Core(p,q,Z);
        }
    }
    return h;
}

MatrixXd Calculate_S(){
    MatrixXd s(8,8);
    for(int p = 0; p<8; p++){
        for(int q = 0; q<8; q++){
            s(p,q) = S(p,q);
        }
    }
    return s;
}


int main(int argc, const char * argv[]) {
    
    // 1. input Z and R
    int Z = 1;
    double epsilon;

    
    MatrixXd h(8,8);
    MatrixXd s(8,8);
    
    h = Calculate_H(Z);
    s = Calculate_S();
    
    //For Hydrogen - Single equation, no SCF
    GeneralizedEigenSolver<MatrixXd> ges;
    ges.compute(h, s);
    VectorXd::Index minVal;
    epsilon = ges.eigenvalues().real().minCoeff(&minVal);
    cout << "E_g = " <<  epsilon << endl;

    return 0;
}


