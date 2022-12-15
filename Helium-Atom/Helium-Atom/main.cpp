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
/*double alpha(int i){
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
        default : x = 100;
    }
    return x;
}*/

double alpha(int i){
    double x;
    
    switch(i){
        case 0 : x = 0.298073;
            break;
        case 1 : x = 1.242567;
            break;
        case 2 : x = 5.782948;
            break;
        case 3 : x = 38.47497;
            break;
        default : x = 100;
    }
    return x;
}

double Gaussian(int i, double r){
    return exp(-alpha(i)*r*r);
}

double S(int p, int q){
    return pow(PI/(alpha(p)+alpha(q)), 1.5);
}

double T(int p, int q){
    return 3*alpha(p)*alpha(q)*pow(PI,1.5)/(pow(alpha(p)+alpha(q), 2.5));
}

double A(int p, int q){
    return -2*PI/(alpha(p)+alpha(q));
}

double H_Core(int p, int q, int Z){
    return T(p,q)+ Z*A(p,q);
}

double Q(int p, int q, int r, int s){
    return 2*pow(PI,2.5)/((alpha(p)+alpha(q))*(alpha(r)+alpha(s))*pow(alpha(p)+alpha(q)+alpha(r)+alpha(s), 0.5));
}

void Normalise(Vector4d &C, Matrix4d S){
    double sum = C.transpose()*S*C;
    C /= pow(sum, 0.5);
}

Matrix4d Calculate_f(Tensor<double, 4> qq, Matrix4d h, Vector4d C){
    
    Matrix4d f(4,4);
    
    for(int p = 0; p<4; p++){
        for(int q = 0; q<4; q++){
            f(p,q) = h(p,q);
            for(int r = 0; r<4; r++){
                for(int s = 0; s<4; s++){
                    f(p,q) += qq(p,q,r,s)*C(r)*C(s);
                }
            }
        }
    }
 
    return f;
}

double Calculate_Eg(Tensor<double, 4> qq, Matrix4d h, Vector4d C){

    double E_g = 2*C.transpose()*h*C;
    
    for(int p = 0; p<4; p++){
        for(int q = 0; q<4; q++){
            for(int r = 0; r<4; r++){
                for(int s = 0; s<4; s++){
                    E_g += qq(p,q,r,s)*C(p)*C(q)*C(r)*C(s);
                }
            }
        }
    }
    return E_g;
}


int main(int argc, const char * argv[]) {

    // 1. input Z and R
    int Z = 2;
    int R = 0;
    
    double epsilon = 100;
    double epsilon_new;
    double delta = 100;
    double precision = 0.000001;
    
    Matrix4d h(4,4);
    Matrix4d s(4,4);
    Matrix4d f(4,4);
    Tensor<double, 4> qq(4, 4, 4, 4);
    
    for(int p = 0; p<4; p++){
        for(int q = 0; q<4; q++){
            h(p,q) = H_Core(p,q,Z);
            s(p,q) = S(p,q);
            for(int r = 0; r<4; r++){
                for(int s = 0; s<4; s++){
                    qq(p,q,r,s) += Q(p,q,r,s);
                }
            }
        }
    }

    
    Vector4d C(1,1,1,1);
    
    
    //SCF - For Helium. Put Z = 2, and interate until convergence
    cout << "epsilon" << "                  " << "C" << endl;
    
    while (delta > precision){
    
    Normalise(C,s);
    
    f = Calculate_f(qq, h, C);
    
    GeneralizedEigenSolver<Matrix4d> ges;
    ges.compute(f, s);
    Vector4d::Index minVal;
    epsilon_new = ges.eigenvalues().real().minCoeff(&minVal);
    delta = abs(epsilon_new - epsilon);
    epsilon = epsilon_new;
    C = ges.eigenvectors().real().col(minVal);
    cout  << epsilon << "       " << C.transpose() <<  endl;

    }

    
    Normalise(C,s);
    cout << "E_g = " <<  Calculate_Eg(qq, h, C) << endl;
 
    
    //For Hydrogen - Single equation, no SCF
    // Change Z to 1, And change the values of alpha
   /* Normalise(C,s);
    GeneralizedEigenSolver<Matrix4d> ges;
    ges.compute(h, s);
    Vector4d::Index minVal;
    epsilon = ges.eigenvalues().real().minCoeff(&minVal);
    cout << "E_g = " <<  epsilon << endl;*/
    
    
    
    
    return 0;
}












/*void Normalise(double (&C)[4], double S[4][4]){
 double sum = 0;
 
 for(int p = 1; p<5; p++){
 for(int q = 1; q<5; q++){
 sum += C[p]*S[p][q]*C[q];
 }
 }
 for(int p = 1; p<5; p++)
 C[p] = C[p]/pow(sum, 1/2);
 }*/


/*
 double h[4][4];
 double s[4][4];
 double qq[4][4][4][4];
 double f[4][4];
 
 for(int p = 0; p<4; p++){
 for(int q = 0; q<4; q++){
 h[p][q] = H_Core(p,q,Z);
 s[p][q] = S(p,q);
 }
 }
 
 double C[4] = {1,1,1,1};
 Normalise(C, s);
 
 for(int p = 0; p<4; p++){
 for(int q = 0; q<4; q++){
 f[p][q] = h[p][q];
 for(int r = 0; r<4; r++){
 for(int s = 0; s<4; s++){
 qq[p][q][r][s] = Q(p,q,r,s);
 f[p][q] += qq[p][q][r][s]*C[r]*C[s];
 }
 }
 }
 }
 
 cout << s[0][0] << " " <<s[0][1] << " " <<s[0][2] << " " <<s[0][3] << endl;
 cout << s[1][0] << " " <<s[1][1] << " " <<s[1][2] << " " <<s[1][3] << endl;
 cout << s[2][0] << " " <<s[2][1] << " " <<s[2][2] << " " <<s[2][3] << endl;
 cout << s[3][0] << " " <<s[3][1] << " " <<s[3][2] << " " <<s[3][3] << endl;
 
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
 default : x = 100;
 }
 return x;
 }
 
 
 */
