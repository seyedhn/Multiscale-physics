using namespace std;

#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace Eigen;

const double PI  =3.141592653589793238463;


//---------------------------------------------------------------------

//For Hydrogen
// The constants in the Gaussian basis. We use 8 bases in total, 4 positioned around each Hydrogen atom.
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

//---------------------------------------------------------------------

//The position of the basis functions around either of the two nuclei
int R(int i){
    
    int x;
    
    switch(i){
        case 0 ... 3 : x = 1;
            break;
        case 4 ... 7 : x = 0;
            break;
        default : x = 100;
    }
    return x;
}

//---------------------------------------------------------------------

//The function F_0 which appears in the nuclear attraction A and electron repulsion Q.
// F_0(0) = 1 according to L'Hopital rule. Implement this manually.
double F_0(double t){
    
    if (t <= 0 && t >= -0.01)
        return 1;
    else
        return 0.5*sqrt(PI/t)*erf(sqrt(t));
}

//---------------------------------------------------------------------

//The overlap matrix of two Gaussian basis functions.
double S(int p, int q){
    
    return pow(PI/(alpha(p)+alpha(q)), 1.5)*
    exp(-alpha(p)*alpha(q)/(alpha(p)+alpha(q))*(R(p) - R(q))*(R(p) - R(q)));
}

//---------------------------------------------------------------------

//The Kinetic Energy term
double T(int p, int q){
    
    return 0.5*pow(PI/(alpha(p)+alpha(q)), 1.5)*
    (alpha(p)*alpha(q)/(alpha(p)+alpha(q)))*
    (6-4*(alpha(p)*alpha(q)/(alpha(p)+alpha(q)))*(R(p) - R(q))*(R(p) - R(q)))*
    exp(-alpha(p)*alpha(q)/(alpha(p)+alpha(q))*(R(p) - R(q))*(R(p) - R(q)));
}

//---------------------------------------------------------------------

//The nuclear attraction term
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

//---------------------------------------------------------------------

//The core Hamiltonian is the sum of Kinetic energy and nuclear attractions
double H_Core(int p, int q, int Z){
    
    return T(p,q) + Z*A(p,q, R(0)) + Z*A(p,q, R(4));
}

//---------------------------------------------------------------------

//The Q matrix of electron-electron repulsion term
double Q(int p, int q, int r, int s){
    
    double A = alpha(p)+alpha(q);
    double B = alpha(r)+alpha(s);
    double R_A = (alpha(p)*R(p) + alpha(q)*R(q))/(alpha(p)+alpha(q));
    double R_B = (alpha(r)*R(r) + alpha(s)*R(s))/(alpha(r)+alpha(s));
    double t = (alpha(p)+alpha(q))*(alpha(r)+alpha(s))/(alpha(p)+alpha(q)+alpha(r)+alpha(s))*(R_A-R_B)*(R_A-R_B);
    
    return 2*sqrt(A*B/(PI*(A+B)))*S(p,q)*S(s,r)*F_0(t);
}

//---------------------------------------------------------------------

//This function normalises the eigenvector C
void Normalise(VectorXd &C, MatrixXd S){
    
    double sum = C.transpose()*S*C;
    C /= pow(sum, 0.5);
}

//---------------------------------------------------------------------

//This function calculates the overlap matrix S, the core Hamiltonian H_core and the exchange matrix Q.
void Calculate_S_H_Q(Tensor<double, 4> &qq, MatrixXd &h, MatrixXd &s, int Z){
    
    for(int p = 0; p<8; p++){
        for(int q = 0; q<8; q++){
            h(p,q) = H_Core(p,q,Z);
            s(p,q) = S(p,q);
            for(int r = 0; r<8; r++){
                for(int s = 0; s<8; s++){
                   // qq(p,q,r,s) = Q(p,q,r,s);
                   qq(p,q,r,s) = 2*Q(p,q,r,s) - Q(p,s,r,q); //The Exchange energy
                }
            }
        }
    }
}

//---------------------------------------------------------------------

//Calculate the Fock matrix from all its components H_core and Q
MatrixXd Calculate_f(Tensor<double, 4> qq, MatrixXd h, VectorXd C){
    
    MatrixXd f(8,8);
    
    for(int p = 0; p<8; p++){
        for(int q = 0; q<8; q++){
            f(p,q) = h(p,q);
            for(int r = 0; r<8; r++){
                for(int s = 0; s<8; s++){
                    f(p,q) += qq(p,q,r,s)*C(r)*C(s);
                }
            }
        }
    }
    
    return f;
}

//---------------------------------------------------------------------

//Calculate the total energy from H_core + 1/2(J-K) + H_NN
double Calculate_Eg(Tensor<double, 4> qq, MatrixXd h, VectorXd C){
    
    double E_g = 2*C.transpose()*h*C;
    
    for(int p = 0; p<8; p++){
        for(int q = 0; q<8; q++){
            for(int r = 0; r<8; r++){
                for(int s = 0; s<8; s++){
                    E_g += qq(p,q,r,s)*C(p)*C(q)*C(r)*C(s);
                }
            }
        }
    }
    return E_g + 1; //+1 is for the nuclear repulsion
}

//---------------------------------------------------------------------

//Alternatively, we can also calculate total energy as 1/2(sum_i epsilon_i + H_core)
double Calculate_Eg2(MatrixXd h, VectorXd C, double epsilon){
    
    return C.transpose()*h*C + epsilon + 1; //+1 is for the nuclear repulsion Z^2/R
}

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------


int main(int argc, const char * argv[]) {
    
    // 1. input Z and R
    int Z = 1;
    
    double epsilon = 100;
    double epsilon_new;
    double delta = 100;
    double precision = 0.0000001;
    
    // 2. Determine the matrices S, H, and Q
    MatrixXd h(8,8);
    MatrixXd s(8,8);
    MatrixXd f(8,8);
    Tensor<double, 4> qq(8, 8, 8, 8);
    
    Calculate_S_H_Q(qq, h, s, Z);
    

    // 3. Make an initial guess for C
    VectorXd C(8);
    C << 1,1,1,1,1,1,1,1;
    VectorXd eps(8);
    
    //4. SCF
    cout << "epsilon" << "                  " << "C" << endl;
    while (delta > precision){
        
        Normalise(C,s);
        
        f = Calculate_f(qq, h, C);
    
        GeneralizedEigenSolver<MatrixXd> ges;
        ges.compute(f, s); //Solve the eigensystem
        VectorXd::Index minVal;
        eps = ges.eigenvalues().real(); //calculate the eigenvalues eps
        epsilon_new = ges.eigenvalues().real().minCoeff(&minVal); //The lowest eigenvalue corresponds to the ground state energy
        delta = abs(epsilon_new - epsilon);
        epsilon = epsilon_new;
        C = ges.eigenvectors().real().col(minVal); //Calculate the eigenvectors C
        
        cout  << epsilon << "       " << C.transpose() <<  endl;
    }
    
    // 5. Calculate and print Energy
    Normalise(C,s);
    cout << endl << C.transpose() << endl;
    cout << "E_g = " <<  Calculate_Eg(qq, h, C) << endl;
    cout << "E_g2 = " <<  Calculate_Eg2(h, C, epsilon) << endl;
    
    return 0;
}



