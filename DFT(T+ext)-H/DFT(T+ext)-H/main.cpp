using namespace std;

#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
//#include <Eigen/Eigenvalues>
//#include <Eigen/Core>
//#include <unsupported/Eigen/CXX11/Tensor>

//using namespace Eigen;
//using Eigen::MatrixXd;
//using Eigen::EigenSolver;

typedef double (*fn)(double);


double* Numerov_in(double u_max, double u_max_1, double h, int r_max, double epsilon){
    
    double* w = new double [r_max];
    double* u = new double [r_max];
    double* f = new double [r_max];
    
    for(int i = r_max; i>0; i--){
        double r = i*0.01;
        f[i] = -2/r-2*epsilon;
    }
    
    
    
    u[r_max] = u_max;
    u[r_max-1] = u_max_1;
    w[r_max] = (1-h*h*f[r_max]/12)*u[r_max];
    w[r_max-1] = (1-h*h*f[r_max-1]/12)*u[r_max-1];
    
    for(int i = r_max-1; i>0; i--){
        w[i-1] = 2*w[i]-w[i+1]+h*h*f[i]*u[i];
        u[i-1] = w[i-1]/(1-h*h*f[i-1]/12);
    }
    
    return u;
}



double Secant(double u_max, double u_max_1, double h, int r_max){
    
    double delta = 0.0000001;
    
    double* epsilon = new double [100];
    double* f = new double [100];
    
    epsilon[0] = -2;
    epsilon[1] = -0.4;
    
    f[0] = Numerov_in(u_max, u_max_1, h, r_max, epsilon[0])[0];
    f[1] = Numerov_in(u_max, u_max_1, h, r_max, epsilon[1])[0];
    
    int i = 1;
    while(abs(f[i]) > delta){
       // cout << epsilon[i] << " " << f[i-1] << " " << f[i] << " " << epsilon[i-1] << endl;
        epsilon[i+1] = (epsilon[i]*f[i-1] - epsilon[i-1]*f[i])/(f[i-1]-f[i]);
        f[i+1] = Numerov_in(u_max, u_max_1, h, r_max, epsilon[i+1])[0];
        cout << epsilon[i+1] << endl;
        i++;
    }
    return epsilon[i];
}




int main(int argc, const char * argv[]) {
    
    int r_max = 7;
    double h = 0.01;
    
    
    double u_max = r_max*exp(-r_max);
    double u_max_1 = (r_max-h)*exp(-(r_max-h));
    
    double epsilon = Secant(u_max, u_max_1, h, r_max/h);
    
    double *u = Numerov_in(u_max, u_max_1, h, r_max/h, epsilon);
    
    for(int i=0; i<700; i++){
        cout << i << " " << u[i] << endl;
    }


    
    return 0;
}
