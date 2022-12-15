#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;

const double PI  =3.141592653589793238463;

//------------------------------------------------------------------------------------

double* Numerov_in(double u_max, double u_max_1, double h, double r_max, double epsilon, int Z, double* V_H, double* V_x){
    
    int n_max = r_max/h;
    
    double* w = new double [n_max+1];
    double* u = new double [n_max+1];
    double* f = new double [n_max+1];
    
    
    for(int i = n_max; i>0; i--){
        double r = i*h;
        f[i] = 2*(-Z/r -epsilon + V_H[i] + V_x[i]);
        
    }
    
    //Boundary Conditions
    u[n_max] = u_max;
    u[n_max-1] = u_max_1;
    w[n_max] = (1-h*h*f[n_max]/12)*u[n_max];
    w[n_max-1] = (1-h*h*f[n_max-1]/12)*u[n_max-1];
    
    for(int i = n_max-1; i>0; i--){
        w[i-1] = 2*w[i]-w[i+1]+h*h*f[i]*u[i];
        u[i-1] = w[i-1]/(1-h*h*f[i-1]/12);
    }
    //At r=0, f(0) diverges, and in above, f[0] = 0 which is wrong. So manually set u[0] = w[0]
    u[0] = w[0];
    return u;
}

//------------------------------------------------------------------------------------

double Bisection(double u_max, double u_max_1, double h, double r_max, int Z, double* V_H, double* V_x, double a, double b){
    
    double delta = 0.000000001;
    double epsn = 0.000000001;
    double maxit = 100;
    
    double u = Numerov_in(u_max, u_max_1, h, r_max, a, Z, V_H, V_x)[0];
    double e = b-a;
    double c;
    
    for(int k=1; k <= maxit; k++){
        
        e*= 0.5;
        c = a + e;
        double w = Numerov_in(u_max, u_max_1, h, r_max, c, Z, V_H, V_x)[0];

        if (fabs(e) < delta || fabs(w) < epsn) return c;
        ((u>0 && w<0) || (u<0 && w>0)) ? (b=c) : (a=c , u=w);
        
    }
    
    return c;
}

//------------------------------------------------------------------------------------

void Normalise(double* &u, double r_max, double h){
    
    int n_max = r_max/h;
    double N;
    double* f = new double [n_max+1];
    
    for(int i =0; i <= n_max; i++)
        f[i] = u[i]*u[i];
    
    //Integrate u^2(r) to get N. Then divide u(r) by Sqrt[N] to normalise it.
    double sum = (f[0]+f[n_max])*0.5;
    for(int i =1; i < n_max; i++)
        sum+= f[i];
    
    N = sum*h;
    
    for(int i =0; i <= n_max; i++)
        u[i] = u[i]/sqrt(N);
    
}

//------------------------------------------------------------------------------------

double* Hartree_Potential(double U_0, double U_1, double h, double r_max, int Z, double u_max, double u_max_1, double epsilon, double* u){
    
    int n_max = r_max/h;
    double* U = new double [n_max+1];
    double* V_H = new double [n_max+1];
    
    //Boundary Conditions
    U[0] = U_0;
    U[1] = U_1;
    
    //Verlet Algorithm to integrate U
    for(int i = 1; i<n_max; i++){
        double r = i*h;
        U[i+1] = 2*U[i] - U[i-1] - h*h*u[i]*u[i]/r;
    }
    
    //The integration constant is alpha*r. Fix alpha such that U(r_max) = Z/r_max
    double alpha = (Z - U[n_max-1])/(n_max*h);
    
    //Add the integration constant to U, and calculate V_H as V_H = 2*U/r.
    //The factor of 2 is because u was the single-orbital wavefunction, but we want the full Hartree potential. For Helium, since both elctrons are in 1s, we can simply multiply by 2.
    for(int i = 0; i<=n_max; i++){
        U[i] += alpha*i*h;
        V_H[i] = 2*U[i]/(h*i);
    }
    
    return V_H;
}

//------------------------------------------------------------------------------------

double Hartree_Energy(double r_max, double h, double* u, double* V_H){
    
    int n_max = r_max/h;
    
    //Use Trapezoidal algorithm to integrate V_H*u^2(r)
    double* f = new double [n_max+1];
    for(int i =1; i <= n_max; i++)
        f[i] = V_H[i]*u[i]*u[i];
    
    f[0] = 0;
    double sum = (f[0]+f[n_max])*0.5;
    for(int i =1; i < n_max; i++)
        sum+= f[i];
    
    return sum*h;
}

//------------------------------------------------------------------------------------

double* Exchange_Potential(double r_max, double h, double* u){
    
    int n_max = r_max/h;
    double* V_x = new double [n_max+1];
    
    //Calculate V_x.
    for(int i =1; i <= n_max; i++)
        V_x[i] = -pow(3*u[i]*u[i]/(2*PI*PI*i*h*i*h), 1.0/3);
    
    //To avoid division by zero, set V_x(0) = V_x(h)
    V_x[0] = V_x[1];
    
    return V_x;
}

//------------------------------------------------------------------------------------

double Exchange_Energy(double r_max, double h, double* u, double* V_x){
    
    int n_max = r_max/h;
    
    //Use Trapezoidal algorithm to integrate V_x*u^2(r)
    double* f = new double [n_max+1];
    for(int i =0; i <= n_max; i++)
        f[i] = V_x[i]*u[i]*u[i];
    
    
    double sum = (f[0]+f[n_max])*0.5;
    for(int i =1; i < n_max; i++)
        sum+= f[i];
    
    return sum*h;
    
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------


int main(int argc, const char * argv[]) {
    
    double r_max = 50;
    double h = 0.01;
    int Z = 2;
    double epsilon;
    
    double E = 100;
    double E_new;
    double precision = 0.000001;
    double delta = 100;
    
    double u_max = r_max*exp(-r_max);
    double u_max_1 = (r_max-h)*exp(-(r_max-h));
    
    double* V_H = new double [r_max/h+1];
    double* u = new double [r_max/h+1];
    double* V_x = new double [r_max/h+1];
    
    for(int i = 0; i<=r_max/h; i++){
        V_H[i] = 0.0;
        V_x[i] = 0.0;
    }
    
    
    while(abs(delta) > precision){
        
        epsilon = Bisection(u_max, u_max_1, h, r_max, Z, V_H, V_x, -3, 0);
        u = Numerov_in(u_max, u_max_1, h, r_max, epsilon, Z, V_H, V_x);
        Normalise(u, r_max, h);
        V_H = Hartree_Potential(0, h, h, r_max, Z, u_max, u_max_1, epsilon, u);
        V_x = Exchange_Potential(r_max, h, u);
        
        E_new = 2*epsilon - Hartree_Energy(r_max, h, u, V_H) - 0.5*Exchange_Energy(r_max, h, u, V_x);
        cout << epsilon << " " << Hartree_Energy(r_max, h, u, V_H) << " " << Exchange_Energy(r_max, h, u, V_x) << " " << E_new << endl;
        
        delta = E_new - E;
        E = E_new;
    }
    
    return 0;
}









/*
 double Secant(double u_max, double u_max_1, double h, double r_max, int Z, const double* const V_H, double* V_x){
 
 
 double delta = 0.0000001;
 
 double* epsilon = new double [100];
 double* f = new double [100];
 
 epsilon[0] = -1.2;
 epsilon[1] = -0.8;
 
 f[0] = Numerov_in(u_max, u_max_1, h, r_max, epsilon[0], Z, V_H, V_x)[0];
 f[1] = Numerov_in(u_max, u_max_1, h, r_max, epsilon[1], Z, V_H, V_x)[0];
 
 int i = 1;
 while(abs(f[i]) > delta){
 epsilon[i+1] = (epsilon[i]*f[i-1] - epsilon[i-1]*f[i])/(f[i-1]-f[i]);
 f[i+1] = Numerov_in(u_max, u_max_1, h, r_max, epsilon[i+1], Z, V_H, V_x)[0];
 // cout << epsilon[i+1] << endl;
 i++;
 }
 return epsilon[i];
 }
 */


/*
 ofstream myFile;
 myFile.open("test.csv");
 for(int i=0; i<700; i++){
 cout << i << " " << U[i] << endl;
 myFile << i << "," << U[i] << endl;
 }*/
