#ifndef _LATTICE_H
#define _LATTICE_H



#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <gmp.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
//#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
//#include <gnuplot-iostream.h>


using namespace Eigen;
using namespace std;

const double KB = 1.38e-23;
//------------------------------------------------------------------------------------

//********************************************************************************
//Classes
//********************************************************************************

/*The class Lattice defines objects for a 2D square Ising model lattice with length L, total number of spins L*L,
and spins facing up (1) or down (-1).*/
class Lattice{
    
private:

    int L;
    int N;
    double J;
    double T;
    MatrixXd lattice;
    
public:
    
    Lattice(int length, double strength, double temperature){
        
        L = length; //Length of the lattice
        N = L*L;    //Number of spins
        J = strength;   //The coupling strength
        T = temperature;
        lattice.resize(L, L); //This is necessary when using MatrixXd objects which defines its size.

        //Creates a lattice with random spins.
        int random;
        for(int i =0; i<L; i++){
            for(int j =0; j<L; j++){
                
                random = ((rand() % 2)-0.5)*2;
                lattice(i,j) = random;
            }
        }

    }

//------------------------------------------------------------------------------------
    
    void set_T(double temp){
        T = temp;
    }

//------------------------------------------------------------------------------------
    
    //Resets the lattice to another random configuration.
    void reset(){
        
        int random;
        for(int i =0; i<L; i++){
            for(int j =0; j<L; j++){
                
                random = ((rand() % 2)-0.5)*2;
                lattice(i,j) = random;
            }
        }
        
    }

//------------------------------------------------------------------------------------
    /*Calculates the local energy of each spin site through the formula E_i = -J*s_i*sum_j s_j
     The algorithm is based on periodic boundary conditions.*/
      double LocalEnergy(int i, int j){
          
          int Xi = i-1;
          int Xj = j-1;
          if(i == 0) Xi = L-1;
          if(j == 0) Xj = L-1;
          
         return -J*lattice(i,j)*(lattice(i,Xj) + lattice(Xi,j) + lattice((i+1)%L,j) + lattice(i,(j+1)%L));

      }

//------------------------------------------------------------------------------------
    
    //Claculates the total energy through adding all the local energies and dividing by 2 to avoid double counting.
    double Energy(){
            double En = 0;
            for(int i =0; i<L; i++){
                for(int j =0; j<L; j++){
                    
                    En += LocalEnergy(i, j);
                }
            }
            return En/2;
    }

//------------------------------------------------------------------------------------
    
    //Calculates the total magnetism by adding all the spins.
    double magnetism(){
        double M = 0;
        for(int i =0; i<L; i++){
            for(int j =0; j<L; j++){
                
                M += lattice(i, j);
            }
        }
        return abs(M);
    }
    
//------------------------------------------------------------------------------------
    //flips the spin s_{ij}
    void flip(int i, int j){
        lattice(i,j) = -lattice(i,j);
    }

 //------------------------------------------------------------------------------------
    
    //Performs a Metropolis Monte Carlo algorithm.
    void MMC(int MAX_STEPS){
        
        int step = 1;
        int success_counts = 0;
        while(step <= MAX_STEPS){
            
            int a = rand() % L;
            int b = rand() % L;


            if (((double) rand() / (RAND_MAX)) <= exp(2/T*(LocalEnergy(a, b)))){
                flip(a,b);
                success_counts++;
            }

            step++;
        }
        
        cout << "The acceptance rate is: " << (double) success_counts / MAX_STEPS << " ";
        
    }

//------------------------------------------------------------------------------------
    
    //Prints the lattice
    void print_lattice(){
        
        for(int i =0; i<L; i++){
            for(int j =0; j<L; j++){
               
                if(lattice(i,j) == 1) cout << 'o';
                else cout << 'x';

            }
            cout << endl;
        }

    }
 
//------------------------------------------------------------------------------------
    //Performs MMC for a range of different temperatures
    void T_range(double a, double b, double c, int steps){
        
        double n = (b-a)/c;

        for(int i =0; i<=n; i++){
                
                 reset();
                 set_T((double)(a+i/n*(b-a)));
                 MMC(steps);

            cout << T << " " << magnetism() << endl;
        
        }
    }
 
    
};




int haha(int a);
void SpecificHeat(double E_0, double E_N, double E_bin_width, double T_0, double T_N, double T_bin_width);

#endif
