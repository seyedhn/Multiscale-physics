#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <gmpxx.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
//#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <gnuplot-iostream.h>


using namespace Eigen;
using namespace std;

const double KB = 1.38e-23;
//------------------------------------------------------------------------------------

vector<string> getData(string FileName){
    
    ifstream file(FileName);
    vector<string> dataList;
    
    string line;
    // Iterate through each line and split the content using delimeter
    while (getline(file, line, ','))
    {
        dataList.push_back(line);
    }
    
    return dataList;
}

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


//************************************************************************************
//Methods
//************************************************************************************

//Checks if a histogram is flat within some precision.
bool isFlat(unordered_map<int,long double> histogram, double precision){
    
    int size = histogram.size();
    double a = histogram[size/2];
    double b = histogram[size/2];
    
    for(int i=0; i<size; i++){
        
        if(histogram[i] < a) a = histogram[i];
        if(histogram[i] > b) b = histogram[i];
        
    }
    
    return (b-a)/(a+0.0001) < precision;

}

//------------------------------------------------------------------------------------

/*Performs a Wang-Landau sampling. Starts with a random configuration, and performs a random
 walk based on a rule. Keeps record of number of times an energy state is visited.*/
void WLSampling(int L, double E_0, double E_N, double bin_width, int INNER_STEPS, int ITERS){
    
    ofstream myfile ("Histogram.csv");
    
    int N = L*L;
    double n_bins = (E_N-E_0)/bin_width; //number of bins
    double f = 2.7;                     //the initial control factor
    int step;
    double E1;                          //energy of the current configuration
    double E2;
    int binE1;
    int binE2;
    //energy of the proposed configuration
    unordered_map<int,long double> lnOmega;
    unordered_map<int,long double> hist;
    
    //initiate the log of DoS Omega. We set all Omega(E) = 1, so log(Omega(E)) = 0 for all energies.
    for(int i =-n_bins/2; i<=n_bins/2; i++)
        lnOmega[i] = 0;

    
    Lattice lat1(L, 1, 0.01);   //Create an initial lattice
    Lattice lat2(L, 1, 0.01);   //Create a proposed lattice
    
    E1 = lat1.Energy()/N;
    binE1 = trunc(E1/bin_width);
    
    int iter = 0;
    while(iter < ITERS){
        
        for(int i =-n_bins/2; i<=n_bins/2; i++)
            hist[i] = 0;        //Keep a histogram of all the visited sites.
        
        step = 1;
        while(step < INNER_STEPS){
            
            int a = rand() % L;
            int b = rand() % L;
            
            lat2 = lat1;
            lat2.flip(a,b);     //Make a random walk from lat1
            
            E2 = double (E1 - 2*lat1.LocalEnergy(a, b)/N);
            binE2 = trunc(E2/bin_width);
            
            cout << f << " " << iter << " " << step << " " << binE1 << " " << binE2 << endl;
            
            //Make a transition to the new state according to p(E1 -> E2) = min (1, Omega(E1) / Omega(E2))
            //Since we're dealing with log(Omega), we subtract the two values.
            if (( log((double) rand() / (RAND_MAX))) <= (double) (lnOmega[binE1] - lnOmega[binE2]) ){
                lat1 = lat2;
                E1 = E2;
                binE1 = binE2;
            }

            lnOmega[binE1] += log(f);       //Update Omega(E) by f*Omega(E). In log, we need to add.
            hist[binE1] += 1;               //Add one to the histogram of visited sites.
            
            step++;
        }
        
        //Once the histogram is flat enough, update f, and repeat everything again (histogram is reset as well).
        //Make consequent iterations slightly longer by increasing the INNER_STEPS.
        f = sqrt(f);
        iter++;
        INNER_STEPS += 50000;
        
    }
    
    for(double i=-n_bins/2; i< n_bins/2; i++){
        myfile << i << "," << hist[i] << "," << lnOmega[i] << endl;
    }

    
}


//------------------------------------------------------------------------------------

void SpecificHeat(double E_0, double E_N, double E_bin_width, double T_0, double T_N, double T_bin_width){
    
    ofstream myfile ("HeatCapacity.csv");
    vector<string> data = getData("Histogram.csv");
    
    double Z = 0;
    double ExpE = 0;
    double ExpE2 = 0;
    double C_v;
    double E_count;
    unordered_map<int,long double> lnOmega;
    
    double j = trunc(E_0/E_bin_width);
    for (int i = 1; i <= data.size(); i+=2){
        lnOmega[j] =  atof(data[i].c_str());
        j++;
    }
    
    //for(int T = trunc(T_0/T_bin_width); T<=trunc(T_N/T_bin_width); T++){
    for(double T = T_0; T<=T_N; T+= T_bin_width){
        
        Z = 0;
        ExpE = 0;
        ExpE2 = 0;
        
        //mpf_t ZZ;
        //mpf_init2(ZZ, 10);
        //mpf_init_set_ui(ZZ, 2);
        for(double E = E_0; E<=E_N; E+= E_bin_width){
        
            E_count = trunc(E/E_bin_width);
            
           // mpf_init_set_ui(ZZ, exp(lnOmega[E_count])*exp(-E/T));
          //  mpf_add_ui(ZZ, ZZ, exp(lnOmega[E_count])*exp(-E/T));
            
            Z += exp(lnOmega[E_count])*exp((-E-2)/T);
            ExpE += E*exp(lnOmega[E_count])*exp((-E-2)/T);
            ExpE2 += E*E*exp(lnOmega[E_count])*exp((-E-2)/T);
            cout << E << " " << exp(lnOmega[E_count]) << " " << Z << " " << (ExpE/Z)*(ExpE/Z) << " " << ExpE2/Z << endl;
        
        }
        C_v = (ExpE2/Z - ExpE*ExpE/(Z*Z))/(T*T);
        
        cout << T << "," <<  Z << "," << (ExpE/Z)*(ExpE/Z) << "," << ExpE2/Z << "," << C_v << endl;

        myfile << T << "," <<  Z << "," << ExpE << "," << ExpE2 << "," << -T*log(Z) << "," <<  C_v << " " <<  endl;
        //myfile << T << "," <<  -T*log(Z) <<  endl;

    }
    
}



//************************************************************************************
//Main
//************************************************************************************

int main(int argc, const char * argv[]) {

    clock_t start, end;
    double cpu_time_used;
    start = clock();
    

    SpecificHeat(-2, 2, 0.05, 0.0001, 0.02, 0.0001);
   // WLSampling(16, -2, 2, 0.05, 200000, 20);
    

    
    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout <<  cpu_time_used << endl;
    
    return 0;
}
