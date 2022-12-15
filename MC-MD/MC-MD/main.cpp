#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <gnuplot-iostream.h>


using namespace Eigen;
using namespace std;

const double PI  = 3.141592653589793238463;


//************************************************************************************
//Methods
//************************************************************************************

//This method returns a random number based on the Gaussian distribution.
double GaussRand(double mu, double sigma){
    
    thread_local double y1;
    thread_local bool generate; //Initial value is false
    generate = !generate;   //Invert the value, so if we generated last time, we use stored value
    //otherwise, we generate
    
    if(!generate)           //If we generated last time, use the stored value.
        return y1*sigma + mu;
    
    double x1, x2;
    x1 = rand() * 1.0/RAND_MAX;
    x2 = rand() * 1.0/RAND_MAX;
    
    double y2;
    y1 = sqrt(-2*log(x1))*cos(2*PI*x2);
    y2 = sqrt(-2*log(x1))*sin(2*PI*x2);
    
    return y2*sigma + mu;
}

//------------------------------------------------------------------------------------

/*The function returns the Lennard-Jones Force Field based on the distance r between two particles.
 An additional division by r is made because we later multiply this by the direction vector between
 the two atoms. We have taken 'sigma' in LJ as the unit of distance, 'epsilon/K_B' as unit of
 temperaute, and sqrt(m*sigma^2/epsilon) as unit of time. */
double LJ(double r){
    
    return 48*pow(r, -14) - 24*pow(r,-8);
    
}

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

/*The class Atom contains information about each individual particle. The private members are position, velocity,
 and the total force acted upon it by all other atoms at any instance. */
class Atom{
    
private:
    
    Vector3d r;     //The position of the atom
    Vector3d v;     //The velocity of the atom
    Vector3d F;     //Sum of forces acting on the particle
    vector<int> NL; //Neighbour List of the closest atoms
    
public:
    
    //Constructor
    Atom(Vector3d pos){
        r = pos;
        v << 0,0,0;
    }
    
    //Position functions
    Vector3d get_r(){
        return r;
    }
    
    double get_r(int n){
        return r(n);
    }
    
    void shift_r(int i, double delta){
        r(i) += delta;
    }
    
    void shift_r(Vector3d delta){
        r += delta;
    }
    
    //Velocity functions
    void set_v(double x, double y, double z){
        v << x,y,z;
    }
    
    Vector3d get_v(){
        return v;
    }
    
    double get_v(int n){
        return v(n);
    }
    
    void shift_v(Vector3d delta){
        v += delta;
    }
    
    void scale_v(double lambda){
        v *= lambda;
    }
    
    //Force functions
    Vector3d get_F(){
        return F;
    }
    
    double get_F(int n){
        return F(n);
    }
    
    void clear_F(){
        F << 0,0,0;
    }
    
    void add_to_F(Vector3d delta){
        F += delta;
    }
    
    void clear_NL(){
        NL.clear();
    }
    void add_to_NL(int j){
        NL.push_back(j);
    }
    
    int get_NL_size(){
        return NL.size();
    }
    
    int get_NL(int j){
        return NL[j];
    }
};


//********************************************************************************

/*The class Cell contains all the information about the simulation environment, as well as all the functions
 for setting up and running the simulation, as well as gathering physical quantities from the simulation.
 The private members of the class are quantities such as cell size, lattice paramter, number of particles,
 density, volume, temperature, and time step. Note that many of these quantities are redundant, but they
 keep the code tidy and easier to track quantities.*/
class Cell{
    
private:
    
    double a;   //Lattice constant
    int M;      //An integer determining the number of particles
    double T;   //Temperaute,
    double h;   //Time step
    double sigma; //Standard deviation in velocity distribution = K_B*T/m
    vector<Atom> atom;  //A list of all atoms in the Cell.
    
    //Redundant quantities
    int N;      //number of particles, related to M through N = 4M^3
    double L;   //cell length, L = Ma
    double V;   //Volume of the cell, V = L^3
    double rho; //Density, rho = N/V
    
    
    
public:
    
    //Constructor
    Cell(double lattice, int box, double temperature, double timestep){
        
        a = lattice;
        M = box;
        T = temperature;
        h = timestep;
        
        L = M*a;
        V = L*L*L;
        N = 4*pow(M, 3);
        rho = N/V;
        
    }
    
    //------------------------------------------------------------------------------------
    // Initialise
    //------------------------------------------------------------------------------------
    
    /*The function sets up the simulation environment, assigning positions and velocities
     to each particle.*/
    void Initialise(){
        
        FCC();
        Initialise_Velocities();
        Update_Forces_NL();
        Update_NL(2.5);
        
    }
    
    //------------------------------------------------------------------------------------
    // Run
    //------------------------------------------------------------------------------------
    
    /*The function let's the particles interact. The positions, velocities and forces
     are updated at each step. Velocities are scaled and measurements are taken after
     every few steps. */
    void Run(int steps){
        
        
        for(int i=0; i<=steps; i++){
            
            Update_r_F_v();
            
            if(i%20 == 0){
                Rescale_Velocities();
                print_Values(i);
                Update_NL(2.5);
            }
        }
        
    }
    
    //------------------------------------------------------------------------------------
    // Measure
    //------------------------------------------------------------------------------------
    
    /*The function measures physical quantities when the system has reach equilibrium.*/
    void Measure(){
        
        Correlation_Function();
        //MSD(100);
    }
    
    
    
    //------------------------------------------------------------------------------------
    // Position functions
    //------------------------------------------------------------------------------------
    
    /*The function FCC creates an initial FCC arrangement of particles.*/
    void FCC(){
        
        double A = 1;        //The Atomic mass of the atoms. We take it as the unit of mass.
        sigma = sqrt(T/A);  //sigma is the standard deviation in Boltzmann distribution of velocities.
        
        //These are the Bravais Lattice Primitive Vectors for the FCC structure.
        Vector3d a0(0,0,0), a1(0.5, 0.5, 0), a2(0.5, 0, 0.5), a3(0, 0.5, 0.5);
        
        /*Assign positions to atoms according to primitive vectors.
         First, assign atoms to each primivite vector, then translate
         by R, and do the same. We end up with 4M^3 lattice positions. */
        for (int x =0; x<M; x++){
            for (int y =0; y<M; y++){
                for (int z =0; z<M; z++){
                    
                    Vector3d R(x, y, z); //Discrete translation operator
                    
                    atom.push_back(Atom(a*(R + a0)));
                    atom.push_back(Atom(a*(R + a1)));
                    atom.push_back(Atom(a*(R + a2)));
                    atom.push_back(Atom(a*(R + a3)));
                }
            }
        }
        
    }
    
    //------------------------------------------------------------------------------------
    
    //Calculates the shortest Position Vector between two particles in a periodic boundary condition.
    Vector3d Pair_Vector(int i, int j){
        
        Vector3d n, shift;
        n << 0,0,0;
        shift << L, L, L;
        
        /*Checks if atoms i and j are less than L/2 apart in any direction.
         If not, we calculate the force between i and the copy of j in a neighbouring
         copy of the cell which is less than L/2 from i. The cutoff L/2 ensures that the
         particle is either in the same cell, or a neighboring cell.*/
        for(int k =0; k<=2; k++)
            n(k) = trunc((atom[i].get_r(k) - atom[j].get_r(k))/(L/2));
        
        return atom[i].get_r() - (atom[j].get_r() + shift.cwiseProduct(n));
        
    }
    
    //------------------------------------------------------------------------------------
    
    /*Checks if an atom has left the box, and will return it back according to periodic BC.
     This method can be made more efficient by checking only the atoms walking near to the boundary.*/
    void Bring_back_to_box(int i){
        
        for(int k=0; k<=2; k++){
            if (atom[i].get_r(k) > L) atom[i].shift_r(k,-L);
            if (atom[i].get_r(k) < 0) atom[i].shift_r(k,L);
        }
        
    }
    
    //------------------------------------------------------------------------------------
    
    /*Calculates the total force on each particle due to all other particles.
     Is based on a Neighbours List technique*/
    void Update_NL(double r_max){
        
        for(int i=0; i<N; i++){
            
            atom[i].clear_NL();
            
            for(int j=0; j<N; j++){
                if (j == i) continue;
                
                if(Pair_Vector(i, j).norm() < r_max)
                    atom[i].add_to_NL(j);
            }
        }
        
    }
    //------------------------------------------------------------------------------------
    // Velocity functions
    //------------------------------------------------------------------------------------
    
    /* This method assigns random velocities to particles according to the
     Boltzmann distribution. The average aggregate velocity is then deducted
     from each velocity to ensure the total momentum is zero.*/
    void Initialise_Velocities(){
        
        //Set the momenta from the Boltzmann distribution. Use GaussRand to generate random speeds.
        for(int i = 0; i < N; i++)
            atom[i].set_v(GaussRand(0, sigma), GaussRand(0, sigma), GaussRand(0, sigma));
        
        Vector3d average = Calculate_Average_Velocity();
        
        //Subtract the average momenta per particle from the momenta of all atoms.
        for(int i = 0; i < N; i++)
            atom[i].shift_v(-average);
        
    }
    
    //------------------------------------------------------------------------------------
    
    /*This method calculates the average velocity of all particles.
     If the system contains atoms of different masses, we need to calculate average momentum.*/
    Vector3d Calculate_Average_Velocity(){
        
        Vector3d average;
        double sum_x = 0;
        double sum_y = 0;
        double sum_z = 0;
        
        //Calculate the average momenta in all 3 directions.
        for(int i = 0; i < N; i++){
            sum_x += atom[i].get_v(0);
            sum_y += atom[i].get_v(1);
            sum_z += atom[i].get_v(2);
        }
        average << sum_x/N, sum_y/N, sum_z/N;
        return average;
    }
    
    //------------------------------------------------------------------------------------
    
    /*Rescales velocities according to 3/2NKT = 1/2mv^2 in order to keep the temperature
     close to the desired value. */
    void Rescale_Velocities(){
        
        double sum = 0;
        double lambda;
        
        for(int i=0; i<N; i++)
            sum += pow(atom[i].get_v().norm(), 2);
        
        lambda = sigma*sqrt((N-1)*3.0/sum);
        
        for(int i=0; i<N; i++)
            atom[i].scale_v(lambda);
        
    }
    //------------------------------------------------------------------------------------
    
    /*Reverse the velocities at some point as a check to see if we arrive
     at the initial configuration.*/
    void Reverse_Velocities(){
        
        for(int i=0; i<N; i++)
            atom[i].scale_v(-1);
    }
    
    
    //------------------------------------------------------------------------------------
    // Force functions
    //------------------------------------------------------------------------------------
    
    /*This function calculates and updates r, F and v according to the Verlet's algorithm.*/
    void Update_r_F_v(){
        
        Vector3d temp_F[N];
        
        /*Updates positions according to v_n & F_n.
         We use velocity-Verlet's algorithm.*/
        for(int i=0; i<N; i++){
            
            atom[i].shift_r(h*atom[i].get_v() + h*h*atom[i].get_F()/2);
            temp_F[i] = atom[i].get_F();                                //Store F_n
            Bring_back_to_box(i);                                       //Check if atom has left box
            
        }
        
        /*Update Forces according to new positions r_n+1.*/
        Update_Forces_NL();
        
        //Update velocities with F_n and F_n+1
        for(int i=0; i<N; i++)
            atom[i].shift_v(h*(atom[i].get_F() + temp_F[i])/2);
        
    }
    
    //------------------------------------------------------------------------------------
    
    /*Calculates the total force on each particle due to all other particles.
     Can be made faster by a factor of 2 if we skip the same calculation of same pairs,
     and only reverse the direction of the force. Can also be made faster by the
     Neighbour List method or a cut-off radius.*/
    void Update_Forces(){
        
        for(int i=0; i<N; i++){
            atom[i].clear_F();
            for(int j=0; j<N; j++){
                if (j == i) continue;
                atom[i].add_to_F(Force(i,j));
            }
        }
        
    }
    
    //------------------------------------------------------------------------------------
    
    /*Calculates the total force on each particle due to all other particles.
     Is based on a Neighbours List technique*/
    void Update_Forces_NL(){
        
        for(int i=0; i<N; i++){
            atom[i].clear_F();
            for(int j=0; j<atom[i].get_NL_size(); j++){
                
                atom[i].add_to_F(Force(i,atom[i].get_NL(j)));
            }
        }
        
    }
    
    
    //------------------------------------------------------------------------------------
    
    //Calculates the LJ Force between two particles.
    Vector3d Force(int i, int j){
        
        Vector3d dir = Pair_Vector(i, j);
        
        return LJ(dir.norm())*dir;
        
    }
    
    //------------------------------------------------------------------------------------
    // Physical Quantities
    //------------------------------------------------------------------------------------
    
    //Calculate the temperature from the velocities using 3/2NKT = 1/2mv^2
    double Temperature(){
        
        double sum = 0;
        for(int i=0; i<N; i++)
            sum += pow(atom[i].get_v().norm(), 2);
        
        return 1.0/(3.0*N)*sum;
        
    }
    
    //------------------------------------------------------------------------------------
    
    /*Calculate the pressure from the virial function (F_ij.r_ij) as:
     P = rho*T + 1/(3V)*virial */
    double Pressure(){
        
        double virial = 0;
        
        for(int i=0; i<N; i++){
            for(int j=i+1; j<N; j++){
                
                virial += Force(i,j).transpose()*Pair_Vector(i,j);
            }
        }
        return (rho*T + 1.0/(3.0*V)*virial);
    }
    
    //------------------------------------------------------------------------------------
    
    //Calculate the potential energy by summing the Lennard Jones potentials between all pairs.
    double U(){
        
        double sum = 0;
        double r;
        for(int i=0; i<N; i++){
            for(int j=i+1; j<N; j++){
                
                r = Pair_Vector(i, j).norm();
                sum += 4*(pow(r,-12) - pow(r, -6));
                
            }
        }
        return sum/N;
    }
    
    
    //Calculate the potential energy by summing the Lennard Jones potentials between all pairs.
    double LocalU(int i){
        
        double sum = 0;
        double r;
        for(int j=0; j<N; j++){
                if (j == i) continue;
                r = Pair_Vector(i, j).norm();
                sum += 4*(pow(r,-12) - pow(r, -6));
            
        }
        return sum/N;
    }
    
    //------------------------------------------------------------------------------------
    
    /*Calculate and write to file the distance between all pairs of particles.
     The correlation function then corresponds to the histogram of particle separations.*/
    void Correlation_Function(){
        
        unordered_map<string,double> hist;
        
        ofstream myfile ("Correlation.csv");
        
        double dist;
        
        for(int i=0; i<N; i++){
            for(int j=i+1; j<N; j++){
                
                dist = Pair_Vector(i, j).norm();
                
                if(dist <= L/2) hist[to_string(roundf(dist*100)/100)] += 1;
                
            }
        }
        
        for(double i=0.01; i<= L/2; i+= 0.01){
            myfile << i << "," << hist[to_string(i)]/(i*i) << endl;
        }
    }
    
    /*This function calculates the Mean Square Displacement of the atoms.
     It runs for several steps and calculate the average square displcaement
     of atoms from their initial positions. Since some atoms leave and enter
     the box, which will affect the mean displacement, this subroutine only
     looks at atoms which are 1/4L away from the Cell sides.*/
    void MSD(int steps){
        
        //Gnuplot gp;
        ofstream myfile2 ("MSD.csv");
        int count = 0;
        Vector3d r_0[N];
        vector<double> y;
        vector<double> x;
        
        for(int i=0; i<N; i++){
            
            //Only deal with the atoms away from the Cell sides
            if(   atom[i].get_r(0) <0.75*L && atom[i].get_r(0) > 0.25*L
               && atom[i].get_r(1) <0.75*L && atom[i].get_r(1) > 0.25*L
               && atom[i].get_r(2) <0.75*L && atom[i].get_r(2) > 0.25*L){
                r_0[i] = atom[i].get_r();   //Record the initial positions r_0
                count++;                    //Count the number of atoms
            }
            else r_0[i](0) = 0;           //All othe particles are assigned r_0(x) = 0
        }
        
        for(int t=1; t<=steps; t++){
            
            double msd = 0;
            
            Update_r_F_v();
            
            for(int i=0; i<N; i++){
                
                if (r_0[i](0) == 0) continue; //Ignore the particles close the sides (which were assigned r_0(x) = 0
                
                msd = pow((atom[i].get_r() - r_0[i]).norm(), 2); //Calculate MSD
                
            }
            
            
            myfile2 << t << "," << msd/(1.0*count) << endl; //Divide by 'count' to get the average, and write to file.
            
            if(t%20 == 0)
                Rescale_Velocities();
            
        }
        
        // gp << "set xrange [0:30]\nset yrange [0:1]\n";
        // '-' means read from stdin.  The send1d() function sends data to gnuplot's stdin.
        // gp << "plot '-' with vectors title 'time', '-' with vectors title 'msd'\n";
        // gp.send1d(x);
        // gp.send1d(y);
        
    }
    
    
    //------------------------------------------------------------------------------------
    // Print methods
    //------------------------------------------------------------------------------------
    
    //Prints physical quantities at step i.
    void print_Values(int i){
        cout << i << " " << Pressure()/rho << " " << U() << endl;
    }
    
    //------------------------------------------------------------------------------------
    
    //prints the position, velocity and acted upon total force of the particle i.
    void print_r_v_F(int i){
        
        cout << "(" << atom[i].get_r(0)
        << "," << atom[i].get_r(1)
        << "," << atom[i].get_r(2)
        << ")" << "  "
        << "(" << h*atom[i].get_v(0)
        << "," << h*atom[i].get_v(1)
        << "," << h*atom[i].get_v(2)
        << ")" << "   "
        << "(" << h*h*0.5*atom[i].get_F(0)
        << "," << h*h*0.5*atom[i].get_F(1)
        << "," << h*h*0.5*atom[i].get_F(2)
        << ")" << endl;
        
    }


    void print_v(int i){
        
        cout
        << "(" << h*atom[i].get_v(0)
        << "," << h*atom[i].get_v(1)
        << "," << h*atom[i].get_v(2)
        << ")" << endl;
        
    }
    
    
    bool HardSphere(int i, double Sphere_Diameter){
        
        bool check = true;
        
        for(int j=0; j<N; j++){
         
            if (j == i) continue;
            
            if(Pair_Vector(i, j).norm() < Sphere_Diameter){
                check = false;
                break;
            }
        }
        return  check;
    }
    
    void RandomWalk(int i, double delta){
        
        Vector3d shift;
        
        double dx = (((double) rand() / (RAND_MAX))-0.5)/sqrt(3)*delta;
        double dy = (((double) rand() / (RAND_MAX))-0.5)/sqrt(3)*delta;
        double dz = (((double) rand() / (RAND_MAX))-0.5)/sqrt(3)*delta;
        
        shift << dx,dy,dz;
        atom[i].shift_r(shift);
        
        
    }
    
    
    
    
};


//************************************************************************************
//Methods
//************************************************************************************



/*Performs a Wang-Landau sampling. Starts with a random configuration, and performs a random
 walk based on a rule. Keeps record of number of times an energy state is visited.*/
void WLSampling(double lattice, int box, double delta, double E_0, double E_N, double bin_width, int INNER_STEPS, int ITERS){
   
    ofstream myfile ("Histogram.csv");
    
    int N = 4*pow(box,3);
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
    for(int i = trunc(E_0/bin_width); i<=trunc(E_N/bin_width); i++)
        lnOmega[i] = 0;
    
   
    Cell cell1(lattice, box, 3.0, 0.004);   //Create an initial lattice
    Cell cell2(lattice, box, 3.0, 0.004);   //Create a proposed lattice
    cell1.Initialise();
    cell2.Initialise();
    
    
    E1 = cell1.U()/N;
    binE1 = trunc(E1/bin_width);
    
    int iter = 0;
    while(iter < ITERS){
        
        for(int i = trunc(E_0/bin_width); i<=trunc(E_N/bin_width); i++)
            hist[i] = 0;        //Keep a histogram of all the visited sites.
        
        step = 1;
        while(step < INNER_STEPS){
            
            int a = rand() % N; //Pick an atom in random
            
            cell2 = cell1;
            cell2.RandomWalk(a, delta);     //Make a random walk from cell1.
            E2 = E1 - cell1.LocalU(a) + cell2.LocalU(a);
            
            while(!cell2.HardSphere(a, 0.9)){
                cell2 = cell1;
                cell2.RandomWalk(a, delta);
             //   cout << "-----------------------------" << cell2.HardSphere(a, 0.9) << endl;
            }
            
            E2 = E1 - cell1.LocalU(a) + cell2.LocalU(a);
            binE2 = trunc(E2/bin_width);
            
           // cout << f << " " << iter << " " << step << " " << E1 << " " << E2 << " " << binE1 << " " << binE2 << endl;
            
            //Make a transition to the new state according to p(E1 -> E2) = min (1, Omega(E1) / Omega(E2))
            //Since we're dealing with log(Omega), we subtract the two values.
            if (( log((double) rand() / (RAND_MAX))) <= (double) (lnOmega[binE1] - lnOmega[binE2]) ){
                cell1 = cell2;
                E1 = E2;
                binE1 = binE2;
            }
            
            lnOmega[binE1] += log(f);       //Update Omega(E) by f*Omega(E). In log, we need to add.
            hist[binE1] += 1;               //Add one to the histogram of visited sites.
            
            step++;
            if (step % 10000 == 0) cout<< step << endl;
        }
        
        //Once the histogram is flat enough, update f, and repeat everything again (histogram is reset as well).
        //Make consequent iterations slightly longer by increasing the INNER_STEPS.
        f = sqrt(f);
        iter++;
        INNER_STEPS += 20000;
        
    }
    
    for(int i = trunc(E_0/bin_width); i<=trunc(E_N/bin_width); i++)
        myfile << i << "," <<  lnOmega[i]  << "," << hist[i] << endl;
    
    
}


void SpecificHeat(double E_0, double E_N, double bin_width, double T_0, double T_N, double T_bin_width){

    ofstream myfile ("HeatCapacity.csv");
    vector<string> data = getData("Histogram2.csv");
    
    double Z = 0;
    double ExpE = 0;
    double ExpE2 = 0;
    double C_v;
    double lnOmega[data.size()];
    
    int j = 0;
    for (int i = 1; i <= data.size(); i+=2){
        lnOmega[j] =  atof(data[i].c_str());
        j++;
    }
    
 //   for(int T = trunc(T_0/T_bin_width); T<=trunc(T_N/T_bin_width); T++){
    double T = 1.1;
        for(int E = trunc(E_0/bin_width); E<=trunc(E_N/bin_width); E++){
        
            Z += lnOmega[E]*exp(-E/T);
            ExpE += E*lnOmega[E]*exp(-E/T);
            ExpE2 += E*E*lnOmega[E]*exp(-E/T);
            cout << ExpE << endl;
        
        }
        C_v = (ExpE2/Z - ExpE*ExpE/(Z*Z))/(T*T);
        myfile << T << "," <<  Z << "," << ExpE << "," << ExpE2 << "," << C_v << endl;
    }
    
    
    
//}









//************************************************************************************
//Main
//************************************************************************************

int main(int argc, const char * argv[]) {
    
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    
  //  Cell cell(2.37, 3, 3.0, 0.004); //gas
   // cell.Initialise();
    

    

 //   WLSampling(2.37, 3, 0.5, -5, 20, 0.05, 200000, 3);

    SpecificHeat(-5, 20, 0.05, 0.1, 1, 0.1);

    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout <<  cpu_time_used << endl;
    
    return 0;
}
