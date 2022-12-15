#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <time.h>

using namespace std;

const double PI  =3.141592653589793238463;

//********************************************************************************
//Classes
//********************************************************************************

/*The Class Atom contains all the information about the atom we are calculating.
 Specifically, these are the atomic number, ionization, and information specific
 to each atom such as the symbol and the orbital energy root-finding domain.
 Note that the latter is not needed if we are using plain-wave bases, as the
 eigenvalues can be found through the diagonalisation of the matrix.*/

class Atom{
    
private:
    
    int Z;          //The Atomic number
    int n_elec;     //The number of electrons
    
public:
 
    Atom(int a){
        Z = a;
        n_elec = a; //By default, we assume the atom is not ionized
        
    }
    
    Atom(int a, int ionization){
        Z = a;
        n_elec = a-ionization;
        
    }
    
    //Returns the symbol of the element from its atomic number
    string Atomic_Symbol(){
        
        string symbol;
        
        switch(Z){
            case 1 : symbol = "H";
                break;
            case 2 : symbol = "He";
                break;
            case 3 : symbol = "Li";
                break;
            case 4 : symbol = "Be";
                break;
            case 5 : symbol = "B";
                break;
            case 6 : symbol = "C";
                break;
            case 7 : symbol = "N";
                break;
            case 8 : symbol = "O";
                break;
            case 9 : symbol = "F";
                break;
            case 10 : symbol = "Ne";
                break;
        }
       
        return symbol;
    }
    
    /*This function returns the search domain of epsilon in finding the root.
    The domain of search is dependant on the system.*/
    double* epsilon_search_domains(){
        
        int n_orbitals = 1;
        
        switch(n_elec){
            case 1 ... 2 : n_orbitals = 1;
                break;
            case 3 ... 4 : n_orbitals = 2;
                break;
            case 5 ... 10 : n_orbitals = 3;
                break;
            case 11 ... 12 : n_orbitals = 4;
                break;
            case 13 ... 18 : n_orbitals = 5;
                break;
        }
        
        double* search_domains = new double[n_orbitals+1];
        search_domains[n_orbitals] = 0;
        
        switch(Z){
            case 1 ... 2 : search_domains[0] = -3;
                break;
            case 3 ... 4 : search_domains[0] = -5;
                           search_domains[1] = -2;
                break;
            case 5 ... 6 : search_domains[0] = -13;
                           search_domains[1] = -4;
                           search_domains[2] = -0.2;
                break;
            case 7 ... 8 : search_domains[0] = -24;
                           search_domains[1] = -12;
                           search_domains[2] = -1;
                break;
            case 9 ... 10 : search_domains[0] = -38;
                            search_domains[1] = -25;
                            search_domains[2] = -1;
                break;
        }
        
        return search_domains;
        
    }
    
    /*Returns the angular momentum depending on which electron
     we're looking at*/
    int l(int a){
        
        int l;
        a++; //Since we begin counting electrons from 0, add 1.
        switch(a){
            case 1 ... 4 : l = 0;
                break;
            case 5 ... 10 : l = 1;
                break;
            case 11 ... 12 : l = 0;
                break;
            case 13 ... 18 : l = 1;
                break;
        }
        
        return l;
        
    }
    
    int get_Z(){
        return Z;
    }
    
    int get_n_elec(){
        return n_elec;
    }
    
};

//********************************************************************************

/* The Mesh class contains information about the space. For a spherically symmetric
 system such as single atoms, we are only concerned with the radial dimension. This
 class defines the fineness of the mesh and the limit of integration. */

class Mesh
{
    
private:
    
    double h;   //The steps in the mesh
    int r_max;  //Maximum distance of integration from the origin
    int n_max;  //Maximum number of steps
    
    
public:
    
    Mesh(){
        
        h = 0.01;
        r_max = 50;
        n_max = r_max/h;
    
    }
    
    double get_h(){
        return h;
    }
    
    int get_n_max(){
        return n_max;
    }
    
    double get_r_max(){
        return r_max;
    }
    
};

//********************************************************************************

/*This is the main class which deals with the wavefunction in hand. Most importantly,
 it has functions for calculating the single-orbital wavefunctions and normalising them,
 the eigenvalues of the Kohn-Sham equation, electron density, Hartree energy and Exchange energy.
 It also has a Self-Consistent Field method for accurately solving the Kohn-Sham equation.*/

class Wavefunction{
    
private:
    
    double* epsilon;    //The array which contains all the eigenvalues of KS equation. Dimension = n_elec
    double* n;          //The vector of charge density n(r). Dimension = n_max
    
    double* V_H;        //Hartree Energy V_H(n(r), r)
    double* V_x;        //Exchange Energy V_x(n(r), r)
    
    double E;           //Total Energy

    
public:
    
    Wavefunction(Atom atom, Mesh mesh){
        
        double h = mesh.get_h();
        double n_max = mesh.get_n_max();
        int n_elec = atom.get_n_elec();
        
        E = 0;

        n = new double [n_max+1];
    
        V_H = new double [n_max+1];
        V_x = new double [n_max+1];
        
        epsilon = new double [n_elec];
        
        //Begin with a guess for n(r), V_H and V_x. We set them all to 0.
        for(int i = 0; i<=n_max; i++){
            n[i] = 0;
            V_H[i] = 0;
            V_x[i] = 0;
        }
        
    }
  
//------------------------------------------------------------------------------------
 
    //This function normalises the single-orbital wavefunction u(r)
    void Normalise(double* &u, Mesh mesh){
        
        int n_max = mesh.get_n_max();
        double N;
        double* f = new double [n_max+1];
        
        for(int i =0; i <= n_max; i++)
            f[i] = u[i]*u[i];
        
        //Integrate u^2(r) to get N. Then divide u(r) by Sqrt[N] to normalise it.
        double sum = (f[0]+f[n_max])*0.5;
        for(int i =1; i < n_max; i++)
            sum+= f[i];
        
        N = sum*mesh.get_h();        //The normalisation factor
        
        for(int i =0; i <= n_max; i++)
            u[i] = u[i]/sqrt(N);     //Divide by sqrt(N) to normalise
        
    }

//------------------------------------------------------------------------------------
    
    /*This function calculates the single-orbital wavefunctions u(r)
    using Numerov's integration algorithm*/
    double* Calculate_u(Atom atom, double epsilon, int l, Mesh mesh){
        
        int n_max = mesh.get_n_max();
        double r_max = mesh.get_r_max();
        double h = mesh.get_h();
        int Z = atom.get_Z();
        
        double* u = new double [n_max+1];
        double* w = new double [n_max+1];
        double* f = new double [n_max+1];
        
        //Define the function f(r) in the Numerov algorithm
        for(int i = n_max; i>0; i--){
            double r = i*h;
            f[i] = 2*(-Z/r -epsilon + l*(l+1)/(2*r*r) + V_H[i] + V_x[i]);
        }
        
        //Boundary Conditions
        u[n_max] = r_max*exp(-sqrt(-2*epsilon)*r_max);
        u[n_max-1] = (r_max-h)*exp(-sqrt(-2*epsilon)*(r_max-h));
        w[n_max] = (1-h*h*f[n_max]/12)*u[n_max];
        w[n_max-1] = (1-h*h*f[n_max-1]/12)*u[n_max-1];
        
        //Perform Numerov integration
        for(int i = n_max-1; i>0; i--){
            w[i-1] = 2*w[i]-w[i+1]+h*h*f[i]*u[i];
            u[i-1] = w[i-1]/(1-h*h*f[i-1]/12);
        }
            
        //At r=0, f(0) diverges, and in above, f[0] = 0 which is wrong. So manually set u[0] = w[0]
        u[0] = w[0];
    
        Normalise(u, mesh);     //Normalise u(r)
        
        return u;
    }

    
//------------------------------------------------------------------------------------
 
    /*This is a root-finding Bisection algorithm in order to find the
    eigenvalues of the KS equation, epsilon_i. We used the boundary condition
     at r_max to find u(r), now we use boundary condition at u(0) to find
     epsilon_i. Since u(0) = 0, we use this info to find all the roots which
     u(r=0, epsilon) = 0.*/
    double Find_epsilon(Atom atom, double a, double b, int l, Mesh mesh){
        
        double deltas = 0.000000001;
        double epsn = 0.00000001;
        double maxit = 100;

        double u_0 = Calculate_u(atom, a, l, mesh)[0];
        double e = b-a;
        double c;;
        
        for(int k=1; k <= maxit; k++){
            
            e*= 0.5;
            c = a + e;
            double w = Calculate_u(atom, c, l, mesh)[0];

            if (fabs(e) < deltas || fabs(w) < epsn) break;    //terminate if e is too small, or w has converged to 0.
            ((u_0>0 && w<0) || (u_0<0 && w>0)) ? (b=c) : (a=c , u_0=w); //if u*w < 0, then move b to c, otherwise move a.

        }

        return c;

    }
    
//------------------------------------------------------------------------------------
    
    /*This function is for atoms which have more than one electron orbital, basically
     anything larger than Helium. Thus we need to calculate the epsilon_i for each of
     these orbitals. Note that the electrons in the same n and l orbital are degenerate.
     We have manually entered the number of electrons and the angular momentum l for each
     orbital. This is not ideal, but gets the job done for single atoms.*/
    void Find_all_epsilons(Atom atom, Mesh mesh){
        
        int n_elec = atom.get_n_elec();
        double* search_domain = atom.epsilon_search_domains();
        
        int n_elec_in_orbital[5] = {2,2,6,2,6};
        int i =0;
        int k = 0;
        
        while(i < n_elec){  //i goes through all the electrons
            
            int j = 0;
            
             /*j measure the number of electrons in an orbital. The loop terminates once
             the epsilon of all the electrons in the orbital is calculated. If the orbital
             is not full, the loop terminates if i reaches maximum number of electrons.
             Note that since the electrons in the same orbital are degenerate, we only need
             to calculate the epsilon once, and the rest would be the same.
             k is a measure of the orbital. Every time the loop is terminated, k goes
             to the next orbital. Note that the domain of search and l depend on the orbital.*/
            while(j < n_elec_in_orbital[k] && i < n_elec){
                
                if(j == 0) epsilon[i] = Find_epsilon(atom, search_domain[k], search_domain[k+1], atom.l(i), mesh);
                else epsilon[i] = epsilon[i-1];
                j++;
                i++;
            }
            k++;
        }
    }

//------------------------------------------------------------------------------------
    
    /*Calculates the charge density from u_i(r) from
     n(r) = sum_i u_i^2/r^2 */
    void Calculate_n(Atom atom, Mesh mesh){
        
        
        int n_max = mesh.get_n_max();
        double h = mesh.get_h();
        int n_elec = atom.get_n_elec();
        double* u;
        
        //reset n(r) to 0 every time a new calculation of n(r) is done.
        for(int j = 0; j<=n_max; j++){
            n[j] = 0;
        }
        
        for(int i=0; i < n_elec; i++){
            
            u = Calculate_u(atom, epsilon[i], atom.l(i), mesh);
            
            for(int j = 1; j<=n_max; j++)
                 n[j] += u[j]*u[j]/(j*h*j*h);

        }
        
    }

//------------------------------------------------------------------------------------
    
    /*Method to calculate V_H(r) from n(r). We use Verlet algorithm
     and the boundary conditions on U(r) = rV_H(r). */
    void Calculate_Hartree_Potential(Atom atom, Mesh mesh){
        
        if(atom.get_Z() == 1) return;
        
        int n_max = mesh.get_n_max();
        double h = mesh.get_h();
        double* U = new double [n_max+1];
        
        //Boundary Conditions
        U[0] = 0;
        U[1] = h;
        
        //Verlet Algorithm to integrate U
        for(int i = 1; i<n_max; i++){
            double r = i*h;
            U[i+1] = 2*U[i] - U[i-1] - h*h*n[i]*r;
        }
        
        //The integration constant is alpha*r. Fix alpha such that U(r_max) = q_max/r_max
        double alpha = (atom.get_n_elec() - U[n_max-1])/(n_max*h);
        
        /*Add the integration constant to U, and calculate V_H as V_H = U/r.*/
        for(int i = 1; i<=n_max; i++){
            double r = i*h;
            U[i] += alpha*r;
            V_H[i] = U[i]/r;
        }
        V_H[0] = V_H[1];

    }

//------------------------------------------------------------------------------------
    
    /*Calculates exchange potenital using Local Density Approximation (LDA) */
    void Calculate_Exchange_Potential(Atom atom, Mesh mesh){
        
        if(atom.get_Z() == 1) return;
        
        int n_max = mesh.get_n_max();
        
        //Calculate V_x. Skip 0 to avoid division by zero.
        for(int i =1; i <= n_max; i++)
            V_x[i] = -pow(3*n[i]/(4*PI*PI), 1.0/3);
        
        //Set V_x(0) = V_x(h)
        V_x[0] = V_x[1];

    }
    
//------------------------------------------------------------------------------------
 
    /*Integrate V_H(r)n(r) to get the Hartree energy */
    double Hartree_Energy(Mesh mesh){
        
        int n_max = mesh.get_n_max();
        double h = mesh.get_h();
        
        //Use Trapezoidal algorithm to integrate
        double* f = new double [n_max+1];
        for(int i =0; i <= n_max; i++){
            double r = i*h;
            f[i] = V_H[i]*n[i]*r*r;
        }
        

        double sum = (f[0]+f[n_max])*0.5;
        for(int i =1; i < n_max; i++)
            sum+= f[i];
        
        return sum*mesh.get_h();
    }
    
//------------------------------------------------------------------------------------
    
    /*Integrate V_x(r)n(r) to get the Exchange energy */
    double Exchange_Energy(Mesh mesh){
        
        int n_max = mesh.get_n_max();
        double h = mesh.get_h();
        
        //Use Trapezoidal algorithm to integrate
        double* f = new double [n_max+1];
        for(int i =0; i <= n_max; i++){
            double r = i*h;
            f[i] = V_x[i]*n[i]*r*r;
        }
        
        double sum = (f[0]+f[n_max])*0.5;
        for(int i =1; i < n_max; i++)
            sum+= f[i];
        
        return sum*mesh.get_h();
        
    }

//------------------------------------------------------------------------------------
    
    /*Calculate the total energy by using:
    E = sum_i epsilon_i - 0.5*int[V_H(r)n(r)] + E_xc - int[V_x(r)n(r)]*/
    void Calculate_E(Atom atom, Mesh mesh){
        
        E = 0;
        
        for (int i =0; i<atom.get_n_elec(); i++)
            E += epsilon[i];
        
        E +=  -0.5*Hartree_Energy(mesh) - 0.25*Exchange_Energy(mesh);
        
    }
    
//------------------------------------------------------------------------------------
    
    /*Prints the total energy, hartree energy and exchange energy.
     Gives the option to print the sum of epsilon_i, or print each
     epsilon_i separately. Can use the arguments 'total' or 'each'
     To control this. */
    void Print_Energy(Atom atom, Mesh mesh, string which_eps){
    
    cout << atom.Atomic_Symbol() << "  " << E << "         " << Hartree_Energy(mesh) << "         " << Exchange_Energy(mesh) << "      ";
    
    if(which_eps == "total"){
        
        double sum = 0;
        
        for (int i =0; i<atom.get_n_elec(); i++)
            sum += epsilon[i];
        
        cout << sum << endl;
    }
        
    if(which_eps == "each"){
        
        for (int i =0; i<atom.get_n_elec(); i++)
            cout << epsilon[i] << "    ";
    
        cout << endl;
    }

}

//------------------------------------------------------------------------------------
    
    /*This is the most important method of the class which glues all the parts together
     and solves the Kohn-Sham equation self-consistently until the energy is converged.
     It begins finding the eigenvalues (epsilons), calculating the eigenfunctions u_i(r)
     from those eigenvalues and thus calculating the charge density n(r), calculating
     Hartree energy potential and exchange potential from n(r), and finally calculating
     the total electronic energy. The process iterates, where the calculated V_H(r) and
     V_x(r) are fed back into Kohn-Sham equation for the new calculation of epsilon_i
     and u_i(r) (remember that V_H(r) and V_x(r) are private members of the class
     Wavefunction. The iteration stops when the change in total energy is negligible.*/
    void SCF(Atom atom, Mesh mesh){
        
        double E_Precision = 0.00001;
        double delta = 1;
        int max_iter = 100;
        int iter = 0;
        double dummy_E;
        
        while(abs(delta) > E_Precision && iter < max_iter){
        
            Find_all_epsilons(atom, mesh);
            Calculate_n(atom, mesh);
            Calculate_Hartree_Potential(atom, mesh);
            Calculate_Exchange_Potential(atom, mesh);
            
            dummy_E = E;
            Calculate_E(atom, mesh);
            delta = E - dummy_E;

            //Print_Energy(atom, mesh, "total");
            iter++;
        }
        
    }
    
};


//************************************************************************************
//Methods
//************************************************************************************

void Print_Header(Mesh mesh){
    
    cout << "Atom  " << "Total Energy  " << "Hartree Energy  " << "Exchange Energy  " << "Total epsilon" << endl;
    
}

//------------------------------------------------------------------------------------

void Print_Header(Atom atom, Mesh mesh, string which_eps){
    
    cout << "Atom  " << "Total Energy  " << "Hartree Energy  " << "Exchange Energy  ";
    
    if(which_eps == "total")
        cout << "Total epsilon" << endl;
    
    if(which_eps == "each"){
        
        for (int i =0; i<atom.get_n_elec(); i++)
            printf("epsilon %d    ", i);
        cout << endl;
    }
}

//------------------------------------------------------------------------------------

/*The method calculates and prints the electronic energies of single atoms.
 The choice of 'up to' calculate the energies of atoms from Hydrogen up to
 the given atomic number. The 'only' argument only calculates the energy of
 the atom with the given atomic number.*/
void Calculate_Atomic_Energies(string which_atoms, int a){
    
    int from;
    if(which_atoms == "up to") from = 1;
    if(which_atoms == "only")  from = a;
        
    Mesh mesh;                              //Create a mesh of the radial coordinate.
    Print_Header(mesh);
       
    for(int i = from; i <=a; i++){
 
        Atom The_Atom(i);                   //Create an Atom object with atomic number i
        Wavefunction Psi(The_Atom, mesh);   //Create a wavefunction Psi for The Atom, defined over the mesh.
        Psi.SCF(The_Atom, mesh);            //Solve the problem self-consistently to obtain the energy.
        Psi.Print_Energy(The_Atom, mesh, "total");
    }

}

//************************************************************************************
//Main
//************************************************************************************

int main(int argc, const char * argv[]) {
    
    Calculate_Atomic_Energies("up to", 10);

    return 0;
}
