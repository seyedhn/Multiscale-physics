#include "Lattice.h"

//************************************************************************************
//Main
//************************************************************************************


int main(int argc, const char * argv[]) {

    
    //std::cout << haha(5);
    
    
    
    
    
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    

    SpecificHeat(-2, 2, 0.05, 0.001, 0.05, 0.001);
   // WLSampling(16, -2, 2, 0.05, 200000, 20);
    

    
    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout <<  cpu_time_used << endl;
    
    return 0;
}
