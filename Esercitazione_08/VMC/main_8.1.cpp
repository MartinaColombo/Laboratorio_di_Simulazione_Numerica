#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;


int main(){

   Random rnd;
   int seed[4];
   int p1, p2;
   int hit=0, hit_m=0, hit_s=0;
   int M=pow(10,4); //N throws
   int N_blocks=pow(10,2); //N blocks
   
   double x=0.5, x_new=0; 
   double delta=1.5, beta=1./1.25;
   double psi_new=0, psi_old=0, psi_trial=0, deriv_psi=0, potential=0; 
   double accept=0; 
   double rand=0, sum_H=0, H=0, H_new=0;
   double mean_H=0, mean2_H=0;
   double mean_prog_H=0, mean2_prog_H=0;
   double mu=0.1, sigma=1;
   //double h=1., mass=1.;
   
   ofstream output;
 /**********************************************************************************/    
   ifstream Primes("Primes");
   
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   
   ifstream input("seed.in");
   string property;
   
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  /**********************************************************************************/  
   output.open("Hamiltonian.dat");
   for(int i=0; i<N_blocks; i++){
     H=0; 
     for(int j=1; j<=M/N_blocks; j++){       
       psi_old=pow(exp(-pow(x-mu,2)/(2*pow(sigma,2))) + exp(-pow(x+mu,2)/(2*pow(sigma,2))),2);
       x_new=rnd.Rannyu(x-delta, x+delta);       
       psi_new=pow(exp(-pow(x_new-mu,2)/(2*pow(sigma,2))) + exp(-pow(x_new+mu,2)/(2*pow(sigma,2))),2);    
       accept=min(1.,psi_new/psi_old);                    
       rand=rnd.Rannyu();     
       if(rand<=accept){
         x=x_new;
         hit++;               
       }
       deriv_psi=1./pow(sigma,2) *(exp(-pow(x-mu,2)/(2*pow(sigma,2)))*(pow(x-mu,2)/pow(sigma,2)-1) + exp(-pow(x+mu,2)/(2*pow(sigma,2)))*(pow(x+mu,2)/pow(sigma,2) -1));
       psi_trial=exp(-pow(x-mu,2)/(2*pow(sigma,2)))+exp(-pow(x+mu,2)/(2*pow(sigma,2)));  
       potential=pow(x,4)-5*1./2 * pow(x,2);
      // potential=0.5*pow(x,2);   
       H+= -1.*0.5*deriv_psi/psi_trial + potential;     
       }     
     //Blocking    
     mean_H+=H*1./(M/N_blocks);
     mean2_H+=pow(H*1./(M/N_blocks),2);
     
     mean_prog_H=mean_H*1./(i+1);
     mean2_prog_H=mean2_H*1./(i+1);
     if(i==0){
       output<<i<<" "<<mean_prog_H<<" "<<0.0<<endl;    
     }
     else{
       output<<i<<" "<<mean_prog_H<<" "<<sqrt((mean2_prog_H-mean_prog_H*mean_prog_H)/i)<<endl;   
     }       
   }
   output.close();     
   cout<<"Percentage of acceptance x: "<<hit*1.*100/M<<endl;

return 0;
} 

