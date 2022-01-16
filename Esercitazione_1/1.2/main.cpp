#include <iostream>
#include <fstream>
#include <string>
#include "random.h"


//Esercizio 1.3


using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   int v[4]={1,2,10,100}; //N
   int throws=pow(10,4);
   double r=0; //random number
   double exp=0, lor=0;
   //parameters
   int lambda=1;
   int gamma=1;

   double sum_r=0, sum_exp=0, sum_lor=0;
   ofstream out_flat;
   ofstream out_exp;
   ofstream out_lorentz;
   
  
   
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

   
   
   
   out_flat.open("flat.dat");
   out_exp.open("exp.dat");
   out_lorentz.open("lorentz.dat");
  
   for(int j=0; j<throws;j++){ 
      
       r = rnd.Rannyu();
       exp=rnd.Exp(lambda);    
       lor=rnd.Lorentz(gamma); 
       
       sum_r=r;
       sum_exp=exp;
       sum_lor=lor; 
       
       out_flat<<r<<" ";
       out_exp<<sum_exp<<" ";
       out_lorentz<<sum_lor<<" ";
       
       sum_r=0;
       sum_exp=0;
       sum_lor=0; 
       
       for(int i=1; i<4; i++){
         for(int z=0; z<v[i]; z++){
         r = rnd.Rannyu();
         exp=rnd.Exp(lambda);    
         lor=rnd.Lorentz(gamma);
       
         sum_r+=r;
         sum_exp+=exp;
         sum_lor+=lor;  
        
         } 
         
         out_flat<<sum_r/v[i]<<" ";
         out_exp<<sum_exp/v[i]<<" ";
         out_lorentz<<sum_lor/v[i]<<" ";
         
         sum_r=0;
         sum_exp=0;
         sum_lor=0;                               
       } 
         out_flat<<endl;
         out_exp<<endl;
         out_lorentz<<endl;      
  
   }
   out_flat.close();
   out_exp.close();
   out_lorentz.close();   

   rnd.SaveSeed();
   return 0;
}

