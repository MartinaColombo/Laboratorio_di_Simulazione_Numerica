#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   int M=pow(10,5); //Monte Carlo steps
   int N=pow(10,2); //num blocks
   int L=M/N; //blocks lenght 
   int N_hit=0; //times  
   
   double x=0, theta=0; 
   double w=0, z=0;
   double len=3.5, d=4; //needle lenght and distance between lines
   double sum_pi=0, sum2_pi=0;
   double mean_pi=0;
   double mean2_pi=0;
   double error_prog_pi=0;
   
   ofstream output;   
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
   
   
   output.open("pi.dat");
      
   for(int i=1;i<=N; i++){
      N_hit=0;
      for(int k=0; k<L; k++){
         x=rnd.Rannyu(0, d*1./2);
         do{
            w=rnd.Rannyu(0, 1);
            z=rnd.Rannyu(0, 1);          
          }while(w*w+z*z>1);         
          theta=acos(w/sqrt(w*w+z*z));          
          if(x<=((len*1./2)*sin(theta))){ //hit condition
            N_hit++;
          }                  
      }
      sum_pi+=(2*len*1./d)*(L*1./N_hit);
      sum2_pi+=pow((2*len*1./d)*(L*1./N_hit),2);
   
      mean_pi=sum_pi/(i); //<pi>
      mean2_pi=sum2_pi/(i);//<pi2>
      if(i-1==0){
         error_prog_pi=0;
      }
      else{
         error_prog_pi=sqrt((mean2_pi-mean_pi*mean_pi)/(i));
      }     
      output<<i<<" "<<mean_pi<<" "<<error_prog_pi<<endl;     
   }  

   output.close();


   rnd.SaveSeed();
   return 0;
}
