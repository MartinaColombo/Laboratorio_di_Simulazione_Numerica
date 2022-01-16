//Esercizio 3 
//Andamento del prezzo di un opzione europea tramite GBM

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include "random.h"


using namespace std;
 
int main (){

   Random rnd;
   int seed[4];
   int p1, p2;
     
   int M=pow(10,4); //number of MC steps
   int N=pow(10,2); //numbers of blocks
   int L=M/N; //lenght of a block
   int time_intervals = pow(10,2); 
   int S_0=100;  //asset price in t=0
   int T=1; //delivery time
   int K=100;  //strike price
   double t_step=T*1./K; //lenght of one subinterval 
   double r=0.1; //risk-free interest rate
   double vol=.25; //volatility
   
   double S=0, C=0, P=0;
   double sum_progC=0, sum_progP=0;
   double mean_progC=0, mean_progP=0;
   double mean_prog2C=0, mean_prog2P=0;
   double mean_C=0, mean_P=0;
   double mean_2C=0, mean_2P=0;
   double error_progC=0, error_progP=0;
   double W=0, Z=0;
   
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


//Esercizio 3.1

   output.open("direct_samp.dat");
   for(int i=1; i<=N; i++){
      sum_progC=0;
      sum_progP=0;
      for(int k=0; k<L; k++){
         W=rnd.Gauss(0,T);
         S=S_0*exp((r-pow(vol,2)/2)*T+vol*W);  
         C=exp(-r*T)*max(0.0,S-K); //Call
         P=exp(-r*T)*max(0.0,K-S); //Put
         sum_progC+=C;
         sum_progP+=P;
      } 
      //Blocking
      mean_C+= sum_progC/L; //A
      mean_2C+=pow(sum_progC/L,2); //A^2
      mean_progC=mean_C/i; //<A>  
      mean_prog2C=mean_2C/i;  //<A^2>
     
      mean_P+= sum_progP/L;
      mean_2P+= pow(sum_progP/L,2);
      mean_progP=mean_P/i; 
      mean_prog2P=mean_2P/i; 
     
      if(i-1==0){
         error_progC=0;
         error_progP=0;    
      }
      else{
         error_progC=sqrt((mean_prog2C-pow(mean_progC,2))/(i-1));
         error_progP=sqrt((mean_prog2P-pow(mean_progP,2))/(i-1));       
      }       
      output<<i<<" "<<mean_progC<<" "<<error_progC<<" "<<mean_progP<<" "<<error_progP<<endl;
   } 
   output.close();


//Esercizio 3.2

   S=0;
   C=0;
   P=0;
   sum_progC=0;
   sum_progP=0;
   mean_progC=0;
   mean_progP=0;
   mean_prog2C=0;
   mean_prog2P=0;
   mean_C=0;
   mean_P=0;
   mean_2C=0;
   mean_2P=0;
   error_progC=0;
   error_progP=0;
   

   output.open("disc_samp.dat");
   for(int i=1; i<=N; i++){
      sum_progC=0;
      sum_progP=0;
      for(int k=0; k<L; k++){
         S=S_0;
         for(int j=1; j<time_intervals; j++){
           Z=rnd.Gauss(0,1);
           S=S*exp((r-pow(vol,2)/2)*t_step+vol*Z*sqrt(t_step)); 
           }     
          
         C=exp(-r*T)*max(0.0,S-K); //Call
         P=exp(-r*T)*max(0.0,K-S); //Put
         sum_progC+=C;
         sum_progP+=P;
      } 
      //Blocking
      mean_C+= sum_progC/L; //A
      mean_2C+=pow(sum_progC/L,2); //A^2
      mean_progC=mean_C/i; //<A>  
      mean_prog2C=mean_2C/i;  //<A^2>
      
      mean_P+= sum_progP/L;
      mean_2P+= pow(sum_progP/L,2);
      mean_progP=mean_P/i; 
      mean_prog2P=mean_2P/i; 
     
      if(i-1==0){
         error_progC=0;
         error_progP=0;    
      }
      else{
         error_progC=sqrt((mean_prog2C-pow(mean_progC,2))/(i-1));
         error_progP=sqrt((mean_prog2P-pow(mean_progP,2))/(i-1));       
      }       
      output<<i<<" "<<mean_progC<<" "<<error_progC<<" "<<mean_progP<<" "<<error_progP<<endl;
   } 
   output.close();

   rnd.SaveSeed();
   return 0;
}

