//Esercizio 2.1 
//Risoluzione di un integrale tramite metodi MC con Sampling uniforme e Importance sampling


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   int M=pow(10,4); //MC steps
   int N=pow(10,2);  //Num blocks
   int L=M/N; //lenght of a block
   
   double c=M_PI/2; 
   double x=0, y=0; 
   double sum_int=0, sum_g_int=0;
   double sum_prog_int=0, sum_prog2_int=0;   
   double sum_prog_g_int=0, sum_prog2_g_int=0;
   double error_prog_int=0, error_prog_g_int=0;
   double mean_int=0, mean2_int=0;
   double mean_g_int=0, mean2_g_int=0;
      
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


   output.open("file.dat");
      
   for(int i=1;i<=N; i++){
      sum_int=0;
      sum_g_int=0;
        for(int k=0; k<L; k++){
           x=rnd.Rannyu();	//uniform sampling in [0,1]
           y=rnd.Dist();	//importance sampling with p(x)=2(1-x) 
           sum_int+=cos(c*x);  
           sum_g_int+=cos(c*y)/(1-y);                   
        }
        //Data blocking
        sum_prog_int+=sum_int*c/L;
        sum_prog2_int+=pow(sum_int*c/L,2);     
        mean_int=(sum_prog_int)/(i); 
        mean2_int=(sum_prog2_int)/(i);
        
        sum_prog_g_int+=sum_g_int*c/(2*L);
        sum_prog2_g_int+=pow(sum_g_int*c/(2*L),2);     
        mean_g_int=(sum_prog_g_int)/(i); 
        mean2_g_int=(sum_prog2_g_int)/(i);        
   
        if(i-1==0){
           error_prog_int=0;
           error_prog_g_int=0;
        }
        else{
        error_prog_int=sqrt((mean2_int-mean_int*mean_int)/(i-1));
        error_prog_g_int=sqrt((mean2_g_int-mean_g_int*mean_g_int)/(i-1));
        }     
        output<<i<<" "<<mean_int<<" "<<error_prog_int<<" "<<mean_g_int<<" "<<error_prog_g_int<<endl;     
   }  

   output.close();

   rnd.SaveSeed();
   return 0;
}

