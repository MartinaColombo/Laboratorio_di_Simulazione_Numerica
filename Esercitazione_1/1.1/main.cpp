
//Esercizio 1.1 
//Stima del valor medio di r, della sua deviazione standard e del chi quadro

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "random.h"

using namespace std;
 
int main (){

   Random rnd;
   int seed[4];
   int p1, p2;
   int M=pow(10,4); //number of throws
   int N=pow(10,2); //number of blocks
   int L=M/N; //number of throws in each block

   double r=0; //random variable
   double sum_rm=0, sum_rs=0; //sums in each block 
   double sum_prog_rm=0, sum_prog_rs=0; //progressive sums
   double sum2_prog_rm=0, sum2_prog_rs=0; //squared progressive sums  
   double mean_rm=0, mean_rs=0; //mean value of a block
   double mean2_rm=0, mean2_rs=0; //squared mean value of a block
   double error_prog_rm=0, error_prog_rs=0; //standard deviation
          
   
   ifstream Primes("Primes");
   ofstream output;
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
   
//Mean and standard deviation
   
   output.open("file.dat"); 
   
   for(int i=1;i<=N; i++){
     sum_rm=0;
     sum_rs=0;
     for(int k=0; k<L; k++){
        r=rnd.Rannyu();
        sum_rm+=r;
        sum_rs+=(r-0.5)*(r-0.5);
        }
	//Blocking
        sum_prog_rm+=(sum_rm/L);
        sum2_prog_rm+=(sum_rm/L)*(sum_rm/L);
     
        sum_prog_rs+=(sum_rs/L);
        sum2_prog_rs+=(sum_rs/L)*(sum_rs/L);     
      
        mean_rm=sum_prog_rm/(i);
        mean2_rm=sum2_prog_rm/(i);
     
        mean_rs=sum_prog_rs/(i);
        mean2_rs=sum2_prog_rs/(i);
      
        if(i-1==0){ //when we have a single block we set the std equal to zero
          error_prog_rm=0;
          error_prog_rs=0;
        }
        else{
          error_prog_rm=sqrt((mean2_rm-mean_rm*mean_rm)/(i-1));
          error_prog_rs=sqrt((mean2_rs-mean_rs*mean_rs)/(i-1));
        }
        output<<i<<" "<<mean_rm<<" "<<error_prog_rm<<" "<<mean_rs<<" "<<error_prog_rs<<endl;    
   }  
      
   output.close();
   
/**************************************************************************************************/   
//Chi square

   int tries=pow(10,2);
   N=pow(10,4); //number of throws 
   M=pow(10,2); //number of sub-intervals which divide [0,1]
   vector<int> m(M); 
   int mean=N/M; //num of expected events in each sub-interval after N throws
   double chi=0; //chi square
   double inf=0, sup=0;

   //we generate 'N' times a random number between [0,1] we divide this interval in M sub-intervals, we pick up the interval contiaining the random number and we record this result in an histogram m[M], we repeat all 'tries' times   
 
   output.open("chi_square.dat");
   for(int k=0; k<tries; k++){
     for(int i=0; i<N; i++){
       r=rnd.Rannyu(); 
       for(int j=0; j<M; j++){
         inf=(j*1./M);
         sup=(j+1)*1./M;
         if(r>=inf && r<sup){ //pick up the interval
             m[j]+=1; //increase the corresponding bin i.e. the number of observed events in each sub-intervals after N throws    
         } 
         else{
           continue;
         }       
       }     
     } 
         
     for(int l=0; l<M; l++){
       chi+=pow(m[l]-mean, 2)/mean; //evaluate chi square
     }  
     output<<k<<" "<<chi<<endl;   
     fill(m.begin(), m.end(), 0); //clean the histogram
     chi=0;    
   }
   output.close();
   rnd.SaveSeed();
   return 0;
}
