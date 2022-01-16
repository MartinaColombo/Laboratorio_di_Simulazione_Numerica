//Esercizio 2.3 
//3D Random Walk (RW) in un reticolo cubico e nel continuo 

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   int a=1; //step length
   int N_walk=pow(10,4); //number of RWs
   int N_step=pow(10,2); //number of steps in a RW
   double pos[3]={0,0,0}; //position vector
   int direct=0; //direction along axes x,y,z
   double way=0;  //forward or backward 
   double dist=0; //|r_n|^2
   double error_mean_pos=0;
   double theta=0, phi=0;
   
   vector <double> mean_pos(N_step); //vector of the mean position after n steps <|r_n|^2>
   vector <double> mean2_pos(N_step); //vector of the mean square position after n steps <|r_n|^2*|r_n|^2>
   
   
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

//Cubic lattice
   
   mean_pos[0]=0;
   mean2_pos[0]=0;
   for(int i=0; i<N_walk; i++){
      pos[0]=0; //always start from the origin
      pos[1]=0;
      pos[2]=0;
      for(int j=0; j<N_step-1; j++){
         direct=int(rnd.Rannyu(0,3));
         if(rnd.Rannyu()<0.5){
            way=1;
         }
         else{
            way=-1;
         }
         pos[direct]+=a*way;
         dist=pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
         mean_pos[j+1]+=dist;
         mean2_pos[j+1]+=dist*dist;      
      } 
   }

   output.open("RW_lattice.dat");
   for(int i=0; i<N_step; i++){
      mean_pos[i]=mean_pos[i]/N_walk;
      mean2_pos[i]=mean2_pos[i]/N_walk;
      if(i==0){
         error_mean_pos=0;
      }
      else{
      error_mean_pos=0.5*sqrt((mean2_pos[i]-pow(mean_pos[i],2))/((N_walk-1)*mean_pos[i]));
      }
      output<<i<<" "<<sqrt(mean_pos[i])<<" "<<error_mean_pos<<endl;
   }

   output.close();


//Continuum

   fill(mean_pos.begin(), mean_pos.end(), 0);     //fill the vectors with zeros
   fill(mean2_pos.begin(), mean2_pos.end(), 0);
   
   mean_pos[0]=0;
   mean2_pos[0]=0;
   for(int i=0; i<N_walk; i++){
      pos[0]=0;  //always start from the origin
      pos[1]=0;
      pos[2]=0;
      for(int j=0; j<N_step-1; j++){
         theta=rnd.Theta();	//sampling solid angle
         phi=rnd.Rannyu(0,2*M_PI);
         pos[0]+=a*sin(theta)*cos(phi);
         pos[1]+=a*sin(theta)*sin(phi);
         pos[2]+=a*cos(theta);
         dist=pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
         mean_pos[j+1]+=dist;
         mean2_pos[j+1]+=dist*dist;      
      } 
   } 

   output.open("RW_continuum.dat");
   for(int i=0; i<N_step; i++){
      mean_pos[i]=mean_pos[i]/N_walk;
      mean2_pos[i]=mean2_pos[i]/N_walk;
      if(i==0){
         error_mean_pos=0;
      }
      else{
      error_mean_pos=0.5*sqrt((mean2_pos[i]-pow(mean_pos[i],2))/((N_walk-1)*mean_pos[i]));
      }
      output<<i<<" "<<sqrt(mean_pos[i])<<" "<<error_mean_pos<<endl;
   }

   output.close();

   rnd.SaveSeed();
   return 0;
}

