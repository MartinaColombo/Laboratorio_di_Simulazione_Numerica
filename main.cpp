//Esercizio 5
//Applicazione dell'algoritmo di Metropolis per il calcolo del valor medio e della varianza del raggio per il ground state e lo stato 2p di un atomo di idrogeno

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;
 
int main (){

   Random rnd;
   int seed[4];
   int p1, p2;
   int hit_gs=0, hit_2p=0;
   int M=pow(10,6); //N throws
   int N_blocks=pow(10,2); //N blocks
   
   double x_gs=0.05, y_gs=0.05, z_gs=0.05; //starting point a_0
   double x_new_gs=0, y_new_gs=0, z_new_gs=0;
   double x_2p=0.05, y_2p=0.05, z_2p=0.05; //starting point a_0
   double x_new_2p=0, y_new_2p=0, z_new_2p=0;   
   double r_bohr=0.052917;
   double delta=0.06;
   double psi_gs=0, psi_2p=0, psi_gs_new=0, psi_2p_new=0; 
   double accept_gs=0, accept_2p=0; 
   double rand=0, sum_r_gs=0, sum_r_2p=0;
   double mean_r_gs=0, mean2_r_gs=0, mean_prog_gs=0, mean2_prog_gs=0;
   double mean_r_2p=0, mean2_r_2p=0, mean_prog_2p=0, mean2_prog_2p=0;
   
   ifstream Primes("Primes");
   
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
   ofstream output_gs;
   ofstream output_2p;
   ofstream output;
   
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
   
   /*************************************    Sampling uniforme ground state    *************************************/  

   output_gs.open("r_gs.dat");
   output.open("Config_gs.dat");
   for(int i=0; i<N_blocks; i++){
     sum_r_gs=0; 
     for(int j=1; j<=M/N_blocks; j++){       
       psi_gs=pow(r_bohr, -3)*1/M_PI*exp(-2*sqrt(x_gs*x_gs+y_gs*y_gs+z_gs*z_gs)/r_bohr);

       x_new_gs=rnd.Rannyu(x_gs-delta, x_gs+delta);
       y_new_gs=rnd.Rannyu(y_gs-delta, y_gs+delta);
       z_new_gs=rnd.Rannyu(z_gs-delta, z_gs+delta);
       
       psi_gs_new=pow(r_bohr, -3)*1/M_PI*exp(-2*sqrt(x_new_gs*x_new_gs+y_new_gs*y_new_gs+z_new_gs*z_new_gs)/r_bohr);       
       accept_gs=min(1.,psi_gs_new/psi_gs);                    
       rand=rnd.Rannyu();      
       if(rand<=accept_gs){
         x_gs=x_new_gs;
         y_gs=y_new_gs;
         z_gs=z_new_gs;        
         hit_gs++;               
       }
       sum_r_gs+=sqrt(x_gs*x_gs+y_gs*y_gs+z_gs*z_gs);   
       output<<x_gs<<" "<<y_gs<<" "<<z_gs<<endl;        
     }     
     //Blocking gs     
     mean_r_gs+=sum_r_gs/(M/N_blocks);
     mean2_r_gs+=pow(sum_r_gs/(M/N_blocks),2);
     
     mean_prog_gs=mean_r_gs/(i+1);
     mean2_prog_gs=mean2_r_gs/(i+1);     
     if(i==0){
       output_gs<<i<<" "<<mean_prog_gs<<" "<<0.0<<endl;    
     }
     else{
       output_gs<<i<<" "<<mean_prog_gs<<" "<<sqrt((mean2_prog_gs-mean_prog_gs*mean_prog_gs)/i)<<endl;   
     } 
              
   }
   output_gs.close(); 
   output.close();    
   cout<<"Percentage of acceptance gs: "<<hit_gs*1.*100/M<<endl;
   
   /*************************************    Sampling uniforme 2p    *************************************/       
   
   output_2p.open("r_2p.dat"); 
   output.open("Config_2p.dat");
   for(int i=0; i<N_blocks; i++){
     sum_r_2p=0; 
     for(int j=1; j<=M/N_blocks; j++){   
       psi_2p=pow(r_bohr, -5)*1/(32*M_PI)*z_2p*z_2p*exp(-sqrt(x_2p*x_2p+y_2p*y_2p+z_2p*z_2p)/(r_bohr)); 
       
       x_new_2p=rnd.Rannyu(x_2p-(2.5*delta), x_2p+(2.5*delta));
       y_new_2p=rnd.Rannyu(y_2p-(2.5*delta), y_2p+(2.5*delta));
       z_new_2p=rnd.Rannyu(z_2p-(2.5*delta), z_2p+(2.5*delta));     
       
       psi_2p_new=pow(r_bohr, -5)*1/(32*M_PI)*z_new_2p*z_new_2p*exp(-sqrt(x_new_2p*x_new_2p+y_new_2p*y_new_2p+z_new_2p*z_new_2p)/(r_bohr));
       
       accept_2p=min(1.,psi_2p_new/psi_2p);
       rand=rnd.Rannyu();              
       if(rand<=accept_2p){
         x_2p=x_new_2p;
         y_2p=y_new_2p;
         z_2p=z_new_2p;        
         hit_2p++;               
       }
       sum_r_2p+=sqrt(x_2p*x_2p+y_2p*y_2p+z_2p*z_2p);  
       output<<x_2p<<" "<<y_2p<<" "<<z_2p<<endl;              
     }           
     //Blocking 2p         
     mean_r_2p+=sum_r_2p/(M/N_blocks);
     mean2_r_2p+=pow(sum_r_2p/(M/N_blocks),2);
     
     mean_prog_2p=mean_r_2p/(i+1);
     mean2_prog_2p=mean2_r_2p/(i+1);     
     if(i==0){
       output_2p<<i<<" "<<mean_prog_2p<<" "<<0.0<<endl;    
     }
     else{
       output_2p<<i<<" "<<mean_prog_2p<<" "<<sqrt((mean2_prog_2p-mean_prog_2p*mean_prog_2p)/i)<<endl;   
     } 
   }
   output_2p.close(); 
   output.close();      
   cout<<"Percentage of acceptance 2p: "<<hit_2p*1.*100/M<<endl;
   
   /*************************************    Sampling gauss ground state    *************************************/         
   sum_r_gs=0;
   mean_r_gs=0;
   mean2_r_gs=0;
   mean_prog_gs=0;
   mean2_prog_gs=0;
   x_gs=0.05;
   y_gs=0.05;
   z_gs=0.05;
   delta=0.6;
   x_new_gs=0;
   y_new_gs=0;
   z_new_gs=0;
   hit_gs=0;
   
   output_gs.open("r_gs_Gauss.dat");
   output.open("Config_gs_gauss.dat");
   for(int i=0; i<N_blocks; i++){
     sum_r_gs=0; 
     for(int j=1; j<=M/N_blocks; j++){       
       psi_gs=pow(r_bohr, -3)*1/M_PI*exp(-2*sqrt(x_gs*x_gs+y_gs*y_gs+z_gs*z_gs)/r_bohr);

       x_new_gs=rnd.Gauss(x_gs, delta/15);
       y_new_gs=rnd.Gauss(y_gs, delta/15);
       z_new_gs=rnd.Gauss(z_gs, delta/15);
       
       psi_gs_new=pow(r_bohr, -3)*1/M_PI*exp(-2*sqrt(x_new_gs*x_new_gs+y_new_gs*y_new_gs+z_new_gs*z_new_gs)/r_bohr);       
       accept_gs=min(1.,psi_gs_new/psi_gs);                    
       rand=rnd.Rannyu();      
       if(rand<=accept_gs){
         x_gs=x_new_gs;
         y_gs=y_new_gs;
         z_gs=z_new_gs;        
         hit_gs++;               
       }
       sum_r_gs+=sqrt(x_gs*x_gs+y_gs*y_gs+z_gs*z_gs); 
      output<<x_gs<<" "<<y_gs<<" "<<z_gs<<endl;                     
     }     
     //Blocking gs     
     mean_r_gs+=sum_r_gs/(M/N_blocks);
     mean2_r_gs+=pow(sum_r_gs/(M/N_blocks),2);
     
     mean_prog_gs=mean_r_gs/(i+1);
     mean2_prog_gs=mean2_r_gs/(i+1);     
     if(i==0){
       output_gs<<i<<" "<<mean_prog_gs<<" "<<0.0<<endl;    
     }
     else{
       output_gs<<i<<" "<<mean_prog_gs<<" "<<sqrt((mean2_prog_gs-mean_prog_gs*mean_prog_gs)/i)<<endl;   
     }          
   }
   output_gs.close();  
   output.close();   
   cout<<"Percentage of acceptance gs (Gauss): "<<hit_gs*1.*100/M<<endl;
   
    /*************************************    Sampling gauss 2p    *************************************/         
   sum_r_2p=0;
   mean_r_2p=0;
   mean2_r_2p=0;
   mean_prog_2p=0;
   mean2_prog_2p=0;
   x_2p=0.05;
   y_2p=0.05;
   z_2p=0.05;
   delta=0.6;
   x_new_2p=0;
   y_new_2p=0;
   z_new_2p=0;
   hit_2p=0;   
   
   output_2p.open("r_2p_Gauss.dat"); 
   output.open("Config_2p_gauss.dat");   
   for(int i=0; i<N_blocks; i++){
     sum_r_2p=0; 
     for(int j=1; j<=M/N_blocks; j++){   
       psi_2p=pow(r_bohr, -5)*1/(32*M_PI)*z_2p*z_2p*exp(-sqrt(x_2p*x_2p+y_2p*y_2p+z_2p*z_2p)/(r_bohr)); 
       
       x_new_2p=rnd.Gauss(x_2p, delta/6.5);
       y_new_2p=rnd.Gauss(y_2p, delta/6.5);
       z_new_2p=rnd.Gauss(z_2p, delta/6.5);     
       
       psi_2p_new=pow(r_bohr, -5)*1/(32*M_PI)*z_new_2p*z_new_2p*exp(-sqrt(x_new_2p*x_new_2p+y_new_2p*y_new_2p+z_new_2p*z_new_2p)/(r_bohr));
       
       accept_2p=min(1.,psi_2p_new/psi_2p);
       rand=rnd.Rannyu();              
       if(rand<=accept_2p){
         x_2p=x_new_2p;
         y_2p=y_new_2p;
         z_2p=z_new_2p;        
         hit_2p++;               
       }
       sum_r_2p+=sqrt(x_2p*x_2p+y_2p*y_2p+z_2p*z_2p);  
      output<<x_2p<<" "<<y_2p<<" "<<z_2p<<endl;              
     }           
     //Blocking 2p         
     mean_r_2p+=sum_r_2p/(M/N_blocks);
     mean2_r_2p+=pow(sum_r_2p/(M/N_blocks),2);
     
     mean_prog_2p=mean_r_2p/(i+1);
     mean2_prog_2p=mean2_r_2p/(i+1);     
     if(i==0){
       output_2p<<i<<" "<<mean_prog_2p<<" "<<0.0<<endl;    
     }
     else{
       output_2p<<i<<" "<<mean_prog_2p<<" "<<sqrt((mean2_prog_2p-mean_prog_2p*mean_prog_2p)/i)<<endl;   
     } 
   }
   output_2p.close(); 
   output.close();      
   cout<<"Percentage of acceptance 2p: "<<hit_2p*1.*100/M<<endl;   
     
   
   rnd.SaveSeed();
   return 0;
}

