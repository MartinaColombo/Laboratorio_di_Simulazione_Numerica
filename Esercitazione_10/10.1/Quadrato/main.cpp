#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <random>
#include "funzioni.h"

using namespace std;


 int main(){

  int length_cromosoma=32;   //lunghezza del numero di geni che formano un cromosoma (=numero totale di città nel percorso)
  int tries= 7000;    //numero di mutazioni effettuate prima di un cambio della temperatura
  double T=15;    //Temperatura iniziale
  double step=0.015;    //Step tra due temperature successive 
  double beta=0;
  double accept=0, r=0, k=0;
  int counter=0, steps=T/step;
 
  ofstream output;
  ofstream out;
  ofstream out_accept;
  ifstream input;
  
  random_device rand;
  default_random_engine gen(rand());
  uniform_real_distribution<> unif_dist(0., 1.);  

  vector<vector<double>> matrix; //matrice che contiene le posizioni e il nome delle città (x,y,allele) 
  
  input.open("Best_cromosoma_iniziale.dat"); //Leggo le posizioni delle città utilizzate nell'esercitazione 7
  
  int dim=3;
  vector<double> temporary;
  double count=0, element=0;
   
  while(count!=length_cromosoma){
    while(temporary.size()!=dim){
      input>>element;
      temporary.push_back(element);
    }
    matrix.push_back(temporary);
    temporary.clear();
    ++count;
  }  
     
  input.close(); 
  
  Cromosoma crom(matrix, length_cromosoma); //inizializzo un cromosoma con le posizioni e città create 
  Cromosoma temp;
  
  cout<<endl<<"Distanza iniziale ";
  crom.PrintDistanza(); 
  cout<<endl<<endl;    
  
  output.open("Cromosoma_iniziale.dat");
 
  output << crom.GetDistanza() <<endl;
  for(int i=0; i<length_cromosoma; i++){
    output << crom.GetGene()[i].GetX() << " " << crom.GetGene()[i].GetY() << " " << crom.GetGene()[i].GetAllele() << endl;
  } 
  output << crom.GetGene()[0].GetX() << " " << crom.GetGene()[0].GetY() << " " << crom.GetGene()[0].GetAllele() << endl;  
  output.close();
  
  output.open("temp&distance.dat");
  out.open("iteration&distance.dat");
  out_accept.open("temp&acceptance.dat"); 
  
  for(int i=1; i<=steps; i++){
    if(i%100==0){
      cout<<"Iterazione "<<i<<endl;
    }
    beta=1./(1.*T);
    k=0;
     
     for(int j=0; j<tries; j++){  
       temp = crom; //tengo memoria del cromosoma precedente
        
       r=unif_dist(gen); 
       if(r<0.35){
         temp.Pair_Permutation();  
       }    
       r=unif_dist(gen); 
       if(r<0.25){
         temp.Reverse_Permutation();  
       }     
       r=unif_dist(gen);
       if(r<0.20){
         temp.Contiguous_Permutation(); 
       }
       r=unif_dist(gen);
       if(r<0.30){
         temp.Swap_Permutation(); 
       }
       
       //Metropolis
       r=unif_dist(gen);
       accept=min(1.,exp(-beta*(temp.GetDistanza()-crom.GetDistanza())));
         
       if(r<=accept){
        crom=temp; //accetto la mossa
        k++;     
       }       
     }
     out_accept << T <<" "<<k<<endl;    
     //if(counter%10 == 0){         
     out<<counter<<" "<<crom.GetDistanza()<<endl;
     //}  
     output<<T<<" "<<crom.GetDistanza()<<endl; 	       
     counter++; 
     T=T-step;       
  }
  
  output.close();
  out.close();
  out_accept.close();  
  
  cout<<endl<<"Distanza finale ";
  crom.PrintDistanza(); 
  cout<<endl<<endl;
   
  output.open("Cromosoma_finale.dat");
  output << crom.GetDistanza() <<endl;
  for(int i=0; i<length_cromosoma; i++){
    output << crom.GetGene()[i].GetX() << " " << crom.GetGene()[i].GetY() << " " << crom.GetGene()[i].GetAllele() << endl;
  }
  output << crom.GetGene()[0].GetX() << " " << crom.GetGene()[0].GetY() << " " << crom.GetGene()[0].GetAllele() << endl;   
  output.close(); 
 
 
 return 0;
 }
