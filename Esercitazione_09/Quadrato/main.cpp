#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include "random.h"
#include "funzioni.h"

using namespace std;

int main(){

 int length_population=100; //lunghezza della popolazione ovvero numero di cromosomi totali
 int length_cromosoma=32;   //lunghezza del numero di geni che formano un cromosoma (=numero totale di città nel percorso)
 int num_generations= 500;  //numero delle generazioni 
 double distanza=0;
 
 ofstream output;
 ofstream out;
 ofstream out_gen_first;
 ofstream out_gen_sec;
 ofstream out_gen_th;

 vector<vector<double>> matrix; //matrice che contiene le posizioni e il nome delle città (x,y,allele) 
 vector <Cromosoma> crom;  //vettore di classi di tipo cromosoma = popolazione, ovvero un sottoinsieme dei possibili percorsi tra città
 
 CreatePosition_InSquare(matrix, length_cromosoma); //creo le posizioni delle città in un quadrato

//inizializzo il vettore popolazione con posizioni e nomi delle città

 for(int i=0; i<length_population; i++){    
   ShufflePosition(matrix);//Shuffle delle righe della matrice posizioni
   crom.push_back(Cromosoma(matrix, length_cromosoma));//creazione della popolazione
 }
 
 Fitness(crom); //Rioridino la popolazione iniziale
 
 cout<<"Best cromosoma iniziale "<<endl;
 crom[0].PrintCromosoma();
 cout<<endl<<"Distanza "<<endl;
 crom[0].PrintDistanza(); 

 NewGeneration new_gen(crom); //Creo un elemento della classe NewGeneration
 
 
 //Stampo su file il cromosoma iniziale migliore
 
 output.open("Best_cromosoma_iniziale.dat");
 
 output << new_gen.GetCromosoma()[0].GetDistanza() <<endl;
 for(int i=0; i<length_cromosoma; i++){
   output << new_gen.GetCromosoma()[0].GetGene()[i].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[i].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[i].GetAllele() << endl;
 } 
 //Riscrivo la prima città perché il percorso deve tornare al punto di partenza
 output << new_gen.GetCromosoma()[0].GetGene()[0].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetAllele(); 
 
 output.close();  
 

 output.open("L2.dat");
 out.open("<L2>.dat"); 
 out_gen_first.open("Best_gen_10.dat");
 out_gen_sec.open("Best_gen_50.dat");
 out_gen_th.open("Best_gen_100.dat");
 
 for(int i=0; i<num_generations; i++){ //loop sulle generazioni
   new_gen.Evolution();  //Creazione di una nuova generazione intera
   output << i << " " << new_gen.GetCromosoma()[0].GetDistanza() << endl; //Stampo su file il valore istantaneo della distanza dei 						cromosomi migliori ad ogni generazione
   
   for (int j=0; j<length_population/2; j++){
     distanza+=new_gen.GetCromosoma()[j].GetDistanza();   
   }
   out << i << " " << distanza*1./(1.*length_population/2.) << endl;  //Stampo su file il valor medio della distanza dei cromosomi ad 					           ogni generazione
   distanza=0;  
   
   if(i==10){
      for(int j=0; j<length_cromosoma; j++){
        out_gen_first << new_gen.GetCromosoma()[0].GetGene()[j].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[j].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[j].GetAllele() << endl;
      } 
 //Riscrivo la prima città perché il percorso deve tornare al punto di partenza
      out_gen_first << new_gen.GetCromosoma()[0].GetGene()[0].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetAllele(); 
   
   }
   
   if(i==50){
     for(int j=0; j<length_cromosoma; j++){
       out_gen_sec << new_gen.GetCromosoma()[0].GetGene()[j].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[j].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[j].GetAllele() << endl;
     } 
 //Riscrivo la prima città perché il percorso deve tornare al punto di partenza
     out_gen_sec << new_gen.GetCromosoma()[0].GetGene()[0].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetAllele(); 
     
   }   
      
   if(i==100){
     for(int j=0; j<length_cromosoma; j++){
       out_gen_th << new_gen.GetCromosoma()[0].GetGene()[j].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[j].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[j].GetAllele() << endl;
     } 
 //Riscrivo la prima città perché il percorso deve tornare al punto di partenza
     out_gen_th << new_gen.GetCromosoma()[0].GetGene()[0].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetAllele(); 
     
   }  
   
 } 
 
 output.close();
 out.close();
 out_gen_first.close();
 out_gen_sec.close();
 
  
 cout<<"Best cromosoma finale "<<endl;
 
 new_gen.Fitness();//Rioridino la popolazione finale
 new_gen.GetCromosoma()[0].PrintCromosoma();
 cout<<endl<<"Distanza "<<endl;
 new_gen.GetCromosoma()[0].PrintDistanza();

 //Stampo su file il cromosoma finale migliore
 
 output.open("Best_cromosoma_finale.dat");
 output << new_gen.GetCromosoma()[0].GetDistanza() <<endl;
 for(int i=0; i<length_cromosoma; i++){
   output << new_gen.GetCromosoma()[0].GetGene()[i].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[i].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[i].GetAllele() << endl;
 } 
 //Riscrivo la prima città perché il percorso deve tornare al punto di partenza
 output << new_gen.GetCromosoma()[0].GetGene()[0].GetX() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetY() << " " << new_gen.GetCromosoma()[0].GetGene()[0].GetAllele(); 
 
 output.close();  
  
 return 0;
}
