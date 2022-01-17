#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <random>
#include "random.h"
#include "funzioni.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[]){


 int size, rank;
 MPI_Init(&argc,&argv);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 MPI_Status stat1, stat2;
 //MPI_Request req;
 double tstart = MPI_Wtime();


 int length_population=100; //lunghezza della popolazione ovvero numero di cromosomi totali
 int length_cromosoma=32;   //lunghezza del numero di geni che formano un cromosoma (=numero totale di città nel percorso)
 int num_generations = 500;

 vector<vector<double>> matrix; //matrice che contiene le posizioni e il nome delle città (x,y,allele)
 vector <Cromosoma> crom;  //vettore di classi di tipo cromosoma = popolazione, ovvero un sottoinsieme dei possibili percorsi tra città
 ifstream in;

 vector <int> allele1(length_cromosoma);
 vector <int> allele2(length_cromosoma);

 int first=0, second=0;
 int itag=1, itag2=2;

 Random rnd;
 int seed[4];
 int p1, p2;
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

 //Ogni nodo inizializza la stessa popolazione di individui

 in.open("Posizioni.dat");
 int dim=3;
 vector<double> temporary;
 double count=0, element=0;

 while(count!=length_cromosoma){
   while(temporary.size()!=dim){
     in>>element;
     temporary.push_back(element);
   }
   matrix.push_back(temporary);
   temporary.clear();
   ++count;
 }

 in.close();

 Cromosoma sample(matrix, length_cromosoma);

 //inizializzo il vettore popolazione con posizioni e nomi delle città

 for(int i=0; i<length_population; i++){
   ShufflePosition(matrix);//Shuffle delle righe della matrice posizioni
   crom.push_back(Cromosoma(matrix, length_cromosoma));//creazione della popolazione
 }

 Fitness(crom); //Rioridino la popolazione iniziale

  //cout<<"Distanza del cromosoma mugliore del nodo "<<rank<<endl;
  //crom[0].PrintDistanza();

 NewGeneration* new_gen;
 new_gen = new NewGeneration(crom);

 for(int i=1; i<=num_generations; i++){ //loop sulle generazioni
   //new_gen.Evolution();  //Creazione di una nuova generazione intera
   new_gen -> Evolution();

    if(i%20==0){//Scambio del cromosoma migliore ogni 10 generazioni
      first=rnd.Rannyu_INT(0,3);
      second=rnd.Rannyu_INT(0,3);

      while(first==second){
        second=rnd.Rannyu_INT(0,3);
      }

      //cout<<"Primo e secondo "<<first<<" "<<second<<endl;

      //cout << "Lo scambio avviene tra il nodo " << first << " e il nodo " << second << endl;

      if(rank==first){

        for(int k=0; k<length_cromosoma; k++){
          allele1[k]=new_gen->GetCromosoma()[0].GetGene()[k].GetAllele();
          allele2[k]=new_gen->GetCromosoma()[0].GetGene()[k].GetAllele();
        }
        /*cout<<"Sono il nodo "<<rank<<" stampo gli alleli migliori "<<endl;
        for(int j=0; j<length_cromosoma; j++){
          cout<<allele1[j]<<" "<<allele2[j]<<endl;
        }*/

      }

      if(rank==second){

        for(int k=0; k<length_cromosoma; k++){
          allele1[k]=new_gen->GetCromosoma()[0].GetGene()[k].GetAllele();
          allele2[k]=new_gen->GetCromosoma()[0].GetGene()[k].GetAllele();
        }
       /* cout<<"Sono il nodo "<<rank<<" stampo gli alleli migliori "<<endl;
        for(int j=0; j<length_cromosoma; j++){
          cout<<allele1[j]<<" "<<allele2[j]<<endl;
        }*/

      }


      /*if(rank==first){
        cout<<"Sono il nodo "<<rank<<" invio a "<<second<<" il cromosoma"<<endl;
        //new_gen.GetCromosoma()[0].PrintCromosoma();
        new_gen -> GetCromosoma()[0].PrintCromosoma();
      }
      if(rank==second){
        cout<<"Sono il nodo "<<rank<<" invio a "<<first<<" il cromosoma"<<endl;
        //new_gen.GetCromosoma()[0].PrintCromosoma();
        new_gen -> GetCromosoma()[0].PrintCromosoma();
      }  */

      if(rank==first){
        MPI_Send(&allele1[0], length_cromosoma, MPI_INTEGER, second, itag, MPI_COMM_WORLD);
        MPI_Recv(&allele2[0], length_cromosoma, MPI_INTEGER, second, itag2, MPI_COMM_WORLD, &stat2);
      }
      else if(rank==second){
        MPI_Recv(&allele1[0], length_cromosoma, MPI_INTEGER, first, itag, MPI_COMM_WORLD, &stat1);
        MPI_Send(&allele2[0], length_cromosoma, MPI_INTEGER, first, itag2, MPI_COMM_WORLD);
      }

      if(rank==first){
        Cromosoma crom_first(sample, allele2);
        //new_gen.SetCromosoma(0, crom_first);
        new_gen -> SetCromosoma(0, crom_first);
      }
      if(rank==second){
        Cromosoma crom_second(sample, allele1);
        //new_gen.SetCromosoma(0, crom_second);
        new_gen -> SetCromosoma(0, crom_second);
      }

     /* if(rank==first){
        cout<<"Sono il nodo "<<rank<<" come miglior cromosoma ora ho "<<endl;
        //new_gen.GetCromosoma()[0].PrintCromosoma();
        new_gen -> GetCromosoma()[0].PrintCromosoma();
        cout<<endl<<endl;
      }
      if(rank==second){
        cout<<"Sono il nodo "<<rank<<" come miglior cromosoma ora ho "<<endl;
        //new_gen.GetCromosoma()[0].PrintCromosoma();
        new_gen -> GetCromosoma()[0].PrintCromosoma();
        cout<<endl<<endl;
      } */
     }
     Measure(*new_gen, rank, i, length_population, length_cromosoma, num_generations);
 }


 cout<<endl;
 cout<<"Best cromosoma finale ottenuto dal nodo "<<rank<<" ha distanza: ";
 //new_gen.Fitness();//Rioridino la popolazione finale
 //new_gen.GetCromosoma()[0].PrintDistanza();
 new_gen -> Fitness();
 new_gen -> GetCromosoma()[0].PrintDistanza();
 cout<<endl<<endl;

 double tend = MPI_Wtime();
 double dt = tend - tstart;
 cout<<"Rank "<<rank<<", tempo di esecuzione: "<<dt<<endl<<endl;

 delete new_gen;
 MPI_Finalize();

 return 0;
}
