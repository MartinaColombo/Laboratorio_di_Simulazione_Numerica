#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <random>
#include "random.h"
#include "funzioni.h"

using namespace std;

//Classi

/******************************************************* Gene  *******************************************************/

 Gene::Gene(double x, double y, int allele){
   _x=x; 
   _y=y; 
   _allele=allele;
 }
 
 Gene::Gene(const Gene &originale){
   _x=originale.GetX();
   _y=originale.GetY();
   _allele=originale.GetAllele();
 }
 
 Gene& Gene::operator=(const Gene &originale){
   _x=originale.GetX();
   _y=originale.GetY();
   _allele=originale.GetAllele();
   return *this;
 }
 
 bool Gene::operator==(const Gene& confronto)const{
   if((_x==confronto.GetX()) &&
      (_y==confronto.GetY()) &&
      (_allele==confronto.GetAllele()))
     return true;   
   else
     return false; 
 }

 bool Gene::operator!=(const Gene& confronto)const{
   if((_x!=confronto.GetX()) &&
      (_y!=confronto.GetY()) &&
      (_allele!=confronto.GetAllele()))
     return true;   
   else
     return false;  
 }

 void Gene::SetPosition(double x, double y){
   _x=x; 
   _y=y;
 }

 void Gene::SetAllele(int allele){
   _allele=allele;
 }
 
 void Gene::PrintGene()const{
   cout << _x << " " << _y << " " << _allele <<endl;
 }
 
 
 /**************************************************  Cromosoma  *****************************************************/

 Cromosoma::Cromosoma(vector<vector<double>> matrix, int length){
   _length=length;
   _gene.resize(length);
   
   for(int i=0; i<length; i++){
      int j=0;
     _gene[i].SetPosition(matrix[i][j], matrix[i][j+1]);  //matrix è una matrice lengthx3 data da (x,y,allele)
     _gene[i].SetAllele(matrix[i][j+2]);  //il vettore di alleli ha come prima e ultima città 1
   }
    SetDistanza();
 }
 
 Cromosoma::Cromosoma(vector <Gene> gene){
   _gene=gene;
   _length=gene.size(); 
   
   double distanza=0;
   
   for(int i=0; i<gene.size(); i++){ 
   
    int j=i+1; 
    if(i==gene.size()-1){ 
      distanza+=(pow(gene[i].GetX()-gene[0].GetX(),2) + pow(gene[i].GetY()-gene[0].GetY(),2));
    }
    else{ 
      distanza+=(pow(gene[i].GetX()-gene[j].GetX(),2) + pow(gene[i].GetY()-gene[j].GetY(),2));      
    }   
   }
   _distanza=distanza;
 }
 
 Cromosoma::Cromosoma(const Cromosoma &originale){  
   _length=originale.GetLength();
   _distanza=originale.GetDistanza();
   _gene=originale.GetGene();
   
 }

 Cromosoma& Cromosoma::operator=(const Cromosoma &originale){  
   _length=originale.GetLength();
   _distanza=originale.GetDistanza();
   _gene=originale.GetGene(); 
   return *this;  
 }

 bool Cromosoma::operator==(const Cromosoma& confronto)const{
   if((_length==confronto.GetLength()) &&
      (_distanza==confronto.GetDistanza()) &&
      (_gene==confronto.GetGene()))
     return true;   
   else
     return false;  
 }

 bool Cromosoma::operator!=(const Cromosoma& confronto)const{
   if((_length!=confronto.GetLength()) &&
      (_distanza!=confronto.GetDistanza()) &&
      (_gene!=confronto.GetGene()))
     return true;   
   else
     return false;  
 }

 void Cromosoma::SetLength(int length){ 
   _length=length;
 }
 
 void Cromosoma::PrintCromosoma() const{
  for(auto i:_gene){
      i.PrintGene();    
  }
 }
 
 void Cromosoma::SetDistanza(){
   
   double distanza=0;
   
   for(int i=0; i<_length; i++){ 
    int j=i+1; 
    if(i==_length-1){ 
      distanza+=(pow(_gene[i].GetX()-_gene[0].GetX(),2) + pow(_gene[i].GetY()-_gene[0].GetY(),2));
    }
    else{ 
      distanza+=(pow(_gene[i].GetX()-_gene[j].GetX(),2) + pow(_gene[i].GetY()-_gene[j].GetY(),2));      
    }   
   }
   _distanza=distanza;
 }

 void Cromosoma::PrintDistanza()const{
   cout<<_distanza <<endl;
 }
 
/**********************************************  NewGeneration  **********************************************/ 
 
  
 NewGeneration::NewGeneration(vector <Cromosoma> &crom){
   _crom=crom;
   SetLengthPopulation(crom.size());
   SetLengthCromosoma(crom[0].GetLength());
 }
 
 NewGeneration::NewGeneration(const NewGeneration &originale){
   _crom=originale.GetCromosoma(); 
 }
 
 NewGeneration& NewGeneration::operator=(const NewGeneration &originale){
   _crom=originale.GetCromosoma(); 
   return *this; 
 }

 void NewGeneration::SetLengthPopulation(int length_population){
   _length_population=length_population;
 }
 
 
 void NewGeneration::SetLengthCromosoma(int length_cromosoma){
   _length_cromosoma=length_cromosoma;
 }
 
 void NewGeneration::PrintGeneration()const{
   for(auto i:_crom){
     i.PrintCromosoma(); 
     i.PrintDistanza();
     cout<<endl<<endl;   
   } 
 }
 
//In fitness ordino dal peggiore al minore secondo la norma 

   void NewGeneration::Fitness(){
   //riordino tutta la popolazione di cromosomi dal migliore al peggiore tramite una lambda function
   sort(_crom.begin(),_crom.end(),[](const Cromosoma &lhs,const Cromosoma &rhs){
     return lhs.GetDistanza() < rhs.GetDistanza();
   }); 
 }
 
//Funzione che seleziona l'indice dei genitori migliori dal punto di vista della distanza
 
 int NewGeneration::Selezione(){
   
  double r=0;
  double p=0;
  int index=0;
   
  random_device rnd;
  default_random_engine gen(rnd());
  uniform_real_distribution<> unif_r(0.,1.);
  uniform_real_distribution<> unif_p(0.5,11.); 
  r=unif_r(gen);
  p=2.5;     
  index=int(GetLengthPopulation()*(pow(r,p)));
      
  return index;   
 }
 
 void NewGeneration::Crossover(int m, int d){
 
   vector <Gene> mum = _crom[m].GetGene();
   vector <Gene> dad = _crom[d].GetGene();
   vector <Gene> temp_one = _crom[m].GetGene();
   vector <Gene> temp_two = _crom[d].GetGene();
   
   random_device rnd;
   default_random_engine gen(rnd());
   uniform_int_distribution<> unif(1,mum.size()-2);
   
   int cut =  unif(gen);
   temp_one.erase(temp_one.begin()+cut, temp_one.end());
   temp_two.erase(temp_two.begin()+cut, temp_two.end());
   
//Primo figlio   
   
   for(int i=1; i<dad.size(); i++){  
     if(find(temp_one.begin(), temp_one.end(), dad[i]) != temp_one.end()){  //Se dad[i] è presente nel vettore non fare niente   
       continue;
     }     
     else{  //Se dad[i] non è presente nel vettore mettilo in coda al figlio
       temp_one.push_back(dad[i]);     
     }   
   }
   
//Secondo figlio   
      
   for(int i=1; i<mum.size(); i++){  
     if(find(temp_two.begin(), temp_two.end(), mum[i]) != temp_two.end()){  //Se mum[i] è presente nel vettore non fare niente
       continue;
     }     
     else{  //Se mum[i] non è presente nel vettore mettilo in coda al figlio
       temp_two.push_back(mum[i]);     
     }   
   }
   
   Cromosoma son_one(temp_one);
   Cromosoma son_two(temp_two);
   
   _crom.push_back(son_one);
   _crom.push_back(son_two);

 }
 
 
 void NewGeneration::Pair_Permutation(){ //Swap of two cities (except for the first one)
 
   random_device rnd;
   default_random_engine gen(rnd());
   uniform_int_distribution<> index_population(0,GetLengthPopulation()-1); 
   uniform_int_distribution<> index_gene(1,_crom[0].GetLength()-1);

   int first=0, second=0, index=0;
   first=index_gene(gen);

   second=index_gene(gen);
   
   index=index_population(gen); //Scelgo l'indice del cromosoma da mutare
   while(first==second){
     second=index_gene(gen);        
   }
   
   vector <Gene> temp =_crom[index].GetGene();  
   
   swap(temp[first], temp[second]);//Scambio due città aventi indici first e second
   Cromosoma Permutate(temp);
   Permutate.SetDistanza();
   _crom[index]=Permutate;      
 }
 
 void NewGeneration::Swap_Permutation(){ //Scambio di m città contigue con altre m città contigue all'interno del cromosoma
   int m=0, index=0;
   random_device rnd;
   default_random_engine gen(rnd());
   
   uniform_int_distribution<> index_population(0,GetLengthPopulation()-1); 
   index=index_population(gen);//index of the cromosome 
   
   
   uniform_int_distribution<> index_m(1, _crom[index].GetLength()/2 -1); 
   m=index_m(gen); //number of contiguous cities to swap
   
   vector <Gene> temp =_crom[index].GetGene();  
   
   rotate(temp.begin() + 1, temp.begin() + m, temp.end());
   
   Cromosoma Permutate(temp);
   Permutate.SetDistanza();
   
   _crom[index]=Permutate;  
 
 }
 

 void NewGeneration::Contiguous_Permutation(){

   int m=0, index=0;
   random_device rnd;
   default_random_engine gen(rnd());
   
   uniform_int_distribution<> index_population(0,GetLengthPopulation()-1); 
   index=index_population(gen);//index of the cromosome 
   
   
   uniform_int_distribution<> index_m(1, _crom[index].GetLength()/2 -1); 
   m=index_m(gen); //number of contiguous cities to swap

   
   vector <Gene> temp =_crom[index].GetGene();  
   
   rotate(temp.begin() + 1, temp.begin() + m, temp.begin() + 2*m);
   
   Cromosoma Permutate(temp);
   Permutate.SetDistanza();
   
   _crom[index]=Permutate;
 }
 
 
 void NewGeneration::Reverse_Permutation(){
 
   random_device rnd;
   default_random_engine gen(rnd());
   uniform_int_distribution<> index_population(0,GetLengthPopulation()-1); 
   uniform_int_distribution<> index_gene(1,_crom[0].GetLength()-1); 
    
   int first=0, second=0, index=0;
   
   first=index_gene(gen);
   second=index_gene(gen);
   
   while(first>=second){
     first=index_gene(gen);  
     second=index_gene(gen);            
   }
   
   index=index_population(gen); //Scelgo l'indice del cromosoma da mutare
   vector <Gene> temp =_crom[index].GetGene();

   
   reverse(temp.begin()+first, temp.begin()+second+1);//Scambio due città aventi indici first e second
   Cromosoma Permutate(temp);
   Permutate.SetDistanza();
   
   _crom[index]=Permutate;   
 
 
 }

 void NewGeneration::Evolution(){
 
   int old_size=_crom.size();
   int mum=0, dad=0; //genero gli indici dei cromosomi genitori per il crossover
   random_device rnd;
   default_random_engine gen(rnd());
   uniform_real_distribution<> unif(0.,1.);   
   double r=0;
   
   for(int i=0; i<old_size/2; i++){ //Costruisco la nuova generazione
     
     r=unif(gen);
     mum=Selezione();
     dad=Selezione();
     while(mum==dad){
       dad=Selezione();  
     }  
     
     if(r<0.8){    
       Crossover(mum,dad); 
     }     
     else{     
        _crom.push_back(_crom[mum]);
        _crom.push_back(_crom[dad]);    
     }   
   }
   _crom.erase(_crom.begin(), _crom.begin()+old_size);  
 
   //Mutazioni  
   
  for(int i=0; i< _length_population;i++){
  
    r=unif(gen); 
    if(r<0.10){
      Pair_Permutation();  
    }    
    r=unif(gen); 
    if(r<0.10){
      Reverse_Permutation();  
    }     
    r=unif(gen);
    if(r<0.10){
      Contiguous_Permutation(); 
    }
    r=unif(gen);
    if(r<0.10){
      Swap_Permutation(); 
    }
  }
   Fitness(); //per prima cosa riordino la popolazione in base alla fitness
 }


 
/**********************************************  Funzioni utili  **********************************************/
 
 
//Funzione per creare posizioni generiche in un quadrato di lato 1 con rispettivonome della città [1,2,3,4,5]
 
 void CreatePosition_InSquare(vector<vector<double>> &matrix, int length_cromosoma){

   //cout<<"Creating positions "<<endl;
   double x=0, y=0;
   
   random_device rnd;
   default_random_engine gen(rnd());
   uniform_real_distribution<> unif(-1.,1.);   
   
   for(int i=0; i<length_cromosoma; i++){
     x=unif(gen);
     y=unif(gen);
     matrix.push_back({x, y,(i+1)*1.});   //(x,y,allele)         
   }   
 }
 
//Funzione per disordinare le città ad esclusione della prima

 void ShufflePosition(vector<vector<double>> & matrix){
   random_shuffle(matrix.begin()+1,matrix.end()); 
 }

//Funzione per riordinare le città in funzione della distanza quadratica 

 void Fitness(vector <Cromosoma> &_crom){
   //riordino tutta la popolazione di cromosomi dal migliore al peggiore tramite una lambda function
   sort(_crom.begin(),_crom.end(),[](const Cromosoma &lhs,const Cromosoma &rhs){
     return lhs.GetDistanza() < rhs.GetDistanza();
   }); 
 }
 

 
