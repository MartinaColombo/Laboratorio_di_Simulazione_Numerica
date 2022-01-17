#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <random>
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
 
 
  void Cromosoma::Pair_Permutation(){ //Swap of two cities (except for the first one)
 
   random_device rnd;
   default_random_engine gen(rnd());
   uniform_int_distribution<> index_gene(1,_length-1);

   int first=0, second=0;
   first=index_gene(gen);
   second=index_gene(gen);
   
   while(first==second){
     second=index_gene(gen);        
   }
   
   vector <Gene> temp =_gene;  
   
   swap(temp[first], temp[second]);//Scambio due città aventi indici first e second
   _gene=temp;
   SetDistanza();
    
 }
 
 void Cromosoma::Swap_Permutation(){ //Scambio di m città contigue con altre m città contigue all'interno del cromosoma
   int m=0;
   random_device rnd;
   default_random_engine gen(rnd());   
   
   uniform_int_distribution<> index_m(1, _length/2 -1); 
   m=index_m(gen); //number of contiguous cities to swap
   
   vector <Gene> temp =_gene;  
   
   rotate(temp.begin() + 1, temp.begin() + m, temp.end());
   _gene=temp;
   SetDistanza();   
 }
 

 void Cromosoma::Contiguous_Permutation(){

   int m=0;
   random_device rnd;
   default_random_engine gen(rnd());
   
   uniform_int_distribution<> index_m(1, _length/2 -1); 
   m=index_m(gen); //number of contiguous cities to swap
   
   vector <Gene> temp =_gene;  
   
   rotate(temp.begin() + 1, temp.begin() + m, temp.begin() + 2*m);
   _gene=temp;
   SetDistanza();
 }
 
 
 void Cromosoma::Reverse_Permutation(){
 
   random_device rnd;
   default_random_engine gen(rnd());
   uniform_int_distribution<> index_gene(1,_length-1); 
    
   int first=0, second=0;
   
   first=index_gene(gen);
   second=index_gene(gen);
   
   while(first>=second){
     first=index_gene(gen);  
     second=index_gene(gen);            
   }
   
   vector <Gene> temp =_gene;
   
   reverse(temp.begin()+first, temp.begin()+second+1);//Scambio due città aventi indici first e second
   _gene=temp;
   SetDistanza(); 
 
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
 

 
