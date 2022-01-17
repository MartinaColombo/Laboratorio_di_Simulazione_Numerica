#ifndef __funzioni_h__
#define __funzioni_h__
#include <iostream>

using namespace std;

 class Gene{

  public:
 
    Gene(){_allele=0, _x=0, _y=0;};     //costruttore di default
    Gene(double, double, int);          //costruttore con parametri
    ~Gene(){};                          //distruttore
    Gene(const Gene &);                 //copy constuctor, copio le informazioni di un gene su un altro nuovo gene
    Gene& operator=(const Gene &);      //copy assignment
    bool operator==(const Gene&)const; //overload == operator
    bool operator!=(const Gene&)const; //overload != operator
   
    void SetPosition(double, double);    //setto le posizioni x e y
    void SetAllele(int);	           //setto il nome del gene
   
    double GetX()const{return _x;} //accedo alla coordinata x 
    double GetY()const{return _y;} //accedo alla coordinata y
    int GetAllele()const{return _allele;} //accedo all'allele ovvero il nome della città
    
    void PrintGene()const; //stampo il gene ovvero le sue coordinate x e y e il suo nome
   
  private:
 
    int _allele;
    double _x;
    double _y;

 };


 class Cromosoma{

  public:
   
    Cromosoma(){_length=0, _distanza=0;};    //costruttore di default
    Cromosoma(vector<vector<double>>, int);  //costruttore con parametri
    Cromosoma(vector <Gene>);
    ~Cromosoma(){}; 		//distruttore
    Cromosoma(const Cromosoma &);	//copy constuctor
    Cromosoma& operator=(const Cromosoma &); //copy assignment
    bool operator==(const Cromosoma&)const; //overload == operator
    bool operator!=(const Cromosoma&)const; //overload != operator   
    void SetLength(int);		//setto la lunghezza del cromosoma
    double GetLength()const{return _length;}  //leggo la lunghezza del cromosoma
    void PrintCromosoma()const;	//stampo il cromosoma
    void SetDistanza();		//calcolo la distanza tra le componenti del cromosoma
    double GetDistanza()const{return _distanza;}   //ottengo la distanza

    vector <Gene> GetGene()const{return _gene;}    //ottengo il vettore di geni (=cromosoma)
   
    void PrintDistanza()const;		      //stampo la distanza

  private:
   
    int _length; //lunghezza del cromosoma
    vector <Gene> _gene; //vector di geni
    double _distanza;
 };
 
 
 class NewGeneration{
  
  public:
  
  NewGeneration(){}     //costruttore di default
  NewGeneration(vector <Cromosoma> &);  //costruttore con parametri
  ~NewGeneration(){};		//distruttore
  NewGeneration(const NewGeneration &);      //copy constuctor
  NewGeneration& operator=(const NewGeneration &);  //copy assignment 
  void SetLengthPopulation(int);
  int GetLengthPopulation()const{return _length_population;}
  void SetLengthCromosoma(int);
  int GetLengthCromosoma()const{return _length_cromosoma;}  
  
  vector <Cromosoma> GetCromosoma()const{return _crom;}  //ottengo il vettore di cromosomi (=popolazione)
  
  void PrintGeneration()const;
  
  void Fitness();
  int Selezione();
  void Crossover(int, int);
  void Pair_Permutation(); 
  void Swap_Permutation();
  void Contiguous_Permutation(); 
  void Reverse_Permutation();
  void Evolution();
  
  private:
  
  vector <Cromosoma> _crom; 
  int _length_population;
  int _length_cromosoma;
 };
 
 //Funzioni utili
 
 void CreatePosition_OnCircumference(vector<vector<double>> &, int);
 void ShufflePosition(vector<vector<double>> &);
 void Fitness(vector <Cromosoma> &);
 
#endif
 
