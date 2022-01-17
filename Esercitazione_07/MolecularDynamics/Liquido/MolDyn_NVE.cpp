#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 

  Input();             //Inizialization
  int nconf = 1;
 
  for(int iblk=1; iblk<=nblk; iblk++){ 
    Reset(iblk);  
    for(int istep=1; istep <= nstep; ++istep){
      Move();           //Move particles with Verlet algorithm
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
         //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!           
         nconf += 1;    
      }    
    }
     Averages(iblk);   //Print results for current block         
  }
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;
  
  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl;   
  
  

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4; //Correzioni viriale
  n_props = 5; //Number of observables #4
  
  
//measurement of g(r)
  igofr = 5; //#4
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;  
  
  if(restart==0){
//Read initial configuration 
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;   
  }
  ReadConf.close(); 
    cout << "Read initial configuration from file config.0 " << endl << endl;
//Prepare initial velocities
    cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i){
      vx[i] = rand()/double(RAND_MAX) - 0.5;
      vy[i] = rand()/double(RAND_MAX) - 0.5;
      vz[i] = rand()/double(RAND_MAX) - 0.5;

      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
      double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
     }    
  }
  else if(restart==1){
//Read initial configuration
  cout << "Read initial configuration from file last.final " << endl << endl;
  ReadConf.open("last.final");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i] >> xold[i] >> yold[i] >> zold[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
      
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box; 
  }
  ReadConf.close(); 
  
    Move();
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){ //v[t+delta/2]=(r(t+dt)-r(t))/2dt
      vx[i]=(x[i]-xold[i])/(2*delta);
      vy[i]=(y[i]-yold[i])/(2*delta);
      vz[i]=(z[i]-zold[i])/(2*delta);
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }  
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for (int i=0; i<npart; ++i){
      vx[i] *=fs;
      vy[i] *=fs;
      vz[i] *=fs;
    
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);    
  
    }
  }    
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew;

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  //int bin;
  double v, w, t, vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("output_epot.dat", ios::app);
  Ekin.open("output_ekin.dat", ios::app);
  Temp.open("output_temp.dat", ios::app);
  Etot.open("output_etot.dat", ios::app);
  Press.open("output_press.dat", ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;
  
//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;  

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     
//update of the histogram of g(r)
     for(int l=igofr; l<igofr+nbins; l++){
       if(dr >= (l-5)*bin_size && dr < (l-4)*bin_size){
         walker[l]+=2;
       }
     }       

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
//Potential energy and virial
       v += vij;
       w += wij;       
     }
    }          
  }
  
  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;  
  
//Kinetic energy
  for (int i=0; i<npart; ++i) t += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  
  walker[ik]=0.5*t;
  walker[it]=2./3. * t;
  walker[ie]=walker[ik]+walker[iv];
   
    stima_pot = walker[iv]/(double)npart + vtail; //Potential energy per particle
    stima_kin = walker[ik]/(double)npart; //Kinetic energy per particle
    stima_temp = walker[it]/(double)npart; //Temperature
    stima_etot = walker[ie]/(double)npart; //Total energy per particle
    stima_press = rho*temp + (walker[iw] + ptail*npart)/vol;

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press <<endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
  ofstream WriteOut;
  ofstream WriteToUse;
    if(restart==0){
    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    for (int i=0; i<npart; ++i){
      WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box <<" "<< xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }
    WriteConf.close();   
    }
    else if(restart==1){
    cout << "Print final configuration to file last.final " << endl << endl;
    WriteOut.open("last.final");
    WriteToUse.open("to_use.final");
    for (int i=0; i<npart; ++i){
      WriteOut << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box <<" "<< xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
      
      WriteToUse << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
      
    }
    WriteOut.close();  
    WriteToUse.close();  
    }
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

/************************************************************/


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
    
   double r, delta_v;
   ofstream Gofr, Gave, Epot, Press, Ekin, Etot, Temp;
    
    cout << "Block number " << iblk << endl;
    
    Temp.open("output.temp.0",ios::app);
    Etot.open("output.etot.0",ios::app);
    Ekin.open("output.ekin.0",ios::app);
    Epot.open("output.epot.0",ios::app);
    Press.open("output.press.0",ios::app);
    Gofr.open("output.gofr.0",ios::app);
    Gave.open("output.gave.0",ios::app);      
    
    stima_temp = blk_av[it]/blk_norm/(double)npart;
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_kin = blk_av[ik]/blk_norm/(double)npart; 
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
    
    stima_etot = blk_av[ie]/blk_norm/(double)npart; 
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
    
    stima_press = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_press;
    glob_av2[iw] += stima_press*stima_press;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);
    
    for(int l=igofr; l<n_props; l++){
      delta_v=(pow((l-3)*bin_size,3)-pow((l-4)*bin_size,3)); 
      stima_gdir = blk_av[l]/blk_norm/(4*M_PI/3*delta_v*rho*npart);
      glob_av[l] += stima_gdir;
      glob_av2[l] += stima_gdir*stima_gdir;
      err_gdir = Error(glob_av[l], glob_av2[l], iblk);
      r = (l-4) * bin_size + bin_size/2;     
      Gofr << iblk <<  " " << r <<  " " << glob_av[l]/(double)iblk << endl;       
        if(iblk==nblk){
           Gave << r <<  " " << glob_av[l]/(double)iblk<<" "<<err_gdir<< endl; 
        }
    }
//Potential energy per particle
    Epot << iblk << " "<< stima_pot << " "<< glob_av[iv]/(double)iblk <<" "<< err_pot << endl;
    
    Temp << iblk << " "<< stima_temp << " "<< glob_av[it]/(double)iblk << " " << err_temp << endl;
    
    Ekin << iblk <<  " "<< stima_kin << " "<< glob_av[ik]/(double)iblk <<  " "<< err_kin << endl;
    
    Etot << iblk << " "<< stima_etot << " " << glob_av[ie]/(double)iblk << " " << err_etot << endl;
    
//Pressure
    Press << iblk <<  " " << stima_press << " " << glob_av[iw]/(double)iblk << " " << err_press << endl;

//g(r)
 

    cout << "----------------------------" << endl << endl;


    Temp.close();
    Etot.close();
    Ekin.close();
    Epot.close();
    Press.close();
    Gofr.close();
    Gave.close();
}



double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
