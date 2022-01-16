#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 

  Input();             //Inizialization
  int nconf = 1, n_blocks=pow(10,2), count=0;
  double sum_pot=0, sum_kin=0, sum_temp=0, sum_etot=0;
  double mean_pot=0, mean_kin=0, mean_temp=0, mean_etot=0;
  double mean2_pot=0, mean2_kin=0, mean2_temp=0, mean2_etot=0;
  double mean_prog_pot=0, mean_prog_kin=0, mean_prog_temp=0, mean_prog_etot=0;
  double mean2_prog_pot=0, mean2_prog_kin=0, mean2_prog_temp=0, mean2_prog_etot=0;
    
  ofstream Epot, Ekin, Etot, Temp;
  ofstream Ave_pot, Ave_kin, Ave_temp, Ave_etot;
  
  Epot.open("output_epot.dat");
  Ekin.open("output_ekin.dat");
  Temp.open("output_temp.dat");
  Etot.open("output_etot.dat");  
  
  Ave_pot.open("ave_epot.out");
  Ave_kin.open("ave_ekin.out");
  Ave_temp.open("ave_temp.out");
  Ave_etot.open("ave_etot.out");  
 
  for(int i=0; i<n_blocks; i++){ 
    sum_pot=0;
    sum_kin=0;
    sum_temp=0;
    sum_etot=0;  
    count=0;      
    cout << "Block number "<< i+1 <<endl;
    for(int istep=1; istep <= nstep/n_blocks; ++istep){
       Move();           //Move particles with Verlet algorithm
    // if(istep%iprint == 0)// cout << "Number of time-steps: " << istep << endl;
       if(istep%10 == 0){
          Measure();     //Properties measurement
	  //Printing the instant value of observables in output files
          Epot << stima_pot  << endl;
          Ekin << stima_kin  << endl;
          Temp << stima_temp << endl;
          Etot << stima_etot << endl;          
          
          sum_pot+=stima_pot;
          sum_kin+=stima_kin;
          sum_temp+=stima_temp;
          sum_etot+=stima_etot;                               
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!       
          nconf += 1;
          count++;          
       }    
    }
    //Blocking
    mean_pot+=sum_pot/count;
    mean_kin+=sum_kin/count;
    mean_temp+=sum_temp/count;
    mean_etot+=sum_etot/count; 
    
    mean2_pot+=pow(sum_pot/count,2);
    mean2_kin+=pow(sum_kin/count,2);
    mean2_temp+=pow(sum_temp/count,2);
    mean2_etot+=pow(sum_etot/count,2); 
             
    mean_prog_pot=mean_pot/(i+1); 
    mean_prog_kin=mean_kin/(i+1);
    mean_prog_temp=mean_temp/(i+1);
    mean_prog_etot=mean_etot/(i+1);  
    
    mean2_prog_pot=mean2_pot/(i+1); 
    mean2_prog_kin=mean2_kin/(i+1);
    mean2_prog_temp=mean2_temp/(i+1);
    mean2_prog_etot=mean2_etot/(i+1);
    
    //Printing average values in output files     
    if(i==0){
      Ave_pot<< i<< " "<<mean_prog_pot<<" "<<0.0<<endl;
      Ave_kin<< i<< " "<<mean_prog_kin<<" "<<0.0<<endl;
      Ave_temp<< i<< " "<<mean_prog_temp<<" "<<0.0<<endl;
      Ave_etot<< i<< " "<<mean_prog_etot<<" "<<0.0<<endl;
    }else{
      Ave_pot<< i<< " "<<mean_prog_pot<<" "<<sqrt((mean2_prog_pot-mean_prog_pot*mean_prog_pot)/(i))<<endl;
      Ave_kin<< i << " "<<mean_prog_kin<<" "<<sqrt((mean2_prog_kin-mean_prog_kin*mean_prog_kin)/(i))<<endl;
      Ave_temp<< i << " "<<mean_prog_temp<<" "<<sqrt((mean2_prog_temp-mean_prog_temp*mean_prog_temp)/(i))<<endl;
      Ave_etot<< i << " "<<mean_prog_etot<<" "<<sqrt((mean2_prog_etot-mean_prog_etot*mean_prog_etot)/(i))<<endl;      
    }      
  }
 
  ConfFinal();         //Write final configuration to restart
    
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  
  Ave_pot.close();
  Ave_kin.close();
  Ave_temp.close();
  Ave_etot.close(); 

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf, ReadConf_final;
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
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables 
  
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
  cout << "Restart enabled" <<endl;
  cout << "Read initial configuration from file old.0 " << endl;
  cout << "Read configuration at time t-dt from file old.final " << endl << endl;
  ReadConf.open("old.0");  //Reading configuration at actual time t
  ReadConf_final.open("old.final"); //Reading configuration at time t-dt
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    ReadConf_final >> xold[i] >> yold[i] >> zold[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
      
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box; 
  }
  ReadConf.close(); 
  ReadConf_final.close();
  
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
  double v, t, vij;
  double dx, dy, dz, dr;
  /*ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat", ios::app);
  Ekin.open("output_ekin.dat", ios::app);
  Temp.open("output_temp.dat", ios::app);
  Etot.open("output_etot.dat", ios::app);*/

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

   /*Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();*/

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
  ofstream WriteOut;
  WriteOut.open("old.0"); //Writing config at time t
  WriteConf.open("old.final"); //Writing configuration at time t-dt
    
  cout << "Print previous configuration to old.final " << endl;
  cout << "Print final configuration to old.0" << endl << endl;  
    
  for (int i=0; i<npart; ++i){
    WriteOut << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close(); 
  WriteOut.close();      
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
