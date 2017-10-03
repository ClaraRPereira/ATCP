#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <ctime>

#include "particles.h"

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

using namespace std;

// Variáveis iniciadas no particles.cpp
extern particles part;           // Initializing one set of particles at time t
extern particles npart;          // Set new set of particles at time t+dt 
extern int NPart;                // Number of particles 

vector<double> grid;             // vector for storing equilibrium positions on the x axis
double dt, tmin, tmax;           // Time grid
double t;                        // Time variable (evolving time)  

vector<double> mass(NPart);      // Vector for storing particle masses
vector<double> pos(NPart);       // Positions
vector<double> vel(NPart);       // Velocities

FILE *data, *energies;           // Create files to store the particle data and the energies

double Vt;                       // Max velocity for uniform distribution   


//variáveis para imprimir as energias
double E_kin;			        
double E_pot;
double E_tot;

double L = 4;  //size of box

bool colision = false;    // Checks if a colision was found
bool print    = true;     // Variable to decide if i want to print stuff

void initial_conditions();      // Inicia as partículas com as condições iniciais
void record_trajectories();     
void record_energies();

void PrintParStatsAll();      // Print position, velocity, and momentum of every particle
void PrintParStats(int i);
void PrintParOrder();

void loop(double dtime , int k);
void energy( double time );

double func();

int main(){

  srand((int) time(0));  // Seeding the random distribution 

  cout << "\n \t  ****** 1D PLASMA MODEL ****** \n" << endl;

  NPart = 100; // Number of particles

  int k=0;
  double tc2;   // position of crossing
  Vt = 5; // Max velocity
  // Time parameters
  tmin = 0.0;
  tmax = 1;
  dt =0.1;
  t = tmin; // Initial time

  // Create the files to store the data
  data = fopen( "DATA", "w" );
  energies = fopen( "ENERGY", "w" );


  //mass.reserve(NPart);

  initial_conditions(); // Call function to initialize the particles' variables (position, velocity and momenta)

  // Print, to the terminal, positions, velocities & momenta.
  if(print){

    PrintParStatsAll();
    PrintParOrder();
  }

  // Dynamics Iteration
  while ( t < tmax ){

    k=k+1;
    cout << " TEMPO " << t;

    // Compute energies at current timestep
    energy(t);

    // Compute positions and velocities at current timestep and determine crossing positions
    tc2 = func();

    if (colision) 
      break; // PARAR O LOOP SE JÁ ENCONTREI A COLISÃO

    // Go to the next timestep
    t = t + dt;

      // Write the positions of the particles in the file "DATA" if print_trajectory ==1.
    if (print == 1){

      record_trajectories( );
      record_energies( );
    }

  }

  if(print){

    PrintParOrder();
  }

  // Close trajectories file
  fclose( data );
  fclose( energies );

  return 0;
}

void
PrintParStatsAll(){

  // Print, to the terminal, positions, velocities & momenta.
  cout << " Nº Partículas:  " << NPart << " \n " << endl;

  for (int i = 0; i < NPart; ++i){

      PrintParStats(i);
  }

}

void
PrintParStats(int i){

  cout << " Partícula " << i << " \t Position : " << part.x[i] << " |  Velocity :  " ;
  cout << part.vx[i] << "   |  Momentum : " << part.mx[i] << " |" << endl;  
}

void 
PrintParOrder(){

  cout << " >>> Ordenação das Partículas : " ;
  for (int i = 0; i < NPart; ++i){

    cout << part.num[i] << " ,";
  }

  cout << endl;
}

// Write the positions, velocities and momenta of the particles in the file " DATA ".
void record_trajectories( ){

  int n=NPart;

  for (int i=0 ; i<n ; i++ ){

    fprintf( data, "%f %f %f %f\n", t, part.x[i], part.vx[i], part.mx[i] ); // Escrever os valores para o ficheiro "DATA"
  }
}

// Write the energies of the system in the file " ENERGY ".
void record_energies(){

  fprintf( energies, "%f, %f, %f, %f\n",t, E_kin, E_tot, E_pot);
}

void initial_conditions(){ // Gonna Define the initial Conditions

  vector<double> Px; 

  double m = 1; //Let's start by defining all masses as 1
  //double sigma=0.5, n0=0.7; //Valores aleatórios para a carga por unidade de área, sigma, e para a density of neutralizing background charges 
  //double Wp= 4*M_PI*sigma*sigma*n0/m; //Plasma Frequency
  double r1;     // auxiliary variable to generate a random
  double axis = 0; // Start of x axis
  double spc  = 0.1; //defining intersheet spacing
  
  for (int i = 0 ; i< NPart ; i++ ){

    part.num.push_back(i);
    npart.num.push_back(i);
    //Defining the x positions  
    //part.x[i]+=spc;
    grid.push_back(axis); 
    pos.push_back(axis);
    //cout << " pos: " << pos[i] << endl;
    axis += spc;

    // random velocities according to uniform distribution
    r1 = (double)rand() / (double)RAND_MAX;  // generating rando between 0 and 1
    vel.push_back(Vt * r1);
    // cout << " vel: " << vel[i] << endl;

    // Defining the masses
    mass.push_back(m);
  }

  part.set_values(pos,vel,mass);
}

void 
loop( double dtime , int k){

  int a=0; // Store position of crossing particle
  double dtt=dtime;
  //double m = 1; //Let's start by defining all masses as 1   ************ PERIGO porque usavas o m para o for ***********
  //double sigma=0.5, n0=0.7; //Valores aleatórios para a carga por unidade de área, sigma, e para a density of neutralizing background charges 
  //double Wp= 4*M_PI*sigma*sigma*n0/m; //Plasma Frequency
  double Wp = 0.5;
  vector <double> X;
  npart.x.reserve( NPart);
  npart.vx.reserve( NPart);
  X.reserve( NPart);
  double t_c1, Delta_c2;// Variables for the crossing times, tc1 and tc2

  double c1;
  vector<double> vec_cross;
  vector<double> d;
  vector<double> c;

  if (k==3){
    for (int i = 0;  i < NPart; ++i){

      pos[i] = part.x[i];
      vel[i] = part.vx[i];
      part.x[i] = npart.x[i];
      part.vx[i]= npart.vx[i];
    } 
  }

  for (int i = 0; i < NPart; ++i)
  {
    X[i]=part.x[i]-grid[i];
  }

    // IF THERE AREN'T ANY CROSSINGS JUST NORMALLY CALCULATING TIME EVOLUTION OF "HARMONIC OSCILLATORS"

  for (int i = 0; i < NPart; ++i){

    npart.vx[i]=part.vx[i]*cos(Wp*dtt)-Wp*X[i]*sin(Wp*dtt);
    npart.x[i]=part.x[i]+part.vx[i]*sin(Wp*dtt)-X[i]*(1-cos(Wp*dtt));
  }   
}    

double 
func(){

  int a=0; // Store position of crossing particle
  int k=0;
  int b;
  int b2;
  double dtt=dt;
  double dt2;
  int n=NPart;
  double m=1; //Let's start by defining all masses as 1
  //double sigma=0.5, n0=0.7; //Valores aleatórios para a carga por unidade de área, sigma, e para a density of neutralizing background charges 
  //double Wp= 4*M_PI*sigma*sigma*n0/m; //Plasma Frequency
  double Wp=0.5;
  vector <double> X;
  npart.x.reserve(n);
  npart.vx.reserve(n);
  X.reserve(n);
  double t_c1, Delta_c2;// Variables for the crossing times, tc1 and tc2


  double c1;
  vector<double> vec_cross;
  vector<double> d;
  vector<double> c;
  double temp, temp1;
  
  // int o=0; ???? 

  loop(dtt,0);

  // LOOPS TO LOOK FOR CROSSINGS
  LOOP:  for ( int i=0 ; i<n ; i++ ){
    for ( int  j=i+1 ; j<n ; j++  ){      //j=i+1

  //   if(npart.x[j]<npart.x[i])  xixi=xixi+1;
      if(npart.x[i]>npart.x[j]){ 

        a=1;
        k=k+1;
        b=npart.num[i];
        b2=npart.num[j];

        cout << " CROSSING  entre posições : " <<  b << " e " << b + 1 ;//c.push_back(b);

        d.push_back(i);
        c.push_back(b); 


        t_c1= dtt*(part.x[b2]-part.x[b])/(part.x[b2]-part.x[b]+npart.x[b]-npart.x[b2]);

        cout << " \n ----- 1a APROXIMAÇÃO AO TEMPO DE CROSSING --- TC1 = " << t_c1 << endl; 
        cout << " t " << t;

        temp=part.x[b]+part.vx[b]*sin(Wp*t_c1)-X[b]*(1-cos(Wp*t_c1));
        temp1=part.x[b2]+part.vx[b2]*sin(Wp*t_c1)-X[b2]*(1-cos(Wp*t_c1));
        Delta_c2= (t_c1-t)*(part.x[b2]-part.x[b])/(part.x[b2]-part.x[b]+temp-temp1);

        cout << "\n ----------- TEMPO DE CROSSING FINAL --------- TC2 = " << Delta_c2 << endl;

        dt2=t+Delta_c2;

        vec_cross.push_back(Delta_c2);
        if (dt2>=t+dt) 
          goto exit;
        else {
          loop(dt2,0); 
          loop(t+dt-Delta_c2,3);
          goto LOOP;
        }
      // else if (k%2==0 && k!=0){ cout << "\n ----------- TEMPO DE CROSSING FINAL --------- TC2 = " << Delta_c << endl;}
      }

      else if (a==0) 
        goto exit;// Particle i collides with particle j . 
    }
  }

  exit: //if (a==0) // if there are no crossings we just store old values and refresh the new ones
   // {
  for (int i = 0; i < n; ++i){

    pos[i]=part.x[i];
    vel[i]=part.vx[i];
    part.x[i]=npart.x[i];
    part.vx[i]=npart.vx[i];

  } 
  //  }  

  if (a==1){

    for (int i = 0; i < d.size(); i++){

      c1=npart.num[d[i]];  

      npart.num[d[i]]=npart.num[d[i]+1];
      npart.num[d[i]+1]=c1;

    }
  }

  else if (a==0) 
    cout << " \t \t \t \t Partículas não chocaram " << endl;

  int y=0;
  for (int i = 0; i < c.size(); i++){

    colision = true;
    y=y+1;
    cout << "\n \t EVENTO " << y << ": \t \t Partícula " << c[i] << " chocou com partícula " << c[i]+1 << endl;
    cout << " TEMPO FINAL DE CROSSING  DO EVENTO  " << y << ": tc2= " << vec_cross[i] << endl;
  }

  part.num=npart.num;
  return t_c1;  
}

void energy( double time ){

  double kinetic, potential, pot, etotal;

  double sigma=0.5, n0=0.7; //Valores aleatórios para a carga por unidade de área, sigma, e para a density of neutralizing background charges 
  double E0= 4*M_PI*sigma*n0; // Amplitude campo Eléctrico
  double xij, xij2;

  // Kinetic energy of the system
  /*kinetic = 0.0;
  for ( i=0 ; i<NPart ; i++ ){

    kinetic = kinetic + 0.5 * part.mx[i]*part.vx[i];

  }*/

  for (int i = 0; i < NPart; ++i){

    kinetic = 0.5*part.dir_product(part.mx,part.vx)[i];
  }

  // Potential Energy of the system
  potential = 0.0;
  for (int i = 0; i< NPart; i++){

    for (int j= i + 1; j< NPart; j++){
      xij = part.x[i] - part.x[j];

      xij2 = xij*xij;

      pot = - E0 * ( xij2/2 );
      potential = potential + pot;
    }
  }

  // Total energy of the system
  etotal = kinetic + potential;

  // Print the energies and current timestep
  //cout << " \t " << time << "    "  << kinetic << "    "  << potential << "    "  << etotal  << endl;

  E_kin=kinetic;
  E_pot=potential;
  E_tot=etotal;
}
