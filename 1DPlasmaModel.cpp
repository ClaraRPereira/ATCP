#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <ctime>
#include <algorithm>

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

vector<double> mass;      // Vector for storing particle masses
vector<double> pos;       // Positions
vector<double> vel;       // Velocities
vector <double> X;        // Displacements

FILE *data, *energies;           // Create files to store the particle data and the energies

double Vt;                       // Max velocity for uniform distribution   


//variáveis para imprimir as energias
double E_kin;			        
double E_pot;
double E_tot;

double L=20;  //size of box

bool colision = false;    // Checks if a colision was found
bool print    = true;     // Variable to decide if i want to print stuff

void initial_conditions();     // Inicia as partículas com as condições iniciais
void record_trajectories();    // Writes trajectories to a file
void record_energies();        // Writes energies to a file

void PrintParStatsAll();      // Print position, velocity, and momentum of every particle
void PrintParStats(int i);
void PrintParOrder();

void loop(double dtime , int k);  // Methods  for the main algorithm
void energy( double time );

double func();

int main(){

  srand((int) time(0));  // Seeding the random distribution 

  cout << "\n \t  ****** 1D PLASMA MODEL ****** \n" << endl;

  NPart = 100;   // Number of particles
  Vt = 2;        // Max absolute velocity
  
  // Time parameters
  tmin = 0.0;
  tmax = 0.1;   // Simulation time
  dt = 0.01;    // Time step
  t = tmin;     // Initial time

  // Create the files to store the data
  data = fopen( "DATA", "w" );
  energies = fopen( "ENERGY", "w" );

  initial_conditions(); // Call function to initialize the particles' variables (position, velocity and momenta)

  // Print, to the terminal, positions, velocities & momenta.
  if(print){

    PrintParStatsAll();
    PrintParOrder();
  }

  // Dynamics Iteration
  while ( t < tmax ){

    cout << " TEMPO " << t << endl;

    // Compute energies at current timestep
    energy(t);

    // Compute positions and velocities at current timestep and determine crossing positions
    func();

    if (colision) 
      break; // PARAR O LOOP SE JÁ ENCONTREI A COLISÃO

    // Go to the next timestep
    t = t + dt;
    //PrintParStatsAll();

    // Write the positions of the particles in the file "DATA" if print_trajectory ==1.
    if (print){

      record_trajectories( );
      record_energies( );
    }
  }

  // Close trajectory and energy files
  fclose( data );
  fclose( energies );

  return 0;
}

void
PrintParStatsAll(){            // Prints TO THE TERIMNAL position, velocity and momentum for all particles

  // Print, to the terminal, positions, velocities & momenta.
  cout << " Nº Partículas:  " << NPart << " \n " << endl;

  for (int i = 0; i < NPart; ++i){

    PrintParStats(i);
  }

}

void
PrintParStats(int i){          // Prints TO THE TERIMNAL position, velocity and momentum for particle i

  cout << " Partícula " << i << " \t Position : " << part.x[i] << " |  Velocity :  " ;
  cout << part.vx[i] << "   |  Momentum : " << part.mx[i] << " |" << endl;  
}

void 
PrintParOrder(){              // Prints TO THE TERMINAL particle ordering

  cout << " >>> Ordenação das Partículas : " ;
  for (int i = 0; i < NPart; ++i){

    cout << part.num[i] << " | ";
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

  double m = 1;           // Let's start by defining all masses as 1
  double r1;              // auxiliary variable to generate a random
  double r2;              // auxiliary variable for another random
  double axis = 0;        // Start of x axis
  double spc  = L/NPart;  // Defining intersheet spacing
  
  for (int i = 0 ; i< NPart ; i++ ){  // For loop to set the different initial vectors

    part.num.push_back(i);  // Defining Particle ordering
    npart.num.push_back(i);
    
    //Defining the x positions  
    grid.push_back(axis); 
    pos.push_back(axis);
    axis += spc;

    // random velocities according to uniform distribution
    r1 = (double)rand() / (double)RAND_MAX;  // generating rando between 0 and 1
    r2 = (double)rand() / (double)RAND_MAX;  
    if(r2>=0.5) vel.push_back(Vt*r1);
    else if (r2<0.5) vel.push_back(-Vt*r1);  // Randomly distributing velocities between -Vt and Vt

    X.push_back(0);           // Initializing displacement vector at zero

    // Defining the masses
    mass.push_back(m);
  }

  part.set_values(pos,vel,mass);        // Setting the values for the current particles and the new particles. 
  npart.set_values(pos,vel,mass);
}

void 
loop( double dtime , int a ){

  double Wp = 1;
  double aux;    // auxiliary variable to store old values of stuff

  double media;
/*
  if(a!=0) 
  {
    media=(pos[a+1]+pos[a])*0.5;
    aux=part.vx[a];
    part.vx[a]=part.vx[a+1];
    part.vx[a+1]=aux;
    part.x[a]=media;
    part.x[a+1]=media;
    cout << " nova velocidade da particula " << npart.num[a] << " é " << part.vx[a] << endl;
    cout << " nova velocidade da particula " << npart.num[a+1] << " é " << part.vx[a+1] << endl;
    cout << " nova posiçao da particula " << npart.num[a] << " é " << part.x[a] << endl;
    cout << " nova posiçao da particula " << npart.num[a+1] << " é " << part.x[a+1] << endl;
  } */
  
  for (int i = 0; i < NPart; ++i)
  {
   // X[i]=npart.x[i]-part.x[i];
    X[i]=part.x[i]-pos[i];
  }

  

    // IF THERE AREN'T ANY CROSSINGS JUST NORMALLY CALCULATING TIME EVOLUTION OF "HARMONIC OSCILLATORS"

  for (int i = 0; i < NPart; ++i){

    npart.vx[i]=part.vx[i]*cos(Wp*dtime)-Wp*X[i]*sin(Wp*dtime);
    npart.x[i]=part.x[i]+part.vx[i]*sin(Wp*dtime)-X[i]*(1-cos(Wp*dtime));
    if (npart.x[i]<0) npart.x[i]=0;
  }   

  if(a!=-1) 
  {
    media=(pos[a+1]+pos[a])*0.5;
    aux=npart.vx[a];
    npart.vx[a]=npart.vx[a+1];
    npart.vx[a+1]=aux;
    npart.x[a]=media;
    npart.x[a+1]=media;
    cout << " nova velocidade da particula " << npart.num[a] << " é " << npart.vx[a] << endl;
    cout << " nova velocidade da particula " << npart.num[a+1] << " é " << npart.vx[a+1] << endl;
    cout << " nova posiçao da particula " << npart.num[a] << " é " << npart.x[a] << endl;
    cout << " nova posiçao da particula " << npart.num[a+1] << " é " << npart.x[a+1] << endl;
  }


}    

double 
func(){

  vector<double> b3;
  int b2=0;
  int b=0;
  int n=NPart;
  double dtt=dt;
  double t_c=0; //, t_c2;// Variables for the crossing times, tc1 and tc2
  vector<double> vec_cross;
  double min_tc2=0;
  double c1;
  int num_col=0;

  int col=0;

  double Wp=1;
  double t_c2;
  /*b3.reserve(NPart);
  vec_cross.reserve(NPart);

  npart.x.reserve(n);
  npart.vx.reserve(n);*/
  double store_time=0;


  //vector<double> d;
  //vector<double> c;
  double temp, temp1;
  cout.precision(17);  

  loop(dtt,-1);

  LOOP:
  col=0;
  vec_cross.clear();
  b3.clear();
  for (int i = 0; i < n; ++i)
  {
    for (int j = i+1; j < n; ++j)
    {


      if(npart.x[i]>npart.x[j]  && j==i+1)
      {
        col=1; // Houve pelo menos uma colisão
        num_col=num_col+1; // Quero saber quantas colisões houve
        b=i;
        
        //a=i;
        b2=j;
        

        cout << " CROSSING  entre posições : " <<  b << " e " << b2 ;
        t_c= dtt*(part.x[b2]-part.x[b])/(part.x[b2]-part.x[b]+npart.x[b]-npart.x[b2]);
        cout << " \n ----- 1a APROXIMAÇÃO AO TEMPO DE CROSSING --- TC1 = " << t_c << endl; 

        temp=part.x[b]+part.vx[b]*sin(Wp*t_c)-X[b]*(1-cos(Wp*t_c));
        temp1=part.x[b2]+part.vx[b2]*sin(Wp*t_c)-X[b2]*(1-cos(Wp*t_c));
        if(part.x[b2]-part.x[b]+temp-temp1 > 0.00005)
        {
          t_c2= (t_c)*(part.x[b2]-part.x[b])/(part.x[b2]-part.x[b]+temp-temp1);
        //  cout << " o t do sistema é  " << t << " e o tc2 é " << t_c2 << endl;
        }
        else t_c2=t_c;
        // cout << " t_c" << t_c << endl;
        if ( t_c2 < dt && t_c2>0 )
        {
          b3.push_back(b);
          vec_cross.push_back(t_c2);
          cout << " - SEGUNDA APROXIMAÇÃO AO TEMPO DE CROSSING -- TC2 = " << t_c2  <<" e o time step é : " << dt << endl;
        }
    //   else if (t_c2 > dt ) cout << " Deu merda aqui                                                   111111111111111111111" << endl;
      }
    }
  }
/*
for (int i = 0; i < n; ++i)
{
  cout << " velocidades  1  ----  " << i << " " << npart.vx[i] << endl;
}*/

//cout <<" estou a chegar aqui                             1111111111111111111111111111" << endl;

  if(col!=0)
  {

    min_tc2=vec_cross[0];
    c1=0;

    for (int j = 1; j < vec_cross.size(); ++j)
    {

      if(vec_cross[j]<min_tc2 ) 
      { 
        min_tc2=vec_cross[j];
    //cout << "\n min_tc2 " << min_tc2 << endl;
        c1=j;
      }
    }

    cout << " SELECIONEI A COLISÃO " << npart.num[b3[c1]] << " E " <<  npart.num[b3[c1]+1] << endl;

    cout << " ----------- TEMPO DE CROSSING FINAL --------- TC2 = " << min_tc2 << endl;
    store_time=store_time+min_tc2;
    cout << " TEEEMPO INCREMENTADO  " << store_time << endl;

 /* c2=npart.num[b3[c1]];  
  npart.num[b3[c1]]=npart.num[b3[c1]+1];
  npart.num[b3[c1]+1]=c2;*/

  //time=t+min_tc2;
    loop(min_tc2,b3[c1]);
  }



//cout <<" estou a chegar aqui                         2222222222222222222222222222222222222222 vou incrementar " << time << endl;


  //pos = part.x;
  //vel = part.vx;
  //part.num=npart.num;
 /* for (int i = 0; i < n; ++i)
{
  cout << " velocidades  2---- "<< i << " " << npart.vx[i] << endl;
}*/


  part.x = npart.x;
  part.vx= npart.vx;
  vec_cross.clear();
  b3.clear();  

 // time + 0.001 > dt 
//if (col==0) {loop(dt,0); delta =dt; goto LOOP;}

 //t=t+min_tc2;

  if ( dt - store_time < 0.00001)
  {
    cout << " dt - store_time " << dt- store_time<< endl;
    loop(dt-store_time,-1);
    dtt=dt-store_time;
    goto LOOP;
  } 


  return min_tc2;
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
