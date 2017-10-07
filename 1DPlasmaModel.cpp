#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
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

// Create files to store the particle data and the energies
ofstream data;
ofstream energies;
ofstream energies_normalized;
ofstream fin_vel;
ofstream ini_vel;
ofstream vel1;
ofstream estatisticas_vel;


double Vt;                       // Max velocity for uniform distribution   


//variáveis para imprimir as energias
double E_kin;			        
double E_pot;
double E_tot;
double E_kin_0;			        
double E_pot_0;
double E_tot_0;

double dE=0.8;

double L;  //size of box

bool colision = false;    // Checks if a colision was found
bool print    = false;     // Variable to decide if i want to print stuff

void initial_conditions();     // Inicia as partículas com as condições iniciais
void record_trajectories();    // Writes trajectories to a file
void record_energies();        // Writes energies to a file
void record_energies_normalized();
void write_velocities(int i);

void PrintParStatsAll();      // Print position, velocity, and momentum of every particle
void PrintParStats(int i);
void PrintParOrder();

void loop(double dtime , int k);  // Methods  for the main algorithm
void energy( double time );

void velocity_statistics(double time );

double func();

int main(){





  srand((int) time(0));  // Seeding the random distribution 

  cout << "\n \t  ****** 1D PLASMA MODEL ****** \n" << endl;

  NPart = 1000;   // Number of particles
  L=0.05*NPart;
  Vt = 1;        // Max absolute velocity
  
  // Time parameters
  tmin = 0.0;
  tmax = 2000;   // Simulation time
  dt = 0.001;    // Time step
  t = tmin;     // Initial time

  // Create the files to store the data
  
  
data.open ("data.dat");
energies.open ("energy.dat");
energies_normalized.open ("energy_normalized.dat");
fin_vel.open ("final_velocity.dat");
ini_vel.open ("initial_velocity.dat");
vel1.open ("vel1.dat");
estatisticas_vel.open("stats.dat");

 // Call function to initialize the particles' variables (position, velocity and momenta)
   initial_conditions();


  double op=0;
  // Dynamics Iteration
  while ( t < tmax ){

    cout << " TEMPO DA SIMULAÇÃO " << t << endl;
    op = op + dt;


    // Compute positions and velocities at current timestep and determine crossing positions
    func();

    if (colision) 
      break; // PARAR O LOOP SE JÁ ENCONTREI A COLISÃO

        // Compute energies at current timestep
  
    energy(t);
	velocity_statistics(t);
	
	cout <<"energia_total_inicial " <<E_tot_0 << endl;
	cout <<"energia_pot_inicial " <<E_pot_0 << endl;
	cout <<"energia_cin_inicial " <<E_kin_0 << endl;
	
    if (t == tmin) {

      write_velocities(1);
    }
    // Go to the next timestep
    t = t + dt;
    
    if(op>=1) {

      write_velocities(3);
      op=0;
    }
    cout << "k " << t << "  ss  " << op << endl;
    //PrintParStatsAll();
    
   

    if (t >= tmax) write_velocities(2);
    
    record_energies( );
    record_energies_normalized( );

    // Write the positions of the particles in the file "DATA" if print_trajectory ==1.
    if (print){

      record_trajectories( );
      
    }
  }


	data.close();
	energies.close();
	energies_normalized.close();
	fin_vel.close();
	ini_vel.close();
	vel1.close();


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

for (int i=0 ; i<n ; i++ )
{
	data << t <<"\t"<< part.x[i] <<"\t"<< part.vx[i] << endl; 
}	
  data << "\n" ;
}

// Write the energies of the system in the file " ENERGY ".
void record_energies(){
 
  energies << t <<"\t"<< E_kin <<"\t"<< E_pot <<"\t"<< E_tot << endl;
  
}

void record_energies_normalized(){
 
  energies_normalized << t <<"\t"<< (E_kin-E_kin_0)/E_kin_0 <<"\t"<< (E_pot-E_pot_0)/E_pot_0 <<"\t"<< (E_tot-E_tot_0)/E_tot_0 << endl;
  
}



void write_velocities(int j){

  for (int i=0 ; i<NPart ; i++ ){

   if(j==2) fin_vel <<  part.vx[i] << endl;
   else if(j==1) ini_vel <<  part.vx[i] << endl;
   else if(j==3)  vel1 <<   part.vx[i] << endl;
   
 }

}

void initial_conditions(){// Gonna Define the initial Conditions

  vector<double> Px; 

  double m = 1;           // Let's start by defining all masses as 1
  double r1;              // auxiliary variable to generate a random
  double axis = 0;        // Start of x axis
  double spc  = L/NPart;  // Defining intersheet spacing
  
  for (int i = 0 ; i< NPart ; i++ ){         // For loop to set the different initial vectors

    part.num.push_back(i);                   // Defining Particle ordering
    npart.num.push_back(i);
    
    //Defining the x positions  
    grid.push_back(axis); 
    pos.push_back(axis);
    axis += spc;

    // random velocities according to uniform distribution
    r1 = (double)rand() / (double)RAND_MAX;  // generating rando between 0 and 1
    vel.push_back(Vt-2*Vt*r1);  // Randomly distributing velocities between -Vt and Vt

    X.push_back(0);                          // Initializing displacement vector at zero

    // Defining the masses
    mass.push_back(m);
  }

  part.set_values(pos,vel,mass);             // Setting the values for the current particles and the new particles. 
  npart.set_values(pos,vel,mass);
}

void 
loop( double dtime , int a ){

  double Wp = 1;
  double aux;    // auxiliary variable to store old values of stuff

  double mean;   // AUXILIARY CHEATING VALUE TO CALCULATE MEAN VALUE BETWEEN PARTICLES POSITIONS WHEN THEY COLLIDE
                 // se a conservação de energia for má. em vez de usarmos o valor médio iteramos um tc3 até as posições de igualarem
  //double temp,temp1;
  //double t_c;
  
  for (int i = 0; i < NPart; ++i){   // Loop to compute displacements

    X[i]=part.x[i]-pos[i];
  }
 /*
if(a!=-1){               ISTO ERA CODIGO DE MIM A TENTAR MERDAS PARA OPTIMIZAR A ENERGIA

  while (temp-temp1 > 0.00005*L/NPart  )
  {
    temp=part.x[a]+part.vx[a]*sin(Wp*dtime)-X[a]*(1-cos(Wp*dtime));
    temp1=part.x[a+1]+part.vx[a+1]*sin(Wp*dtime)-X[a+1]*(1-cos(Wp*dtime));
      
    t_c= (dtime)*(part.x[a+1]-part.x[a])/(part.x[a+1]-part.x[a]+temp-temp1);
    dtime=t_c;
  }
}
*/
  // CALCULATING TIME EVOLUTION OF THE "HARMONIC OSCILLATOR" Equation. Solution to differential equation for each particle
  for (int i = 0; i < NPart; ++i){

    npart.vx[i]=part.vx[i]*cos(Wp*dtime)-Wp*X[i]*sin(Wp*dtime);    // npart are the new values for the particles. part are the old values (previous time step).
    npart.x[i]=part.x[i]+part.vx[i]*sin(Wp*dtime)-X[i]*(1-cos(Wp*dtime));
    //if (npart.x[i]<=0) npart.vx[i]=0;        // só está aqui porque ainda não fizemos as condiçoes de fronteira periódicas. 
    if (npart.x[i]>=L || npart.x[i]<=0) npart.vx[i]=-npart.vx[i];        // só está aqui porque ainda não fizemos as condiçoes de fronteira periódicas. 
  }   

  double vel_max;
  double vel_min;

  if(a!=-1){ // if a=-1 means the new npart values will still be tested for collisions and are not meant to be stored yet.
             // When a!=0, a assumes the value of the position of the colliding particle. 
    /*if(npart.vx[a]>npart.vx[a+1]){

     vel_max=npart.vx[a]; 
     vel_min=npart.vx[a+1];
     npart.vx[a]=npart.vx[a]-dE*(vel_max-vel_min);
     npart.vx[a+1]=npart.vx[a+1]+dE*(vel_max-vel_min);
   }
   else if(npart.vx[a]<npart.vx[a+1]){ 

    vel_max=npart.vx[a+1]; 
    vel_min=npart.vx[a];
    npart.vx[a+1]=npart.vx[a+1]-dE*(vel_max-vel_min);
    npart.vx[a]=npart.vx[a]+dE*(vel_max-vel_min);
  }
  else if ( npart.vx[a]==npart.vx[a+1])
  {   */    
    aux=npart.vx[a];
    npart.vx[a]=npart.vx[a+1];             //  Velocities are switched with its neighbour (elastic collision)
    npart.vx[a+1]=aux;
  //}

  mean=(npart.x[a+1]+npart.x[a])*0.5;                                        
  //aux=npart.vx[a];
  //  npart.vx[a]=npart.vx[a+1];             //  Velocities are switched with its neighbour (elastic collision)
  //  npart.vx[a+1]=aux;
    npart.x[a]=mean;                       // X positions of colliding particles are the exact spot where they collide at. (mean position at t=dt)
    npart.x[a+1]=mean;

    /*
    cout << " nova velocidade da particula " << npart.num[a] << " é " << npart.vx[a] << endl;
    cout << " nova velocidade da partícula " << npart.num[a+1] << " é " << npart.vx[a+1] << endl;
    cout << " nova posiçao da partícula " << npart.num[a] << " é " << npart.x[a] << endl;
    cout << " nova posiçao da partícula " << npart.num[a+1] << " é " << npart.x[a+1] << endl;
    */
  }


}    

double 
func(){

  vector<double> col_pos;  // stores positions of colliding particles
  int cpar1=0;             // position of 1st colliding particle
  int cpar2=0;             // position of second colliding particle 
  int n=NPart;             // storing total no. os particles
  int col=0;               // col=0 if there aren't any collisions. col=1 if there are collisions
  int j=0;                 // Iterator
  int num_shocks=0; // NUmber of collisions

  double Wp=1;             // Plasma freq.
  double dtt=dt;           // Time step
  double t_c=0;            // Variable for the crossing time, tc1 
  double t_c2;             // Variable for the crossing time, tc2 
  double min_tc2=0;        // Minimum tc2 (time of first crossing event)
  double minp;             // Position of particle colliding at min_tc2
  double store_time=0;     // Stores sum of time steps made so far
  double temp, temp1;      // Auxiliary 

  vector<double> vec_cross;   // Vector stores tc2 values
  double cross_size;          // size of vector vec_cross (just to avoid warnings and always compare integer values)

  loop(dtt,-1);            // First loop, calculates new npart positions and velocities for my time step dt

  LOOP:
  col=0;
  vec_cross.clear();
  col_pos.clear();
  for (int i = 0; i < n-1; ++i){

    j=i+1;
    if(npart.x[i]>=npart.x[j] ){
      num_shocks=num_shocks+1;
      //cout << "Number collisions" << num_shocks << endl;
      col=1;                                       // Houve pelo menos uma colisão
      cpar1=i;
      cpar2=j;
      //cout << " HOUVE COLISÃO \n CROSSING  entre posições : " <<  cpar1 << " e " << cpar2 ;

      t_c= dtt*(part.x[cpar2]-part.x[cpar1])/(part.x[cpar2]-part.x[cpar1]+npart.x[cpar1]-npart.x[cpar2]);    // Calculating tc1
      //cout << " \n ----- 1a APROXIMAÇÃO AO TEMPO DE CROSSING --- TC1 = " << t_c << endl; 

      temp=part.x[cpar1]+part.vx[cpar1]*sin(Wp*t_c)-X[cpar1]*(1-cos(Wp*t_c));                        // temporary values for tc2
      temp1=part.x[cpar2]+part.vx[cpar2]*sin(Wp*t_c)-X[cpar2]*(1-cos(Wp*t_c));  
      
      if(part.x[cpar2]-part.x[cpar1]+temp-temp1 > 0.00001){                                         // Guarantee code doesn't explode

        t_c2= (t_c)*(part.x[cpar2]-part.x[cpar1])/(part.x[cpar2]-part.x[cpar1]+temp-temp1);              // Calculating tc2
    }
    else t_c2=t_c;

      //if ( t_c2 < dt && t_c2>0 ){
      col_pos.push_back(cpar1);                  // Storing tc2 values and positions of colliding particles
      vec_cross.push_back(t_c2);
      //cout << " - SEGUNDA APROXIMAÇÃO AO TEMPO DE CROSSING -- TC2 = " << t_c2  << endl;
      //}
    }
  }

  cross_size=vec_cross.size();

  if(col!=0){            /// Algorithm to get the smallest value out of vec_cross

    min_tc2=vec_cross[0];
    minp=0;

    for (int j = 1; j < cross_size; ++j){

      if(vec_cross[j]<min_tc2 ){ 

        min_tc2=vec_cross[j];
        minp=j;                      // Also get the position of the collision corresponding to tc2_min
      }
    }

    //cout << " >>>>>>>>>>> SELECIONEI A COLISÃO " << npart.num[col_pos[minp]] << " E " <<  npart.num[col_pos[minp]+1] << endl;
    //cout << " ->>>--------- TEMPO DE CROSSING FINAL --------- TC2 = " << min_tc2 << endl;
    
    store_time=store_time+min_tc2;                                 // storing time to know where i am at
    //cout << " TEEEMPO INCREMENTADO  " << store_time << endl;

    loop(min_tc2,col_pos[minp]);                                     // Loop to advance particles positons up to tc2
  }

  part.x = npart.x;                                                 // storing new advanced values in part.
  part.vx= npart.vx;

  cout << " TEMPO ACTUAL DA SIMULAÇÂO " << t + store_time << endl; 
  //t=t+min_tc2;

  if ( dt - store_time < 0.0000)                          // If i haven't reached t+ dt yet i iterate the remaining time and look for collisions once more
  {
    //cout << " dt - store_time " << dt- store_time<< endl;
    loop(dt-store_time,-1);
    dtt=dt-store_time;
    goto LOOP;
  } 
  
  //cout << " TEMPO ACTUAL DA SIMULAÇÂO " << t + store_time << endl;

  return min_tc2;
}

void energy( double time ){

  double v2, kinetic, potential, pot, etotal;
  double xij, xij2;
  
  etotal=0.0;
  kinetic=0.0;
   potential = 0.0;
  for (int i = 0; i < NPart; ++i){
// kinectic energy of the system 
    v2=part.vx[i]*part.vx[i];
    kinetic = kinetic + 0.5*v2;
    
  // Potential Energy of the system  
    xij = part.x[i] - pos[i];
    xij2 = xij*xij;
    pot = +  ( xij2/2 );
    potential = potential + pot;
  }



  // Total energy of the system
  etotal = kinetic + potential;

  // Print the energies and current timestep
  //cout << " \t " << time << "    "  << kinetic << "    "  << potential << "    "  << etotal  << endl;

  E_kin=kinetic;
  E_pot=potential;
  E_tot=etotal;
  
  // na primeira iteração grava os valores das energias de modo a comparar a sua conservação 
  if(t==tmin){
	  E_kin_0=kinetic;
  E_pot_0=potential;
  E_tot_0=etotal;
	  }
}


void velocity_statistics(double time){

double moment1=0.0, moment2=0.0, moment3=0.0, moment4=0.0; 	
	
double mean; 
double desvpad; 
double skew;	
double kurt;	

 for (int i = 0; i< NPart; i++){ 
	 moment1 += part.vx[i];
	 moment2 += pow(part.vx[i],2);
	 moment3 += pow(part.vx[i],3);
	 moment4 += pow(part.vx[i],4);
	 }

moment1 = moment1/NPart; 
moment2 = moment2/NPart; 
moment3 = moment3/NPart; 
moment4 = moment4/NPart; 

mean = moment1;
desvpad = sqrt(moment2 - moment1*moment1); 
skew = moment3/pow(desvpad,3);
kurt = moment4/pow(desvpad,4)-3;
	
	estatisticas_vel << time << "\t" << mean <<"\t" << desvpad <<"\t" << skew <<"\t" << kurt << endl;
	
}
