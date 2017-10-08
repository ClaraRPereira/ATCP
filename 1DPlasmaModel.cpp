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

# define Wp 1.0

using namespace std;



// Create files to store the particle data,  energies et caetera
ofstream data;
ofstream energies;
ofstream energies_normalized;
ofstream fin_vel;
ofstream ini_vel;
ofstream vel1;
ofstream estatisticas_vel;
ofstream ergodic;
ofstream liouville;
////////////////////////////////////////////////////////////////////////

// DECLARAÇÃO FUNÇÕES
void initial_conditions(int NPart,double L,double Vt);        // Inicia as partículas com as condições iniciais
void time_step_iteration(int NPart,double dtime , int k);     // Methods  for the iteration 
double Dawson_algorithm(int NPart , double dt);				  // Method  for the main algorithm
void record_trajectories(int NPart ,double t, double dt);     // Writes trajectories to a file
void write_velocities(int NPart, int i);					  // Writes the velocities to a file 
void energy(int NPart, double time, double& E_kin,double& E_pot,double& E_tot, double& E_kin_0,double& E_pot_0,double& E_tot_0); //calculates the energie 
void record_energies( double time, double& E_kin,double& E_pot,double& E_tot, double& E_kin_0,double& E_pot_0,double& E_tot_0);  //write the energies to a file 
void record_energies_normalized( double time, double& E_kin,double& E_pot,double& E_tot, double& E_kin_0,double& E_pot_0,double& E_tot_0);//write the normalized energies to a file 
void ergodic_test(int NPart, double t  ,double dt, double&  average);  // compares the time average of v^2 of the central particle with the average v^2 over all  particles 
void velocity_statistics(int NPart,double time);    // calcullates the mean , standard deviation , skewness , and kurtosis of velocity distribution 
////////////////////////////////////////////////////////////////////////

// Variáveis iniciadas no particles.cpp
extern particles part;      // Initializing one set of particles at time t
extern particles npart;     // Set new set of particles at time t+dt 
extern int NPart;           // Number of particles 
vector<double> grid;     	// vector for storing equilibrium positions on the x axis
vector<double> mass;      	// Vector for storing particle masses   
vector<double> pos;       	// Positions
vector<double> vel;       	// Velocities
vector <double> X;        	// Displacements

//// Condições para teste de choques inelasticos 
double dE=0.8;
int Inelastic=0;         // If Inelastic==0 --> fazemos colisões elásticas . Else if Inelastic == 1 fazemos colisões inelásticas.

///// opçoes de output 
bool colision = false;    // Checks if a colision was found
bool print    = true;     // Variable to decide if i want to print stuff


////////////////////////////////////////////////////////////////////////
//					                                                  //
//                       MAIN FUNCTION                                //
//                                                                    //
////////////////////////////////////////////////////////////////////////
int main(){

	srand((int) time(0));  // Seeding the random distribution 

	cout << "\n \t  ****** 1D PLASMA MODEL ****** \n" << endl;

	int NPart = 2500;   // Number of particles
	double dx = 0.05;
	double L=dx*NPart;
	double Vt = 1;        // Max absolute velocity 

	// Time parameters
	double   tmin = 0.0;
	double tmax = 300;   // Simulation time
	double dt = 0.001;    // Time step
	double t = tmin;     // Initial time


	//variáveis para imprimir as energias
	double E_kin , E_pot , E_tot;
	double E_kin_0 , E_pot_0 , E_tot_0;			        

	double average_V2=0.0 ;

	// Create the files to store the data
	data.open ("raw_data.dat");
	energies.open ("energy.dat");
	energies_normalized.open ("energy_normalized.dat");
	fin_vel.open ("final_velocity.dat");
	ini_vel.open ("initial_velocity.dat");
	vel1.open ("vel1.dat");
	estatisticas_vel.open("stats.dat");
	ergodic.open("ergodic.dat");
	liouville.open("liouville.dat");

	// Call Dawson_algorithmtion to initialize the particles' variables (position, velocity and momenta)
	initial_conditions(NPart,L,Vt);


	double op=0;
	// Dynamics Iteration
	while ( t < tmax )
	{
		cout << " TEMPO DA SIMULAÇÃO " << t << endl;
		op = op + dt;
		// Compute positions and velocities at current timestep and determine crossing positions
		Dawson_algorithm(NPart,dt);


		if (colision) 
		break; // PARAR O time_step_iteration SE JÁ ENCONTREI A COLISÃO

		// Go to the next timestep
		t = t + dt;


		if (t == dt)
		{
			write_velocities(NPart,1); // records the initial distribution of velocity 
		} 
		if(op>=1)
		{
			write_velocities(NPart,3); 
			op=0;
		}
		if (t >= tmax) write_velocities(NPart,2); //records the  final distribution of velocity 

		// Write the data files and perform some diagnosis 
		if (print)
		{
			energy(NPart, t,E_kin,E_pot,E_tot,E_kin_0,E_pot_0,E_tot_0);
			record_energies(t,E_kin,E_pot,E_tot,E_kin_0,E_pot_0,E_tot_0 );
			record_energies_normalized(t,E_kin,E_pot,E_tot,E_kin_0,E_pot_0,E_tot_0 );
			record_trajectories(NPart ,t, dt);
			ergodic_test(NPart,t,dt, average_V2);
			velocity_statistics(NPart,t);
		}
		
	}  //final of the while cicle of the algorithm 


	data.close();
	energies.close();
	energies_normalized.close();
	fin_vel.close();
	ini_vel.close();
	vel1.close();
	ergodic.close();
	estatisticas_vel.close();
	liouville.close();
	return 0;
}

// Gonna Define the initial Conditions
void initial_conditions(int NPart,double L,double Vt)
{
	vector<double> Px; 
	double m = 1;           // Let's start by defining all masses as 1
	double r1;              // auxiliary variable to generate a random
	double axis = 0;        // Start of x axis
	double spc  = L/NPart;  // Defining intersheet spacing
	grid.reserve(NPart);
	pos.reserve(NPart);
	vel.reserve(NPart);
	X.reserve(NPart);
	mass.reserve(NPart);
	part.num.reserve(NPart);
	npart.num.reserve(NPart);

	for (int i = 0 ; i< NPart ; i++ ) // For time_step_iteration to set the different initial vectors
	{         
		part.num.push_back(i);                   // Defining Particle ordering	
		npart.num.push_back(i);
		//Defining the x positions  
		grid.push_back(axis); 
		pos.push_back(axis);
		axis += spc;
		// random velocities according to uniform distribution
		r1 = (double)rand() / (double)RAND_MAX;  // generating rando between 0 and 1
		vel.push_back(Vt-2*Vt*r1);  // Randomly distributing velocities between -Vt and Vt
		// Initializing displacement vector at zero
		X.push_back(0);                          
		// Defining the masses
		mass.push_back(m);
	}
	part.set_values(pos,vel,mass);             // Setting the values for the current particles and the new particles. 
	npart.set_values(pos,vel,mass);
}

void time_step_iteration(int NPart, double dtime , int a )
{

	double aux;    // auxiliary variable to store old values of stuff
	double vel_max;   // Store relative values of velocities in case we are considering inelastic collisions
	double vel_min;
	double mean;   // AUXILIARY CHEATING VALUE TO CALCULATE MEAN VALUE BETWEEN PARTICLES POSITIONS WHEN THEY COLLIDE
				
	// defining the cosine and sine beforehand for better calculation time 
	double cosx=cos(Wp*dtime);
	double sinx=sin(Wp*dtime);
	// CALCULATING TIME EVOLUTION OF THE "HARMONIC OSCILLATOR" Equation. Solution to differential equation for each particle
	for (int i = 0; i < NPart; ++i)
	{	  
		X[i]=part.x[i]-pos[i];
		npart.vx[i]=part.vx[i]*cosx-Wp*X[i]*sinx;    // npart are the new values for the particles. part are the old values (previous time step).
		npart.x[i]=part.x[i]+part.vx[i]*sinx-X[i]*(1-cosx);
		//if (npart.x[i]>=L || npart.x[i]<=0) npart.vx[i]=0;        // condições fronteira de parede absorsora 
		//if (npart.x[i]>=L || npart.x[i]<=0) npart.vx[i]=-npart.vx[i];        // condiçoes fronteira de parede reflectora 
	}   
	if(a!=-1)    // if a=-1 means the new npart values will still be tested for collisions and are not meant to be stored yet.
	{ 			 // When a!=0, a assumes the value of the position of the colliding particle. 


		//  Velocities are switched with its neighbour (elastic collision)
		aux=npart.vx[a];
		npart.vx[a]=npart.vx[a+1];            
		npart.vx[a+1]=aux;

		// X positions of colliding particles are the exact spot where they collide at. (mean position at t=dt)
		mean=(npart.x[a+1]+npart.x[a])*0.5;                                        
		npart.x[a]=mean;                       
		npart.x[a+1]=mean;

	}
	if(a!=-1 && Inelastic==1)   // if a=-1 means the new npart values will still be tested for collisions and are not meant to be stored yet.
	{						 	// When a!=0, a assumes the value of the position of the colliding particle. 
		if(npart.vx[a]>npart.vx[a+1])
		{
			vel_max=npart.vx[a]; 
			vel_min=npart.vx[a+1];
			npart.vx[a]=npart.vx[a]-dE*(vel_max-vel_min);
			npart.vx[a+1]=npart.vx[a+1]+dE*(vel_max-vel_min);
		}
		else if(npart.vx[a]<npart.vx[a+1])
		{ 
			vel_max=npart.vx[a+1]; 
			vel_min=npart.vx[a];
			npart.vx[a+1]=npart.vx[a+1]-dE*(vel_max-vel_min);
			npart.vx[a]=npart.vx[a]+dE*(vel_max-vel_min);
		}
		else if ( npart.vx[a]==npart.vx[a+1]){  
		//  Velocities are switched with its neighbour (elastic collision)
		aux=npart.vx[a];
		npart.vx[a]=npart.vx[a+1];            
		npart.vx[a+1]=aux;
		}
		// X positions of colliding particles are the exact spot where they collide at. (mean position at t=dt)
		mean=(npart.x[a+1]+npart.x[a])*0.5;                                        
		npart.x[a]=mean;                      
		npart.x[a+1]=mean;
	}
}    

double Dawson_algorithm( int NPart , double dt )
{

	vector<double> col_pos;  // stores positions of colliding particles
	col_pos.reserve(NPart);
	int cpar1=0;             // position of 1st colliding particle
	int cpar2=0;             // position of second colliding particle 
	int n=NPart;             // storing total no. os particles
	int col=0;               // col=0 if there aren't any collisions. col=1 if there are collisions
	int j=0;                 // Iterator
	int num_shocks=0; // NUmber of collisions

	double dtt=dt;           // Time step
	double t_c=0;            // Variable for the crossing time, tc1 
	double t_c2;             // Variable for the crossing time, tc2 
	double min_tc2=0;        // Minimum tc2 (time of first crossing event)
	double minp;             // Position of particle colliding at min_tc2
	double store_time=0;     // Stores sum of time steps made so far
	double temp, temp1;      // Auxiliary 

	vector<double> vec_cross;   // Vector stores tc2 values
	double cross_size;          // size of vector vec_cross (just to avoid warnings and always compare integer values)
	vec_cross.reserve(NPart);

	time_step_iteration(NPart,dtt,-1);            // First time_step_iteration, calculates new npart positions and velocities for my time step dt

	time_step_iteration:
	col=0;
	vec_cross.clear();
	col_pos.clear();
	
	for (int i = 0; i < n-1; ++i) //ciclo de pesqueisa de choques 
	{
		j=i+1;
		if(npart.x[i]>=npart.x[j] )
		{
			num_shocks=num_shocks+1;
			col=1;                                       // Houve pelo menos uma colisão
			cpar1=i;
			cpar2=j;
			// Calculating tc1
			t_c= dtt*(part.x[cpar2]-part.x[cpar1])/(part.x[cpar2]-part.x[cpar1]+npart.x[cpar1]-npart.x[cpar2]);    
			double cosx=cos(Wp*t_c);
			double sinx=sin(Wp*t_c);
			temp=part.x[cpar1]+part.vx[cpar1]*sinx-X[cpar1]*(1-cosx);                        // temporary values for tc2
			temp1=part.x[cpar2]+part.vx[cpar2]*sinx-X[cpar2]*(1-cosx);  
			if(part.x[cpar2]-part.x[cpar1]+temp-temp1 > 0.00001)  // Guarantee code doesn't explode
			{        // Calculating tc2                                     
				t_c2= (t_c)*(part.x[cpar2]-part.x[cpar1])/(part.x[cpar2]-part.x[cpar1]+temp-temp1);         
			}
			else t_c2=t_c;
			col_pos.push_back(cpar1);                  // Storing tc2 values and positions of colliding particles
			vec_cross.push_back(t_c2);
		}
	}

	cross_size=vec_cross.size();

	if(col!=0)           //Algorithm to get the smallest value out of vec_cross
	{  

		min_tc2=vec_cross[0];
		minp=0;

		for (int j = 1; j < cross_size; ++j)
		{
			if(vec_cross[j]<min_tc2 )
			{ 
				min_tc2=vec_cross[j];
				minp=j;                      // Also get the position of the collision corresponding to tc2_min
			}
		}

		// storing time to know where i am at
		store_time=store_time+min_tc2;                                
		// time_step_iteration to advance particles positons up to tc2
		time_step_iteration(NPart , min_tc2,col_pos[minp]);                                     
	}
	// storing new advanced values in part.
	part.x = npart.x;                                                 
	part.vx= npart.vx;

	// If i haven't reached t+ dt yet i iterate the remaining time and look for collisions once more
	if ( dt - store_time < 0.0000)                          
	{
		time_step_iteration(NPart,dt-store_time,-1);
		dtt=dt-store_time;
		goto time_step_iteration;
	} 
  return min_tc2;
}

void energy(int NPart, double time, double& E_kin,double& E_pot,double& E_tot, double& E_kin_0,double& E_pot_0,double& E_tot_0){

	double v2, kinetic, potential, pot, etotal;
	double xij, xij2;
  
	etotal=0.0;
	kinetic=0.0;
	potential = 0.0;
	for (int i = 0; i < NPart; ++i)
	{
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

	// Saves the energies 
	E_kin=kinetic;
	E_pot=potential;
	E_tot=etotal;
  
	// na primeira iteração grava os valores das energias de modo a comparar a sua conservação 
	if(time==0.001)
	{
		E_kin_0=kinetic;
		E_pot_0=potential;
		E_tot_0=etotal;
	}
}

   
    
void  ergodic_test(int NPart,double time, double deltatime, double& average) // pretende comparar a media temporal com a media nas particulas
{
	double mean = 0.0; // media no conjunto
	double density=0.0; // densidade de particulas no espaço d fase
	// avança a media temporal do quadrado da velocidade da particula central 
	average += part.vx[NPart/2]*part.vx[NPart/2]; 
	// calcula para esse mesmo tempo a média nas particulas da mesma quantidade
	for(int i = 0 ; i < NPart ; ++i)
	{
		mean += part.vx[i]*part.vx[i]/NPart;
		// conta quantas particulas estao numa regiao do espaço de fases 0<X<0.0 e 0<v<0.1
		if( X[i]>=0.0 && X[i]<=0.1 && part.vx[i]>=0.0 && part.vx[i]<=0.1  )
		{
		density += 1.0;	
		}
		
	}
	//escreve os dados 
	liouville << time <<"\t"<< density/0.01 <<endl; 
	ergodic << time <<"\t"<<deltatime*average/time<<"\t"<<mean<<endl;
}


void velocity_statistics(int NPart,double time)
{
	double moment1=0.0, moment2=0.0, moment3=0.0, moment4=0.0; 	
	double mean , desvpad , skew , kurt;	
	//calculos dos momentos de probabilidade para todas as particulas
	for (int i = 0; i< NPart; i++)
	{ 
		moment1 += part.vx[i];
		moment2 += pow(part.vx[i],2);
		moment3 += pow(part.vx[i],3);
		moment4 += pow(part.vx[i],4);
	}
	moment1 = moment1/NPart; 
	moment2 = moment2/NPart; 
	moment3 = moment3/NPart; 
	moment4 = moment4/NPart; 
	//normalização dos momentos
	mean = moment1;
	desvpad = sqrt(moment2 - moment1*moment1); 
	skew = moment3/pow(desvpad,3);
	kurt = moment4/pow(desvpad,4)-3;
	
	estatisticas_vel << time << "\t" << mean <<"\t" << desvpad <<"\t" << skew <<"\t" << kurt << endl;
	
}

    
    


// Write the positions, velocities and momenta of the particles in the file " DATA ".
void record_trajectories(int NPart ,double t, double dt)
{
	for (int i=0 ; i<NPart ; i++ )
	{
		// time ---- position ---- displacement ---- velocity ---- particle id  
		data << t <<"\t"<< part.x[i] <<"\t"<< X[i] <<"\t"<< part.vx[i] <<"\t"<< part.num[i] <<endl; 
	}	
	data <<"\n"<< endl;
}

// Write the energies of the system in the file " ENERGY ".
void record_energies( double time, double& E_kin,double& E_pot,double& E_tot, double& E_kin_0,double& E_pot_0,double& E_tot_0)
{
	energies << time <<"\t"<< E_kin <<"\t"<< E_pot <<"\t"<< E_tot << endl; 
}


void record_energies_normalized( double time, double& E_kin,double& E_pot,double& E_tot, double& E_kin_0,double& E_pot_0,double& E_tot_0)
{
	energies_normalized << time <<"\t"<< (E_kin-E_kin_0)/E_kin_0 <<"\t"<< (E_pot-E_pot_0)/E_pot_0 <<"\t"<< (E_tot-E_tot_0)/E_tot_0 << endl;
}



void write_velocities(int NPart,int j)
{
	for (int i=0 ; i<NPart ; i++ )
	{
		if(j==2) fin_vel <<  part.vx[i] << endl;
		else if(j==1) ini_vel <<  part.vx[i] << endl;
		else if(j==3)  vel1 <<   part.vx[i] << endl;
	}
}    
