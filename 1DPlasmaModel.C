#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <ctime>

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

using namespace std;


class particles
{
public:
  void set_values (vector<double>,vector<double>,vector<double>);
  vector<double> dir_product (vector<double>,vector<double>);
  vector<double> x;
  vector<double> vx;
  vector<double> mx;
};
   
   int NPart;                       // No of particles
   particles part;                  // Initializing one set of particles at time t
   particles npart;                 // Set new set of particles at time t+dt 
   vector<double> grid;             // vector for storing equilibrium positions on the x axis
   double dt, tmin, tmax;           // Time grid
   double t;                        // Time variable (evolving time)  
   vector<double> mass(NPart);      // Vector for storing particle masses
   vector<double> pos(NPart);       // Positions
   vector<double> vel(NPart);       // Velocities
   FILE *data, *energies;                      // Create files to store the particle data and the energies
   double Vt;                       // Max velocity for uniform distribution   
   double E_kin;			        //variáveis para imprimir as energias
   double E_pot;
   double E_tot;
   double L=4;  //size of box

vector<double> particles::dir_product (vector<double> a, vector<double> b) // Defining direct product
{
  vector<double> product(NPart);

  for(int i=0; i<NPart; i++)
  {
    product[i]=a[i]*b[i]; // Multiplying component by component
  }
  return product;
} 

// Write the positions, velocities and momenta of the particles in the file " DATA ".
void record_trajectories( )
{
  int n=NPart;

  for (int i=0 ; i<n ; i++ )
    {
      fprintf( data, "%f %f %f %f\n", t, part.x[i], part.vx[i], part.mx[i] ); // Escrever os valores para o ficheiro "DATA"
    }
}

// Write the energies of the system in the file " ENERGY ".
void record_energies()
{

	fprintf( energies, "%f, %f, %f, %f\n",t, E_kin, E_tot, E_pot);

}

void particles::set_values (vector<double> position, vector<double> velocity, vector<double> mass) {  // Setting positions, velocities and momenta
   x = position;
  vx = velocity;
  mx = part.dir_product(mass,velocity);
}

void initial_conditions() // GOnna Define the initial COnditions
{
  vector<double> Px;
  

  double m=1; //Let's start by defining all masses as 1
  //double sigma=0.5, n0=0.7; //Valores aleatórios para a carga por unidade de área, sigma, e para a density of neutralizing background charges 
  //double Wp= 4*M_PI*sigma*sigma*n0/m; //Plasma Frequency
  double r1;     // auxiliary variable to generate a random
  int i;
  double axis=0; // Start of x axis
  double spc=1; //defining intersheet spacing
  int n=NPart;
  
  for ( i=0 ; i<n ; i++ ) 
    {

   //Defining the x positions 	
      //part.x[i]+=spc;
      grid.push_back(axis);	
      pos.push_back(axis);
      //cout << " pos: " << pos[i] << endl;
      axis+=spc;

   // random velocities according to uniform distribution
      r1 = (double)rand()/(double)RAND_MAX;  // generating rando between 0 and 1
      vel.push_back(Vt*r1);
     // cout << " vel: " << vel[i] << endl;

    // Defining the masses
    	mass.push_back(m);
    }

    part.set_values(pos,vel,mass);
}

int func( double t)
{
	
  int a=0; // Store position of crossing particle
  int k=0;
  int n=NPart;
  double m=1; //Let's start by defining all masses as 1
  //double sigma=0.5, n0=0.7; //Valores aleatórios para a carga por unidade de área, sigma, e para a density of neutralizing background charges 
  //double Wp= 4*M_PI*sigma*sigma*n0/m; //Plasma Frequency
  double Wp=1;
  vector <double> X;
  npart.x.reserve(n);
  npart.vx.reserve(n);
   double Delta_c;// Variables for the crossing times, tc1 and tc2

 

  for (int i = 0; i < n; ++i)
  {
  	X.push_back(part.x[i]-grid[i]);
  }

  // IF THERE AREN'T ANY CROSSINGS JUST NORMALLY CALCULATING TIME EVOLUTION OF "HARMONIC OSCILLATORS"
  
 LOOP:for (int i = 0; i < n; ++i)
   {
   npart.vx[i]=part.vx[i]*cos(Wp*t)-Wp*X[i]*sin(Wp*t);
   npart.x[i]=part.x[i]+part.vx[i]*sin(Wp*t)-X[i]*(1-cos(Wp*t));
    
   }  

 // LOOPS TO LOOK FOR CROSSINGS
for ( int i=0 ; i<n ; i++ )
    {
      for ( int  j=i+1 ; j<n ; j++  )      
       {
           if(npart.x[i]>npart.x[j])      { a=i; goto OUT;}  // Particle i collides with particle j . Tou a usar o goto só para
                                                              // conseguir sair simultaneamente dos dois loops.
           else a=-1;                      // There are no collisions
	   }
    }
  
  if (a==-1) // if there are no crossings we just store old values and refresh the new ones
   {
   	for (int i = 0; i < n; ++i)
   {
   	 pos[i]=part.x[i];
   	 vel[i]=part.vx[i];
     part.x[i]=npart.x[i];
     part.vx[i]=npart.vx[i];
   } 
   }  
  OUT:
  if (a!=-1){
    k=k+1;
	Delta_c= dt*(part.x[a+1]-part.x[a])/(part.x[a+1]-part.x[a]+npart.x[a]-npart.x[a+1]);
    t=Delta_c;
    if (k==1) {goto LOOP;}
    else if (k==2){ cout << " ------------------ TEMPO DE CROSSING ------------- Delta_c2 = " << Delta_c << endl;}
    }


 return a; 	
}

void energy( double time )
{
  double kinetic, potential, pot, etotal;
 
  double sigma=0.5, n0=0.7; //Valores aleatórios para a carga por unidade de área, sigma, e para a density of neutralizing background charges 
  double E0= 4*M_PI*sigma*n0; // Amplitude campo Eléctrico
  double xij, xij2;
  int i, j;
  int n=NPart;

  // Kinetic energy of the system
  /*kinetic = 0.0;
  for ( i=0 ; i<NPart ; i++ )
    {
      kinetic = kinetic + 0.5 * part.mx[i]*part.vx[i];
    }*/
  for (i = 0; i < n; ++i)
  {
  	 kinetic = 0.5*part.dir_product(part.mx,part.vx)[i];
  }
 
  // Potential Energy of the system
  potential = 0.0;
  for ( i=0 ; i<n ; i++ )
    {
      for ( j=i+1 ; j<n ; j++  )       {
            xij = part.x[i] - part.x[j];
            
            xij2 = xij*xij;
            
            pot = - E0 * ( xij2/2 );
            potential = potential + pot;
	}
    }

  // Total energy of the system
  etotal = kinetic + potential;

  // Print the energies and current timestep
    cout << " " << time << "    "  << kinetic << "    "  << potential << "    "  << etotal  << endl;
 
 E_kin=kinetic;
 E_pot=potential;
 E_tot=etotal;


}

int main()
{
	srand((int) time(0));  // Seeding the random distribution 

   cout << "\n \t  ****** 1D PLASMA MODEL ****** \n" << endl;
    	
   
   NPart = 7; // Number of particles
   int n=NPart;
   int cross;   // position of crossing
   Vt=3; // Max velocity
   int k=0;

   // Time parameters
   tmin = 0.0;
   tmax = 10.0;
   dt = 0.5;
   t = tmin; // Initial time

   int print=1; // Variable to decide if i want to pront stuff

   // Create the files to store the data
   data = fopen( "DATA", "w" );
   energies = fopen( "ENERGY", "w" );


   //mass.reserve(NPart);

    initial_conditions(); // Call function to initialize the particles' variables (position, velocity and momenta)

     // Print, to the terminal, positions, velocities & momenta.
    cout << " Nº Partículas:  " << NPart << " \n " << endl;

    for (int i = 0; i < n; ++i)
    {
    	cout << " Partícula " << i << " \t Position : " << part.x[i] << " |  Velocity :  " ;
    	cout << part.vx[i] << "   |  Momentum : " << part.mx[i] << " |" << endl;
    }

    cout << " \n ENERGIAS" << endl;
    cout << " time   kinetic   potential   total " << endl;

    // Dynamics Iteration
  while ( t < tmax )
    {
         k=k+1;
    	cout << " RUN " << k ;
      // Compute positions and velocities at current timestep and determine crossing positions
      cross=func(dt);
      if (cross!=-1) {cout << " \t \t \t \t Partícula " << cross << " choca com partícula " << cross+1 << endl; break;}
      else cout << " \t \t \t \t Partículas não chocaram " << endl;
       
      // Go to the next timestep
      t = t + dt;
      // Write the positions of the particles in the file "DATA" if print_trajectory ==1.
      if (print == 1){
         record_trajectories( );
         record_energies( );}

      // Compute energies at current timestep
      energy(t);

    }


  // Close trajectories file
  fclose( data );
  fclose( energies );
	

	return 0;
}
