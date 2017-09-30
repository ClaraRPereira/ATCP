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
   particles part;                  // Initializing one set of particles
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
  vector<double> old_x;
  vector<double> old_vx;
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
    
/*
  for ( i = 0; i < NPart; ++i) // Applying equations (2) and (3) from Reference [2] 
    {
    	pos[i]=pos[i]*cos(Wp*t)+ vel[i]*sin(Wp*t)/Wp;
    	cout << " POS: " << pos[i] << endl;
    	vel[i]=vel[i]*cos(Wp*t)- pos[i]*sin(Wp*t)*Wp;
    	cout << " VEL: " << vel[i] << endl;
    }  
 */

  for ( i=0 ; i<n; i++ )
    {
      old_x.push_back(pos[i]) ;  // Guardar valores da posiçao e velocidade no tempo t para quando recalcularmos
      old_vx.push_back(vel[i]);  // tudo com o novo timestep (usando o tc1 e mais tarde o tc2 em vez do delta t inicial)
    }

    part.set_values(pos,vel,mass);
}

void func( )
{
  
  
 for (int i = 0; i < NPart; ++i)
   {
   	/* code */
   }  

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
      for ( j=i+1 ; j<NPart ; j++  )       {
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
    	
   // Number of particles
   NPart = 5;
   int n=NPart;
   // Max velocity
   Vt=4;

   // Time parameters
   tmin = 0.0;
   tmax = 10.0;
   dt = 0.5;
   t = tmin; // Initial time

   int print=0; // Variable to decide if i want to pront stuff

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
      // Compute positions and velocities at current timestep
      //func();
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
