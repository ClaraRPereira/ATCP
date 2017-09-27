#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>


using namespace std;


class particles
{
public:
  void set_values (vector<double>,vector<double>,vector<double>);
  vector<double> dir_product (vector<double>,vector<double>);
  vector<double> position;
  vector<double> velocity;
  vector<double> momentum;
};
   
   int NPart;                       // No of particles
   particles part;                  // Initializing one set of particles
   double dt, tmin, tmax;           // Time grid
   vector<double> mass(NPart);      // Vector for storing particle masses
   vector<double> pos(NPart);       // Positions
   vector<double> vel(NPart);       // Velocities
   FILE *fp_traj;                   // Create a file to store the trajectories

vector<double> particles::dir_product (vector<double> a, vector<double> b) // Defining direct product
{
  vector<double> product(NPart);

  for(int i=0; i<NPart; i++)
  {
    product[i]=a[i]*b[i]; // Multiplying component by component
  }
  return product;
} 


void particles::set_values (vector<double> x, vector<double> vx, vector<double> m) {  // Setting positions, velocities and momenta
  position = x;
  velocity = vx;
  momentum = part.dir_product(m,velocity);
}



int main()
{
	cout << "\n \t  ****** 1D PLASMA MODEL ****** \n" << endl;
    	
   // Number of particles
   NPart = 5;

   // Time parameters
   tmin = 0.0;
   tmax = 10.0;
   dt = 0.001;

   //mass.reserve(NPart);

   
   for (int i = 0; i < NPart; ++i)
   {
     mass.assign(NPart,4);
     pos.push_back(3);
     vel.push_back(2);
   }

    //cout << "SIZE: " << mass.size() << endl;
   //part.set_values()
	part.set_values(pos,vel,mass);

    cout << " Nº Partículas:  " << NPart << " \n " << endl;

    for (int i = 0; i < part.position.size(); ++i)
    {
    	cout << " Partícula " << i << " \t Position : " << part.position[i] << " |  Velocity :  " ;
    	cout << part.velocity[i] << " |  Momentum : " << part.momentum[i] << " |" << endl;
    }

    cout << "\n COrreu tudo bem \n" << endl;
	
	return 0;
}
