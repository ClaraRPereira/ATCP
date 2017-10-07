#include "particles.h"

//int NPart = 100; 		         // No of particles
particles part;                  // Initializing one set of particles at time t
particles npart;                 // Set new set of particles at time t+dt 

vector<double> particles::dir_product (vector<double> a, vector<double> b){ // Defining direct product

	vector<double> product((int) a.size());

	for(int i=0; i< (int) a.size(); i++){

    	product[i] = a[i]*b[i]; // Multiplying component by component
    }

    return product;
}

// Setting positions, velocities and momenta
void particles::set_values (vector<double> position, vector<double> velocity, vector<double> mass){ 
	
	x  = position;
	vx = velocity;
	mx = part.dir_product(mass,velocity);
}
