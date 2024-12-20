#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<string.h>
#inlcude "random_particle.h"

typedef struct{
 
  double x;
  double y;
  double z;
  double momentum_x;
  double momentum_y;
  double momentum_z;
     
}particle;





particle * particles =NULL;




void create_particles(int num_particles){


  if(particles!=NULL){
    free(particles);
  }

  particles = malloc(sizeof(particle) * num_particles);
  if(particles==NULL){

     fprintf(stderr, "Memory allocation failed\n");
      exit(EXIT_FAILURE); 
  }

}

void init_particles(int num_particles, int patches){

    for(int i=0; i < num_particles; i++){

    particles[i].id = i+3;
    particles[i].x =0.0+ truncated_normal(0.0, 100.0);
    particles[i].y =0.0+ truncated_normal(0.0, 100.0);
    particles[i].z =0.0+ truncated_normal(0.0, 100.0);
    
    particles[i].mass = 1.0;
   

    
  }


}



double newton_potential(double distance, double cutoff) {
  double newton = 0.0;

    if (distance > cutoff) {

	newton = 1.0 / pow(distance,2.0);

	}


  return newton;
}


double lennard_jones_potential(double distance) {

		double lj_pot_1 = 0.0;
		double lj_pot_2 = 0.0;
		double diff = 1E1;

		double epsilon = 1.1;
		double sigma = 10.0;

		lj_pot_1 = epsilon * (pow((sigma / (distance - diff)), 2.0) - 2.0 * pow((sigma / (distance - diff)), 1.0));
		lj_pot_2 = epsilon * (pow((sigma / (distance)), 2.0) - 2.0 * pow((sigma / (distance)), 1.0));

		if(distance >30.0) {
		  
		        return (lj_pot_2 - lj_pot_1);

		} else {
			return 0.0;
		}
}


double get_distance(particle a, particle b) {
		

     double distance_x = -a.x + (b.x);
     double distance_y = -a.y + (b.y);
     double distance_z = -a.z + (b.z);

     double distance = sqrt( (
			        pow(distance_x, 2.0)
			      + pow(distance_y, 2.0)
			      + pow(distance_z, 2.0)));

     
     return distance; 
		
}

double* get_vector(particle a, particle b) {
    
    double* vector = (double*)malloc(3 * sizeof(double));
    if (vector == NULL) {
        
        return NULL;
    }

   
    double dX = -(a.x) + (b.x);
    double dY = -(a.y) + (b.y);
    double dZ = -(a.z) + (b.z);

    
    vector[0] = dX;
    vector[1] = dY;
    vector[2] = dZ;

    return vector;
}


void update_particle_location(particle* p , momentum * m, double speed) {

    p->x = p->x + m->x * speed;
    p->y = p->y + m->y * speed;
    p->z = p->z + m->z * speed;

}


void update_particle_momentum( momentum*  m,double* vector, double potential, double decay_ratio, double speed) {
  
  m->x = m->x * decay_ratio + (vector[0] * speed * potential);
  m->y = m->y * decay_ratio + (vector[1] * speed * potential);
  m->z = m->z * decay_ratio + (vector[2] * speed * potential);

   
}


int main(){

  init_particles();

}
