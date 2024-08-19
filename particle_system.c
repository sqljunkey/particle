#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<string.h>

#define num_particles 500
#define animation_size 10000
#define n_1 (double)(num_particles+100)
#define momentum_multiplier  1.0
#define potential_multiplier 1.0
#define decay_ratio 1.0

typedef struct{
  double x;
  double y;
  double z;
  double x_momentum;
  double y_momentum;
  double z_momentum;
  //double weight; 
  double mass;
  char *color;
}particle;





typedef struct{

  particle *particles;
   

}frame;





frame    frames[animation_size];
particle particles[num_particles];
particle global_particle;
particle shifted_particle;


double rand_zero_to_one() {
    return rand() / (RAND_MAX + 1.0);
}

char* random_color() {
       
    int random_num = rand() % 3;

   
    switch (random_num) {
        case 0:
            return "GREEN";
        case 1:
            return "BLUE";
        case 2:
            return "RED";
        default:
            return "UNKNOWN"; 
    }
}


double truncated_normal(double mean, double std) {
    
    const double pi = 3.14159265358979323846;
    double x, y, z;

    do {
       
        x = 2.0 * rand_zero_to_one() - 1.0;
        y = 2.0 * rand_zero_to_one() - 1.0;
        z = x * x + y * y;
    } while (z >= 1.0 || z == 0.0);

   
    return mean + std * x * sqrt(-2.0 * log(z) / z);
}


void init_particles(){

    for(int i=0; i < num_particles; i++){

   
    particles[i].x = truncated_normal(0.0, 100.0);
    particles[i].y = truncated_normal(0.0, 100.0);
    particles[i].z = truncated_normal(0.0, 100.0);
    
    particles[i].x_momentum =0.0;
    particles[i].y_momentum =0.0;
    particles[i].z_momentum =0.0;

    particles[i].mass = 1.0;
    particles[i].color = random_color();

    
  }



}


void reset_global(){

  global_particle.x =0.0;
  global_particle.y =0.0;
  global_particle.z =0.0;

  global_particle.x_momentum = 0.0;
  global_particle.y_momentum = 0.0;
  global_particle.z_momentum = 0.0;

  global_particle.mass = 0.0;

}


void calculate_global(){

  reset_global();


  /*
  for(int i=0;i< num_particles;i++){

    global_particle.x +=particles[i].x;
    global_particle.y +=particles[i].y;
    global_particle.z +=particles[i].z;
    global_particle.mass += particles[i].mass;

    }*/


 
}

void subtract_from_global( const particle a){

  double scaling = 1.0;

  shifted_particle.x = (global_particle.x - a.x*scaling)/n_1;
  shifted_particle.y = (global_particle.y - a.y*scaling)/n_1;
  shifted_particle.z = (global_particle.z - a.z*scaling)/n_1;
  shifted_particle.mass = (global_particle.mass-a.mass);
 
  
}

double newton_potential(double distance) {
  double newton = 0.0;

    if (distance > 10.0) {

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

void update_particle_location(particle* p) {
    p->x = p->x + p->x_momentum * momentum_multiplier;
    p->y = p->y + p->y_momentum * momentum_multiplier;
    p->z = p->z + p->z_momentum * momentum_multiplier;
}

void update_particle_momentum(particle* p, double* vector, double potential) {
  
    p->x_momentum = p->x_momentum * decay_ratio + (vector[0] * potential_multiplier * potential);
    p->y_momentum = p->y_momentum * decay_ratio + (vector[1] * potential_multiplier * potential);
    p->z_momentum = p->z_momentum * decay_ratio + (vector[2] * potential_multiplier * potential);

   
}

void print_particle_locations(size_t step) {
  printf("Particle locations: %d\n", step);
    
    for (size_t i = 0; i < num_particles; i++) {
      
        printf("Particle %zu: x = %.6f, y = %.6f, z = %.6f\n"
	       , i + 1, particles[i].x, particles[i].y, particles[i].z);
    }

    printf("\n");
}
void calculate_next_step_naive(){

    
  for(int i=0;i< num_particles; i++){

    for(int j =0; j< num_particles; j++ ){


      if(j!=i){
    
      double potential = lennard_jones_potential( get_distance(particles[j], particles[i]));
      //double potential = newton_potential( get_distance(shifted_particle, particles[i]));
      double *vector = get_vector(particles[i] , particles[j]);

     

    update_particle_momentum(&particles[i], vector, potential);  
    update_particle_location(&particles[i]);
    free(vector);
    }
    }
  
  
  }



}
void calculate_next_step_global(){



  calculate_global();

  
  for(int i=0;i< num_particles; i++){

    subtract_from_global(particles[i]);

    
    double potential = lennard_jones_potential( get_distance(shifted_particle, particles[i]));
    //double potential = newton_potential( get_distance(shifted_particle, particles[i])) * shifted_particle.mass;
    double *vector = get_vector(particles[i] , shifted_particle);

     

    update_particle_momentum(&particles[i], vector, potential);  
    update_particle_location(&particles[i]);
    free(vector);

  
  
  }



}

void copy_particles_to_frames(size_t position) {
    if (position < animation_size) {
        
        frames[position].particles = malloc(num_particles * sizeof(particle));
        if (frames[position].particles == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            exit(EXIT_FAILURE); 
        }

       
        memcpy(frames[position].particles, particles, num_particles * sizeof(particle));
    }
}



void write_particle_positions(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return;
    }

   
    
    for (int i = 0; i < animation_size; i++) {
      
      fprintf(file, "%d\n\n", num_particles);
      
      for(int j =0; j <num_particles; j++ ){

	
       fprintf(file, "%s %.6f %.6f %.6f\n", frames[i].particles[j].color, frames[i].particles[j].x, frames[i].particles[j].y, frames[i].particles[j].z);
       
       
      }

    }

    fclose(file);

    
    for (int i = 0; i < animation_size; i++) {
        
         
       free(frames[i].particles);
    
    }
    
}


int main(int argc, char *argv[] ){
  /*if(){
    printf("Usage: sys num_frames num_particles");
    }*/
  srand(time(NULL));

  init_particles();
  calculate_global();

 struct timeval stop, start;
 gettimeofday(&start, NULL);


 for(int i=0; i < animation_size; i++){

    calculate_next_step_naive();
    copy_particles_to_frames(i);

  }
  gettimeofday(&stop, NULL);
  printf("Process took %lu u-seconds.\n", (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec);

  write_particle_positions("file.xyz");

  

  
  return 0;

}
