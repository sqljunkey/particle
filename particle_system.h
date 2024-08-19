#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<string.h>


typedef struct{
  double x;
  double y;
  double z;
  double x_momentum;
  double y_momentum;
  double z_momentum;
  double mass;
  char *color;
  int *neighbor_list;
  int neighbor_list_size; 
}particle;



particle *PARTICLES;
int NUM_PARTICLES;
int **NEIGHBOR_MATRIX;
int **NEIGHBOR_MATRIX_COPY;
double RATIO;
int NUM_FRAMES;






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
       
        x = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
        y = 2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
        z = x * x + y * y;
    } while (z >= 1.0 || z == 0.0);

   
    return mean + std * x * sqrt(-2.0 * log(z) / z);
}




void init_particles( ){

    for(int i=0; i < NUM_PARTICLES; i++){

   
    PARTICLES[i].x = truncated_normal(0.0, 100.0);
    PARTICLES[i].y = truncated_normal(0.0, 100.0);
    PARTICLES[i].z = truncated_normal(0.0, 100.0);
    
    PARTICLES[i].x_momentum =0.0;
    PARTICLES[i].y_momentum =0.0;
    PARTICLES[i].z_momentum =0.0;

    PARTICLES[i].mass = 1.0;
    PARTICLES[i].color = random_color();
    PARTICLES[i].neighbor_list = NULL;
    PARTICLES[i].neighbor_list_size = 0;

    
  }

    NEIGHBOR_MATRIX = malloc(NUM_PARTICLES * sizeof(int *));
    if (NEIGHBOR_MATRIX == NULL) {

        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
        NEIGHBOR_MATRIX[i] = malloc(NUM_PARTICLES * sizeof(int));
        if (NEIGHBOR_MATRIX[i] == NULL) {
           
            perror("Failed to allocate memory for NEIGHBOR_MATRIX[i]");
            exit(EXIT_FAILURE);
        }
        
       
        for (int j = 0; j < NUM_PARTICLES; j++) {
            NEIGHBOR_MATRIX[i][j] = 0;
        }
    }

    

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





void update_particle_location(particle* p, double multiplier) {
    p->x = p->x + p->x_momentum * RATIO;
    p->y = p->y + p->y_momentum * RATIO;
    p->z = p->z + p->z_momentum * RATIO;
}






void print_particle_locations(size_t step, particle *PARTICLES, int num_PARTICLES ) {
  printf("Particle locations: %d\n", step);
    
    for (size_t i = 0; i < num_PARTICLES; i++) {
      
        printf("Particle %zu: x = %.6f, y = %.6f, z = %.6f\n"
	       , i + 1, PARTICLES[i].x, PARTICLES[i].y, PARTICLES[i].z);
    }

    printf("\n");
}





void add_neighbor(particle *p, int i) {
   
    int new_size = p->neighbor_list_size + 1;
    int *new_list = malloc(new_size * sizeof(int));
    
    if (p->neighbor_list != NULL) {
      
        memcpy(new_list, p->neighbor_list, p->neighbor_list_size * sizeof(int));
        free(p->neighbor_list);
    }
    
   
    new_list[p->neighbor_list_size] = i;
    p->neighbor_list = new_list;
    p->neighbor_list_size = new_size;
}



void free_matrix_copy(){



  for(int i=0; i < NUM_PARTICLES; i++){

    free(NEIGHBOR_MATRIX_COPY[i]);

  }

  free(NEIGHBOR_MATRIX_COPY);
  
}

void create_matrix_copy(){

     NEIGHBOR_MATRIX_COPY= malloc(NUM_PARTICLES * sizeof(int *));
    if (NEIGHBOR_MATRIX_COPY == NULL) {

        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);
    }


        for (int i = 0; i < NUM_PARTICLES; i++) {
        NEIGHBOR_MATRIX_COPY[i] = malloc(NUM_PARTICLES * sizeof(int));
        if (NEIGHBOR_MATRIX_COPY[i] == NULL) {
           
            perror("Failed to allocate memory for NEIGHBOR_MATRIX[i]");
            exit(EXIT_FAILURE);
        }


	 memcpy(NEIGHBOR_MATRIX_COPY[i], NEIGHBOR_MATRIX[i], NUM_PARTICLES * sizeof(int));
       
      
    }


}


void create_naive_neighbor_list(){

  for(int i=0; i< NUM_PARTICLES; i++){

    
    for(int j=0; j< NUM_PARTICLES; j++){

      if(i!=j){

	add_neighbor(&PARTICLES[i],j);
	NEIGHBOR_MATRIX[i][j]=1;

      }
    }

  }

}






void create_distance_based_neighbor_list(double radius){

  for(int i=0; i< NUM_PARTICLES; i++){

    
    for(int j=0; j< NUM_PARTICLES; j++){

      if(i!=j){

	double distance = get_distance(PARTICLES[i], PARTICLES[j]);
	if(distance<=radius){
	  add_neighbor(&PARTICLES[i],j);
	  NEIGHBOR_MATRIX[i][j]=1;
	  
	}

      }
    }

  }

}






void newton_force(int j){

  for(int i = 0; i < PARTICLES[j].neighbor_list_size; i++ ){

     if(NEIGHBOR_MATRIX_COPY[j][i]==1){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[p_j.neighbor_list[i]];

    int loc = p_j.neighbor_list[loc]; 

    double distance  = get_distance(p_j, p_i);
    double potential = newton_potential(distance);
    double *vector   = get_vector(p_j, p_i);
    

    PARTICLES[loc].x_momentum = PARTICLES[i].x_momentum + (vector[0] * potential);
    PARTICLES[loc].y_momentum = PARTICLES[i].y_momentum + (vector[1] * potential);
    PARTICLES[loc].z_momentum = PARTICLES[i].z_momentum + (vector[2] * potential);

    PARTICLES[j].x_momentum = PARTICLES[j].x_momentum + -(vector[0] * potential);
    PARTICLES[j].y_momentum = PARTICLES[j].y_momentum + -(vector[1] * potential);
    PARTICLES[j].z_momentum = PARTICLES[j].z_momentum + -(vector[2] * potential);
    

    NEIGHBOR_MATRIX_COPY[i][j]=0;
    NEIGHBOR_MATRIX_COPY[j][i]=0;
    }
  }

}






void traverse(void (*f)(int)){

  for(int i=0; i < NUM_PARTICLES; i++){

    f(i);

  }

}




void write_particle_positions(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return;
    }

   
    
   
      
      fprintf(file, "%d\n\n", NUM_PARTICLES);
      
      for(int j =0; j <NUM_PARTICLES; j++ ){

	
       fprintf(file, "%s %.6f %.6f %.6f\n", PARTICLES[j].color, PARTICLES[j].x, PARTICLES[j].y, PARTICLES[j].z);
       
       
      }

    

    fclose(file);

    

    
}
