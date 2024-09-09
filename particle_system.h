#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<sys/time.h>
#include<string.h>
#include <omp.h>
#include <stdbool.h>

typedef struct{
  double x;
  double y;
  double z;
  double x_momentum;
  double y_momentum;
  double z_momentum;
  double x_momentum_bond;
  double y_momentum_bond;
  double z_momentum_bond; 
  double mass;
  char *color;
}particle;



particle *PARTICLES;
int NUM_PARTICLES;
bool **NEIGHBOR_MATRIX;
bool **NEIGHBOR_MATRIX_COPY;
bool **BOND_MATRIX;
bool **BOND_MATRIX_COPY;

double GLOBAL_RATIO_MOVE;
double GLOBAL_RATIO_MOMENTUM;
double BOND_RATIO_MOVE;
double BOND_RATIO_MOMENTUM;

int NUM_FRAMES;
double DECAY_RATIO;
double SPREAD;
int VECTOR_SIZE = 3;

void init_matrices(){

  
    /*  neighbor matrix */
   
    NEIGHBOR_MATRIX = malloc(NUM_PARTICLES * sizeof(bool *));
    if (NEIGHBOR_MATRIX == NULL) {
         
        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);
    }

 

    for (int i = 0; i < NUM_PARTICLES; i++) {
        NEIGHBOR_MATRIX[i] = malloc(NUM_PARTICLES * sizeof(bool ));
        if (NEIGHBOR_MATRIX[i] == NULL) {
           
            perror("Failed to allocate memory for NEIGHBOR_MATRIX[i]");
            exit(EXIT_FAILURE);
        }
        
       
        for (int j = 0; j < NUM_PARTICLES; j++) {
            NEIGHBOR_MATRIX[i][j] = false;
        }
    }


        NEIGHBOR_MATRIX_COPY = malloc(NUM_PARTICLES * sizeof(bool *));
    if (NEIGHBOR_MATRIX_COPY == NULL) {
         
        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);
    }

 

    for (int i = 0; i < NUM_PARTICLES; i++) {
        NEIGHBOR_MATRIX_COPY[i] = malloc(NUM_PARTICLES * sizeof(bool));
        if (NEIGHBOR_MATRIX_COPY[i] == NULL) {
           
            perror("Failed to allocate memory for NEIGHBOR_MATRIX[i]");
            exit(EXIT_FAILURE);
        }
        
       
        for (int j = 0; j < NUM_PARTICLES; j++) {
            NEIGHBOR_MATRIX_COPY[i][j] = false;
        }
    }


    /* bond matrix */
    
     BOND_MATRIX = malloc(NUM_PARTICLES * sizeof(bool *));
    if (BOND_MATRIX == NULL) {
         
        perror("Failed to allocate memory for BOND_MATRIX");
        exit(EXIT_FAILURE);
    }

 

    for (int i = 0; i < NUM_PARTICLES; i++) {
        BOND_MATRIX[i] = malloc(NUM_PARTICLES * sizeof(bool));
        if (BOND_MATRIX[i] == NULL) {
           
            perror("Failed to allocate memory for BOND_MATRIX[i]");
            exit(EXIT_FAILURE);
        }
        
       
        for (int j = 0; j < NUM_PARTICLES; j++) {
            BOND_MATRIX[i][j] = false;
        }
    }


    BOND_MATRIX_COPY = malloc(NUM_PARTICLES * sizeof(bool *));
    if (BOND_MATRIX_COPY == NULL) {
         
        perror("Failed to allocate memory for BOND_MATRIX");
        exit(EXIT_FAILURE);
    }

 

    for (int i = 0; i < NUM_PARTICLES; i++) {
        BOND_MATRIX_COPY[i] = malloc(NUM_PARTICLES * sizeof(bool));
        if (BOND_MATRIX_COPY[i] == NULL) {
           
            perror("Failed to allocate memory for BOND_MATRIX[i]");
            exit(EXIT_FAILURE);
        }
        
       
        for (int j = 0; j < NUM_PARTICLES; j++) {
	  BOND_MATRIX_COPY[i][j] = false;
        }
    }



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

 char* beach_ball_color(float x, float y, float z) {

    
    if (y < 10.0 && y > -10.0 ) { 
        return "BLUE";
    } else if (y>=10.0) {
        return "GREEN";
    } else if (y<=-10.0) {
        return "RED";
    } else {
        return "UNKNOWN";
    }
}

void change_colors(){

  for(int i = 0; i < NUM_PARTICLES; i++){
        PARTICLES[i].color = beach_ball_color( PARTICLES[i].x
					  ,PARTICLES[i].y
					  ,PARTICLES[i].z);

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


void init_debug_particles(){

  NUM_PARTICLES = 4;
  
    PARTICLES = malloc(NUM_PARTICLES * sizeof(particle));

    if(PARTICLES == NULL){
      
        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);

    }

    PARTICLES[0].x =1.0;
    PARTICLES[0].y =0.0;
    PARTICLES[0].z =1.0;


    PARTICLES[1].x =1.0;
    PARTICLES[1].y =1.0;
    PARTICLES[1].z =0.0;


    PARTICLES[2].x =0.0;
    PARTICLES[2].y =1.0;
    PARTICLES[2].z =1.0;

    PARTICLES[3].x =1.0;
    PARTICLES[3].y =1.0;
    PARTICLES[3].z =1.0;
    
    for(int i=0; i < NUM_PARTICLES; i++){
    PARTICLES[i].x_momentum =0.0;
    PARTICLES[i].y_momentum =0.0;
    PARTICLES[i].z_momentum =0.0;

    PARTICLES[i].x_momentum_bond = 0.0;
    PARTICLES[i].y_momentum_bond = 0.0;
    PARTICLES[i].z_momentum_bond = 0.0;

    PARTICLES[i].mass = 1.0;
    PARTICLES[i].color =random_color();
    }

    init_matrices();

}

void init_particles( ){


      int num_threads = omp_get_max_threads();
      printf("Setting number of threads to maximum: %d\n", num_threads);

      
      omp_set_num_threads(num_threads);
      double spread = 50.0;


    PARTICLES = malloc(NUM_PARTICLES * sizeof(particle));

    if(PARTICLES == NULL){
      
        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);

    }
    
    for(int i=0; i < NUM_PARTICLES; i++){

   
    PARTICLES[i].x = truncated_normal(0.0, SPREAD);
    PARTICLES[i].y = truncated_normal(0.0, SPREAD);
    PARTICLES[i].z = truncated_normal(0.0, SPREAD);

  
    
    PARTICLES[i].x_momentum =0.0;
    PARTICLES[i].y_momentum =0.0;
    PARTICLES[i].z_momentum =0.0;

    PARTICLES[i].x_momentum_bond = 0.0;
    PARTICLES[i].y_momentum_bond = 0.0;
    PARTICLES[i].z_momentum_bond = 0.0;

    PARTICLES[i].mass = 1.0;
    PARTICLES[i].color =random_color();


    
  }


    init_matrices();

    

}




double spring_potential(double distance){


  double result = 0.0;
  double equilibrium_distance = 50.0;
  double max_distance = 100000.0;


    
    if(distance < 1000000){
        result =distance - equilibrium_distance;
    }


		
    return result;
		


}

double newton_potential(double distance) {
  double newton = 0.0;


  if(distance > 50.0){

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
    
  double* vector = (double*)malloc(VECTOR_SIZE * sizeof(double));
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


double* get_normalized( double *vector){

  double power=0.0;

  for(int i=0; i < VECTOR_SIZE; i++){

    power += pow(vector[i], 2.0);
  }

  double magnitude = sqrt(power);
  
  for(int i=0; i < VECTOR_SIZE; i++){

    vector[i] = vector[i]/magnitude; 
  }

  return vector; 

}





void print_particle_locations(size_t step ) {
  printf("Particle locations: %d\n", step);
    
    for (size_t i = 0; i < NUM_PARTICLES; i++) {
      
        printf("Particle %zu: x = %.6f, y = %.6f, z = %.6f\n"
	       , i , PARTICLES[i].x, PARTICLES[i].y, PARTICLES[i].z);
	printf("\tMomentum: x = %.6f, y = %.6f, z = %.6f\n"
	       , i , PARTICLES[i].x_momentum, PARTICLES[i].y_momentum, PARTICLES[i].z_momentum);
	printf("\tMomentum Bound: x = %.6f, y = %.6f, z = %.6f\n"
	       , i , PARTICLES[i].x_momentum_bond, PARTICLES[i].y_momentum_bond, PARTICLES[i].z_momentum_bond);
			

    }

    printf("\n");
}



void free_matrix_copy(){



  for(int i=0; i < NUM_PARTICLES; i++){

    free(NEIGHBOR_MATRIX_COPY[i]);

  }

  free(NEIGHBOR_MATRIX_COPY);
  
}

void free_particles(){

  free(PARTICLES);

}

void free_matrix(){

    for(int i=0; i < NUM_PARTICLES; i++){

      free(NEIGHBOR_MATRIX[i]);

    }

  free(NEIGHBOR_MATRIX);
  

}

void free_bond_matrix(){

      for(int i=0; i < NUM_PARTICLES; i++){

      free(BOND_MATRIX[i]);

    }

  free(BOND_MATRIX);
}

void free_bond_matrix_copy(){

    for(int i=0; i < NUM_PARTICLES; i++){

      free(BOND_MATRIX_COPY[i]);

    }

  free(BOND_MATRIX_COPY);
}

void make_matrix_copy(){


        #pragma omp parallel for 
        for (int i = 0; i < NUM_PARTICLES; i++) {

	  for(int j=0; j < NUM_PARTICLES; j++){

	    NEIGHBOR_MATRIX_COPY[i][j]=NEIGHBOR_MATRIX[i][j];

	  }
            
    }


}

void make_bond_matrix_copy(){


        #pragma omp parallel for 
        for (int i = 0; i < NUM_PARTICLES; i++) {

	  for(int j=0; j < NUM_PARTICLES; j++){

	    BOND_MATRIX_COPY[i][j]=BOND_MATRIX[i][j];

	  }
            
    }


}


void reset_neighbor_list(){

  for(int i =0; i < NUM_PARTICLES; i++){

    for(int j =0; j < NUM_PARTICLES; j++){

      NEIGHBOR_MATRIX[i][j] = false;

    }
  }
  

}


void create_naive_neighbor_list(){

 
  for(int i=0; i< NUM_PARTICLES; i++){

   
    for(int j=0; j< NUM_PARTICLES; j++){

      if(i!=j){


	NEIGHBOR_MATRIX[i][j]=true;

      }
    }

  }

}

void create_naive_bond_list(){

 
  for(int i=0; i< NUM_PARTICLES; i++){

   
    for(int j=0; j< NUM_PARTICLES; j++){

      if(i!=j){


	BOND_MATRIX[i][j]=true;

      }
    }

  }

}


void create_distance_based_bond_list(double radius){

  for(int i=0; i< NUM_PARTICLES; i++){

    
    for(int j=0; j< NUM_PARTICLES; j++){

      if(i!=j){

	double distance = get_distance(PARTICLES[i], PARTICLES[j]);
	if(distance<=radius){
	 
	  BOND_MATRIX[i][j]=true;
	  
	}

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
	 
	  NEIGHBOR_MATRIX[i][j]=true;
	  
	}

      }
    }

  }

}




void center_force(int j){

  particle center;

  double ratio = 0.0001; 
  
  center.x = 0.0;
  center.y = 0.0;
  center.z = 0.0;

  
  double *vector = get_vector(center, PARTICLES[j]);

  PARTICLES[j].x_momentum = PARTICLES[j].x_momentum + -(vector[0])* ratio;
  PARTICLES[j].y_momentum = PARTICLES[j].y_momentum + -(vector[1])* ratio;
  PARTICLES[j].z_momentum = PARTICLES[j].z_momentum + -(vector[2])* ratio;
}

void newton_force(int j){

 
  for(int i = j+1; i < NUM_PARTICLES; i++ ){


  //if there is an edge connecting the particles;
    
  if(NEIGHBOR_MATRIX_COPY[j][i]){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    double distance  = get_distance(p_j, p_i);
    double potential = newton_potential(distance);
    double *vector   = get_normalized(get_vector(p_j, p_i));
    
   
    PARTICLES[i].x_momentum = PARTICLES[i].x_momentum + -(vector[0] * potential);
    PARTICLES[i].y_momentum = PARTICLES[i].y_momentum + -(vector[0] * potential);
    PARTICLES[i].z_momentum = PARTICLES[i].z_momentum + -(vector[0] * potential);

    PARTICLES[j].x_momentum = PARTICLES[j].x_momentum + (vector[0] * potential);
    PARTICLES[j].y_momentum = PARTICLES[j].y_momentum + (vector[0] * potential);
    PARTICLES[j].z_momentum = PARTICLES[j].z_momentum + (vector[0] * potential);

  
    

    NEIGHBOR_MATRIX_COPY[i][j]=false;
    NEIGHBOR_MATRIX_COPY[j][i]=false;
    }
  }

}

void bond_force(int j){

  
 
  for(int i =j+1; i < NUM_PARTICLES; i++ ){


      
  if(BOND_MATRIX_COPY[j][i]){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    double distance  = get_distance(p_j, p_i);
    double potential = spring_potential(distance);
    double *vector   = get_normalized(get_vector(p_j, p_i));
    
    potential *= BOND_RATIO_MOMENTUM;
    
    PARTICLES[i].x_momentum_bond = PARTICLES[i].x_momentum_bond * DECAY_RATIO + -(vector[0] * potential);
    PARTICLES[i].y_momentum_bond = PARTICLES[i].y_momentum_bond * DECAY_RATIO + -(vector[1] * potential);
    PARTICLES[i].z_momentum_bond = PARTICLES[i].z_momentum_bond * DECAY_RATIO + -(vector[2] * potential);

    PARTICLES[j].x_momentum_bond = PARTICLES[j].x_momentum_bond * DECAY_RATIO + (vector[0] * potential);
    PARTICLES[j].y_momentum_bond = PARTICLES[j].y_momentum_bond * DECAY_RATIO + (vector[1] * potential);
    PARTICLES[j].z_momentum_bond = PARTICLES[j].z_momentum_bond * DECAY_RATIO + (vector[2] * potential);

    /*printf("x: %f y: %f z: %f\n", PARTICLES[i].x_momentum_bond,
	                          PARTICLES[i].y_momentum_bond,
	                          PARTICLES[i].z_momentum_bond);
    printf("vec_x: %f vec_y: %f vec_z: %f\n", vector[0], vector[1], vector[2] );
    printf("potential: %f \n", potential); 
    */                            
   

    BOND_MATRIX_COPY[i][j]=false;
    BOND_MATRIX_COPY[j][i]=false;
    }
  }

  
}




void lennard_force(int j){

  for(int i = j+1; i < NUM_PARTICLES; i++ ){


  //if there is an edge connecting the particles;
    
  if(NEIGHBOR_MATRIX_COPY[j][i]){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    double distance  = get_distance(p_j, p_i);
    double potential = lennard_jones_potential(distance);
    double *vector   = get_normalized(get_vector(p_j, p_i));

    potential *=GLOBAL_RATIO_MOMENTUM;
 
    PARTICLES[i].x_momentum = PARTICLES[i].x_momentum + -(vector[0] * potential);
    PARTICLES[i].y_momentum = PARTICLES[i].y_momentum + -(vector[1] * potential);
    PARTICLES[i].z_momentum = PARTICLES[i].z_momentum + -(vector[2] * potential);

    PARTICLES[j].x_momentum = PARTICLES[j].x_momentum + (vector[0] * potential);
    PARTICLES[j].y_momentum = PARTICLES[j].y_momentum + (vector[1] * potential);
    PARTICLES[j].z_momentum = PARTICLES[j].z_momentum + (vector[2] * potential);


    

    NEIGHBOR_MATRIX_COPY[i][j]=false;
    NEIGHBOR_MATRIX_COPY[j][i]=false;
    
    }
  }

}


void move(int j){

  PARTICLES[j].x =PARTICLES[j].x + GLOBAL_RATIO_MOVE * PARTICLES[j].x_momentum;
  PARTICLES[j].y =PARTICLES[j].y + GLOBAL_RATIO_MOVE * PARTICLES[j].y_momentum;
  PARTICLES[j].z =PARTICLES[j].z + GLOBAL_RATIO_MOVE * PARTICLES[j].z_momentum;


}

void move_bond(int j){

 

  PARTICLES[j].x =PARTICLES[j].x + GLOBAL_RATIO_MOVE * PARTICLES[j].x_momentum;
  PARTICLES[j].y =PARTICLES[j].y + GLOBAL_RATIO_MOVE * PARTICLES[j].y_momentum;
  PARTICLES[j].z =PARTICLES[j].z + GLOBAL_RATIO_MOVE * PARTICLES[j].z_momentum;
  PARTICLES[j].x =PARTICLES[j].x + BOND_RATIO_MOVE * PARTICLES[j].x_momentum_bond;
  PARTICLES[j].y =PARTICLES[j].y + BOND_RATIO_MOVE * PARTICLES[j].y_momentum_bond;
  PARTICLES[j].z =PARTICLES[j].z + BOND_RATIO_MOVE * PARTICLES[j].z_momentum_bond;





}



void traverse(void (*f)(int)){
  
  #pragma omp parallel for schedule(dynamic) 
  for(int i=0; i < NUM_PARTICLES; i++){

    f(i);

  }

}




void write_particle_positions(const char *filename) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return;
    }

   
    
   
      
      fprintf(file, "%d\n\n", NUM_PARTICLES);
      // #pragma omp parallel for schedule(dynamic) 
      for(int j =0; j <NUM_PARTICLES; j++ ){


	  fprintf(file, "%s %.6f %.6f %.6f\n", PARTICLES[j].color, PARTICLES[j].x, PARTICLES[j].y, PARTICLES[j].z);
	  

      }

    

    fclose(file);

    

    
}


void delete_file(const char *filename){

    if (remove(filename) == 0)
        printf("Deleted successfully\n");
    else
        printf("Unable to delete the file\n");
 

}

void free_resources(){

  free_matrix_copy();
  free_matrix();
  free_particles();
  free_bond_matrix();
  free_bond_matrix_copy();

}
