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
  int group_id;
}particle;



particle *PARTICLES;
int NUM_PARTICLES;
bool *NEIGHBOR_MATRIX;
bool *NEIGHBOR_MATRIX_COPY;
bool *BOND_MATRIX;
bool *BOND_MATRIX_COPY;

double GLOBAL_RATIO_MOVE;
double GLOBAL_RATIO_MOMENTUM;
double BOND_RATIO_MOVE;
double BOND_RATIO_MOMENTUM;

double NEWTON_GAP;

int NUM_FRAMES;
double BOND_DECAY_RATIO;
double GLOBAL_DECAY_RATIO;
double SPREAD;
int VECTOR_SIZE = 3;
int BASE_GROUP = -1;
int COLOR_NUM; 




void init(){

  srand(time(NULL));

 int num_threads = omp_get_max_threads();

 printf("Setting number of threads to maximum: %d\n", num_threads);
    
 omp_set_num_threads(num_threads);

 COLOR_NUM = rand()%3;
 
}


/* free stuff */

void free_matrix_copy() {
    free(NEIGHBOR_MATRIX_COPY);
}

void free_particles() {
    free(PARTICLES);
}

void free_matrix() {
    free(NEIGHBOR_MATRIX);
}

void free_bond_matrix() {
    free(BOND_MATRIX);
}

void free_bond_matrix_copy() {
    free(BOND_MATRIX_COPY);
}

void free_matrices(){

  free_matrix_copy();
  free_matrix();
  free_bond_matrix();
  free_bond_matrix_copy();

}

void free_resources(){

  free_matrix_copy();
  free_matrix();
  free_particles();
 
  free_bond_matrix();
  free_bond_matrix_copy();


}

/* end free stuff */ 




void init_matrices() {
  
    NEIGHBOR_MATRIX = malloc(NUM_PARTICLES * NUM_PARTICLES * sizeof(bool));
    if (NEIGHBOR_MATRIX == NULL) {
        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);
    }

    NEIGHBOR_MATRIX_COPY = malloc(NUM_PARTICLES * NUM_PARTICLES * sizeof(bool));
    if (NEIGHBOR_MATRIX_COPY == NULL) {
        perror("Failed to allocate memory for NEIGHBOR_MATRIX_COPY");
        free(NEIGHBOR_MATRIX); 
        exit(EXIT_FAILURE);
    }

    BOND_MATRIX = malloc(NUM_PARTICLES * NUM_PARTICLES * sizeof(bool));
    if (BOND_MATRIX == NULL) {
        perror("Failed to allocate memory for BOND_MATRIX");
        free(NEIGHBOR_MATRIX);
        free(NEIGHBOR_MATRIX_COPY);
        exit(EXIT_FAILURE);
    }

    BOND_MATRIX_COPY = malloc(NUM_PARTICLES * NUM_PARTICLES * sizeof(bool));
    if (BOND_MATRIX_COPY == NULL) {
        perror("Failed to allocate memory for BOND_MATRIX_COPY");
        free(NEIGHBOR_MATRIX);
        free(NEIGHBOR_MATRIX_COPY);
        free(BOND_MATRIX);
        exit(EXIT_FAILURE);
    }

   
    for (int i = 0; i < NUM_PARTICLES * NUM_PARTICLES; i++) {
        NEIGHBOR_MATRIX[i] = false;
        NEIGHBOR_MATRIX_COPY[i] = false;
        BOND_MATRIX[i] = false;
        BOND_MATRIX_COPY[i] = false;
    }
}


char* get_next_color(){


  COLOR_NUM++;

  if(COLOR_NUM>=3){

    COLOR_NUM=0;

  }
  
    switch (COLOR_NUM) {
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

void set_group_id(int old_group_id, int new_group_id){

  for(int i=0; i < NUM_PARTICLES; i++){

    if(PARTICLES[i].group_id==old_group_id)
         PARTICLES[i].group_id=new_group_id;

    }
}




int contains(int* arr, int size, int value) {
    for (int i = 0; i < size; i++) {
        if (arr[i] == value) {
            return 1;
        }
    }
    return 0; 
}


int* get_group_ids(int* count) {
    int* distinct_ids = malloc(NUM_PARTICLES * sizeof(int)); 
    if (distinct_ids == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        *count = 0;
        return NULL;
    }

    int unique_count = 0;

   
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int group_id = PARTICLES[i].group_id;
        if (!contains(distinct_ids, unique_count, group_id)) {
            distinct_ids[unique_count] = group_id;
            unique_count++;
        }
    }

   
    distinct_ids = realloc(distinct_ids, unique_count * sizeof(int));
    if (distinct_ids == NULL && unique_count > 0) {
        fprintf(stderr, "Memory reallocation failed\n");
        *count = 0;
        return NULL;
    }

    *count = unique_count;
    return distinct_ids;
}


void display_groups(){

  int count = 0;
  int* groups= get_group_ids(&count);

  for(int i =0; i < count; i++){

    printf("Group %d Found.\n", groups[i]);

  }

  free(groups); 

}







void assign_group_colors() {
    int num_groups = 0;

    int *groups = get_group_ids(&num_groups); 

    if(num_groups!=0){
   
    char **group_colors = malloc(num_groups * sizeof(char*));
    
  
    for (int id = 0; id < num_groups; id++) {
        group_colors[id] = get_next_color();
        printf("Group %d was assigned color: %s\n", id, group_colors[id]);
    }



    for(int f=0; f < num_groups; f++){
      
      for (int i = 0; i < NUM_PARTICLES; i++) {


	if(PARTICLES[i].group_id == groups[f] ){
           PARTICLES[i].color = group_colors[f];  

	}
     }
    }
       
     


    free(groups); 
    free(group_colors);
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

double random_double() {
    return (double)rand() / (double)RAND_MAX;
}


void init_disc_particles(double x_center, double y_center, double z_height, double disc_radius) {
    PARTICLES = malloc(NUM_PARTICLES * sizeof(particle));

    if (PARTICLES == NULL) {
        perror("Failed to allocate memory for PARTICLES");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
       
        double r = sqrt(random_double()) * disc_radius; 
        double theta = random_double() * 2.0 * M_PI;    

        
        double x_offset = r * cos(theta);
        double y_offset = r * sin(theta);

        
        PARTICLES[i].x = x_center + x_offset;
        PARTICLES[i].y = y_center + y_offset;
        PARTICLES[i].z = z_height;

        
        PARTICLES[i].x_momentum = 0.0;
        PARTICLES[i].y_momentum = 0.0;
        PARTICLES[i].z_momentum = 0.0;

       
        PARTICLES[i].mass = 1.0;
        PARTICLES[i].color = random_color();
	PARTICLES[i].group_id=BASE_GROUP;
    }

   
    init_matrices(); 
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
    PARTICLES[i].group_id = BASE_GROUP;
    }

  

    init_matrices();

}



void add_particles(int particle_num, double x_init_pos, double y_init_pos, double z_init_pos) {


    free_matrices();

  
    particle* temp_particles = realloc(PARTICLES, (NUM_PARTICLES + particle_num) * sizeof(particle));


    if (temp_particles == NULL) {
        perror("Failed to reallocate memory for PARTICLES");
        exit(EXIT_FAILURE);
    }


    PARTICLES = temp_particles;


    for (int i = NUM_PARTICLES; i < NUM_PARTICLES + particle_num; i++) {
        PARTICLES[i].x = truncated_normal(0.0, SPREAD) + x_init_pos;
        PARTICLES[i].y = truncated_normal(0.0, SPREAD) + y_init_pos;
        PARTICLES[i].z = truncated_normal(0.0, SPREAD) + z_init_pos;

        PARTICLES[i].x_momentum = 0.0;
        PARTICLES[i].y_momentum = 0.0;
        PARTICLES[i].z_momentum = 0.0;

        PARTICLES[i].x_momentum_bond = 0.0;
        PARTICLES[i].y_momentum_bond = 0.0;
        PARTICLES[i].z_momentum_bond = 0.0;

        PARTICLES[i].mass = 1.0;
        PARTICLES[i].color = random_color();
	PARTICLES[i].group_id = BASE_GROUP;
    }


    NUM_PARTICLES += particle_num;


    init_matrices();
}


void init_ring_particles(double x_init_pos, double y_init_pos, double z_init_pos, double inner_radius, double outer_radius) {
    PARTICLES = malloc(NUM_PARTICLES * sizeof(particle));

    if (PARTICLES == NULL) {
        perror("Failed to allocate memory for PARTICLES");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
        
        double r = inner_radius + ((double)rand() / RAND_MAX) * (outer_radius - inner_radius); 
       
        double theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;  

       
        double x_offset = r * cos(theta);
        double y_offset = r * sin(theta);

       
        PARTICLES[i].x = x_init_pos + x_offset;
        PARTICLES[i].y = y_init_pos + y_offset;
        PARTICLES[i].z = z_init_pos;

       
        PARTICLES[i].x_momentum = 0.0;
        PARTICLES[i].y_momentum = 0.0;
        PARTICLES[i].z_momentum = 0.0;

        PARTICLES[i].x_momentum_bond = 0.0;
        PARTICLES[i].y_momentum_bond = 0.0;
        PARTICLES[i].z_momentum_bond = 0.0;

       
        PARTICLES[i].mass = 1.0;
        PARTICLES[i].color = random_color();
        PARTICLES[i].group_id = BASE_GROUP;
    }

    init_matrices();
}





void init_lattice_particles(double x_init_pos, double y_init_pos, double z_init_pos) {

   
    int particles_per_axis = (int) ceil(cbrt(NUM_PARTICLES));
    PARTICLES = malloc((NUM_PARTICLES+1) * sizeof(particle));

    int count = 0;
    for (int i = 0; i < particles_per_axis && count < NUM_PARTICLES; ++i) {
        for (int j = 0; j < particles_per_axis && count < NUM_PARTICLES; ++j) {
            for (int k = 0; k < particles_per_axis && count < NUM_PARTICLES; ++k) {
	      PARTICLES[count].x = (i * SPREAD)+x_init_pos;
	      PARTICLES[count].y = (j * SPREAD)+y_init_pos;
	      PARTICLES[count].z = (k * SPREAD)+z_init_pos;

		PARTICLES[count].x_momentum =0.0;
		PARTICLES[count].y_momentum =0.0;
		PARTICLES[count].z_momentum =0.0;

		PARTICLES[count].x_momentum_bond = 0.0;
		PARTICLES[count].y_momentum_bond = 0.0;
		PARTICLES[count].z_momentum_bond = 0.0;

		PARTICLES[count].mass = 1.0;
		PARTICLES[count].color =random_color();
		PARTICLES[i].group_id = BASE_GROUP;
		 
                count++;
                if (count >= NUM_PARTICLES) break;
            }
        }
    }

    free_matrices(); 
    init_matrices();
   
}

void set_momentum(double x, double y, double z){

  for(int i =0; i < NUM_PARTICLES; i++){
    PARTICLES[i].x_momentum=x;
    PARTICLES[i].y_momentum=y;
    PARTICLES[i].z_momentum=z;

  }
 
}


void scale_by_group(int x, int y, int z, double scale_factor, double group_center, int group){


  for(int i=0; i < NUM_PARTICLES; i++){

    if(PARTICLES[i].group_id==group){


      if(x==1){

	if(PARTICLES[i].x<group_center){

	  PARTICLES[i].x+=scale_factor; 

	}

      }

      if(y==1){

	if(PARTICLES[i].y<group_center){

	  PARTICLES[i].y+=scale_factor; 

	}

      }


      if(y==1){

	if(PARTICLES[i].y<group_center){

	  PARTICLES[i].y+=scale_factor; 

	}

      }




      
    }

  }

  
}

void move_particles_by_group(double x, double y, double z, int group){

  for(int i=0; i< NUM_PARTICLES; i++){
    
    if(PARTICLES[i].group_id == group){

      PARTICLES[i].x += x;
      PARTICLES[i].y += y;
      PARTICLES[i].z += z;

    }

  }


}


void set_momentum_by_group(double x, double y, double z, int group){

  for(int i =0; i < NUM_PARTICLES; i++){

    if(PARTICLES[i].group_id==group){

    PARTICLES[i].x_momentum= x;
    PARTICLES[i].y_momentum=y;
    PARTICLES[i].z_momentum=z;
       

    }
  }
 
}


void give_bond_momentum_by_group(double x, double y, double z, int group){

  for(int i =0; i < NUM_PARTICLES; i++){

    if(PARTICLES[i].group_id==group){

    PARTICLES[i].x_momentum_bond= x;
    PARTICLES[i].y_momentum_bond=y;
    PARTICLES[i].z_momentum_bond=z;
       

    }
  }
 
}



void give_random_momentum( double x, double y, double z){
    for(int i =0; i < NUM_PARTICLES; i++){
      PARTICLES[i].x_momentum=   fabs(truncated_normal(0.0, 1.0)) * x;
      PARTICLES[i].y_momentum=   fabs(truncated_normal(0.0, 1.0)) * y;
      PARTICLES[i].z_momentum=   fabs(truncated_normal(0.0, 1.0)) * z;

  }
  
}

void give_random_mass(double center, double std){

    for(int i =0; i < NUM_PARTICLES; i++){
      PARTICLES[i].mass=fabs(truncated_normal(center, std));

  }

}


void set_random_mass_by_group(double center, double std, int group){

    for(int i =0; i < NUM_PARTICLES; i++){

      if(PARTICLES[i].group_id==group)
         PARTICLES[i].mass=fabs(truncated_normal(center, std));

      
  }

}


void set_mass_by_group(double mass, int group){

    for(int i =0; i < NUM_PARTICLES; i++){
      if(PARTICLES[i].group_id == group)
         PARTICLES[i].mass=mass;

  }

}


void init_cylindrical_particles(double x_init_pos, double y_init_pos, double z_init_pos, double radius, double height) {
    PARTICLES = malloc(NUM_PARTICLES * sizeof(particle));

    if (PARTICLES == NULL) {
        perror("Failed to allocate memory for PARTICLES");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
      
        double r = sqrt((double)rand() / RAND_MAX) * radius;  
       
        double theta = ((double)rand() / RAND_MAX) * 2.0 * M_PI;  
       
        double z = ((double)rand() / RAND_MAX) * height;  

       
        double x_offset = r * cos(theta);
        double y_offset = r * sin(theta);

       
        PARTICLES[i].x = x_init_pos + x_offset;
        PARTICLES[i].y = y_init_pos + y_offset;
        PARTICLES[i].z = z_init_pos + z;

      
        PARTICLES[i].x_momentum = 0.0;
        PARTICLES[i].y_momentum = 0.0;
        PARTICLES[i].z_momentum = 0.0;

        PARTICLES[i].x_momentum_bond = 0.0;
        PARTICLES[i].y_momentum_bond = 0.0;
        PARTICLES[i].z_momentum_bond = 0.0;

        
        PARTICLES[i].mass = 1.0;
        PARTICLES[i].color = random_color();
        PARTICLES[i].group_id = BASE_GROUP;
    }

    init_matrices();
}




void init_spherical_particles(double x_init_pos, double y_init_pos, double z_init_pos, double radius) {
    PARTICLES = malloc(NUM_PARTICLES * sizeof(particle));

    if (PARTICLES == NULL) {
        perror("Failed to allocate memory for PARTICLES");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
       
        double r = cbrt((double)rand() / RAND_MAX) * radius;  
        double theta = ((double)rand() / RAND_MAX) * M_PI;   
        double phi = ((double)rand() / RAND_MAX) * 2.0 * M_PI; 

      
        double x_offset = r * sin(theta) * cos(phi);
        double y_offset = r * sin(theta) * sin(phi);
        double z_offset = r * cos(theta);

       
        PARTICLES[i].x = x_init_pos + x_offset;
        PARTICLES[i].y = y_init_pos + y_offset;
        PARTICLES[i].z = z_init_pos + z_offset;

       
        PARTICLES[i].x_momentum = 0.0;
        PARTICLES[i].y_momentum = 0.0;
        PARTICLES[i].z_momentum = 0.0;

        PARTICLES[i].x_momentum_bond = 0.0;
        PARTICLES[i].y_momentum_bond = 0.0;
        PARTICLES[i].z_momentum_bond = 0.0;

        
        PARTICLES[i].mass = 1.0;
        PARTICLES[i].color = random_color();
	PARTICLES[i].group_id = BASE_GROUP;
    }

    init_matrices();
}




void init_particles(double x_init_pos, double y_init_pos, double z_init_pos ){



      


    PARTICLES = malloc(NUM_PARTICLES * sizeof(particle));

    if(PARTICLES == NULL){
      
        perror("Failed to allocate memory for NEIGHBOR_MATRIX");
        exit(EXIT_FAILURE);

    }
    
    for(int i=0; i < NUM_PARTICLES; i++){

   
    PARTICLES[i].x = truncated_normal(0.0, SPREAD) + x_init_pos;
    PARTICLES[i].y = truncated_normal(0.0, SPREAD) + y_init_pos;
    PARTICLES[i].z = truncated_normal(0.0, SPREAD) + z_init_pos;

  
    
    PARTICLES[i].x_momentum =0.0;
    PARTICLES[i].y_momentum =0.0;
    PARTICLES[i].z_momentum =0.0;

    PARTICLES[i].x_momentum_bond = 0.0;
    PARTICLES[i].y_momentum_bond = 0.0;
    PARTICLES[i].z_momentum_bond = 0.0;

    PARTICLES[i].mass = 1.0;
    PARTICLES[i].color =random_color();
    PARTICLES[i].group_id = BASE_GROUP;


    
  }



    init_matrices();



}


double push_potential(double distance){
 double result = 0.0;

  result = 1.0 / pow(distance,2.0);

  return -result;


}

double spring_potential(double distance){


  double result = 0.0;
  double equilibrium_distance = 400.0;
  double max_distance = 100000.0;


    
    if(distance < 1000000){
        result =distance - equilibrium_distance;
    }


		
    return result;
		


}

double newton_potential(double distance) {
  double newton = 0.0;


  if(distance > NEWTON_GAP){

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





void make_matrix_copy() {
    #pragma omp parallel for
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_PARTICLES; j++) {
            NEIGHBOR_MATRIX_COPY[i * NUM_PARTICLES + j] = NEIGHBOR_MATRIX[i * NUM_PARTICLES + j];
        }
    }
}

void make_bond_matrix_copy() {
    #pragma omp parallel for
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_PARTICLES; j++) {
            BOND_MATRIX_COPY[i * NUM_PARTICLES + j] = BOND_MATRIX[i * NUM_PARTICLES + j];
        }
    }
}

void reset_neighbor_list() {
    for (int i = 0; i < NUM_PARTICLES * NUM_PARTICLES; i++) {
        NEIGHBOR_MATRIX[i] = false;
    }
}
void create_linked_neighbor_list() {
    reset_neighbor_list();

    for (int i = 0; i < NUM_PARTICLES; i++) {
        if (i > 0) {
	  // NEIGHBOR_MATRIX[i * NUM_PARTICLES + (i - 1)] = true; 
            NEIGHBOR_MATRIX[(i - 1) * NUM_PARTICLES + i] = true; 
        }
        if (i < NUM_PARTICLES - 1) {
	  //NEIGHBOR_MATRIX[i * NUM_PARTICLES + (i + 1)] = true; 
            NEIGHBOR_MATRIX[(i + 1) * NUM_PARTICLES + i] = true; 
        }
    }
}


void create_naive_neighbor_list() {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_PARTICLES; j++) {
            if (i != j) {
                NEIGHBOR_MATRIX[i * NUM_PARTICLES + j] = true;
            }
        }
    }
}

void create_naive_bond_list() {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_PARTICLES; j++) {
            if (i != j) {
                BOND_MATRIX[i * NUM_PARTICLES + j] = true;
            }
        }
    }
}


void create_neighbor_selection(double x, double y, double z, double radius) {
  
    for (int i = 0; i < NUM_PARTICLES; i++) {
       
        double distance_from_center = sqrt(pow(PARTICLES[i].x - x, 2) +
                                           pow(PARTICLES[i].y - y, 2) +
                                           pow(PARTICLES[i].z - z, 2));

        
        if (distance_from_center <= radius) {

	 
          
            for (int j = 0; j < NUM_PARTICLES; j++) {
                if (i != j) {  
                    double distance_from_center_j = sqrt(pow(PARTICLES[j].x - x, 2) +
                                                         pow(PARTICLES[j].y - y, 2) +
                                                         pow(PARTICLES[j].z - z, 2));

                    
                    if (distance_from_center_j <= radius) {
                        BOND_MATRIX[i * NUM_PARTICLES + j] = true;
                        BOND_MATRIX[j * NUM_PARTICLES + i] = true; 
                    }
                }
            }
        }
    }
}


void create_distance_based_bond_list(double radius) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_PARTICLES; j++) {
            if (i != j) {
                double distance = get_distance(PARTICLES[i], PARTICLES[j]);
                if (distance <= radius) {
		  
                    BOND_MATRIX[i * NUM_PARTICLES + j] = true;
                }
            }
        }
    }
}

void create_lattice_neighbor_list(){

      for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_PARTICLES; j++) {
            if (i != j) {
	    }
	}
      }
      



}


void create_distance_based_neighbor_list(double radius) {
    for (int i = 0; i < NUM_PARTICLES; i++) {
        for (int j = 0; j < NUM_PARTICLES; j++) {
            if (i != j) {
	     
                double distance = get_distance(PARTICLES[i], PARTICLES[j]);

	
                if (distance <= radius) {
                    NEIGHBOR_MATRIX[i * NUM_PARTICLES + j] = true;
		   
                }
            }
        }
    }
}



void center_force(int j){

  particle center;


  
  center.x = 0.0;
  center.y = 0.0;
  center.z = 0.0;

  
  double *vector = get_normalized(get_vector(center, PARTICLES[j]));

  PARTICLES[j].x_momentum = PARTICLES[j].x_momentum + -(vector[0]* GLOBAL_RATIO_MOMENTUM)*GLOBAL_DECAY_RATIO * PARTICLES[j].mass;
  PARTICLES[j].y_momentum = PARTICLES[j].y_momentum + -(vector[1]* GLOBAL_RATIO_MOMENTUM)*GLOBAL_DECAY_RATIO * PARTICLES[j].mass;
  PARTICLES[j].z_momentum = PARTICLES[j].z_momentum + -(vector[2]* GLOBAL_RATIO_MOMENTUM)*GLOBAL_DECAY_RATIO * PARTICLES[j].mass;

  free(vector);
}

void newton_force(int j){

 
  for(int i = j+1; i < NUM_PARTICLES; i++ ){



    
    if(NEIGHBOR_MATRIX_COPY[j * NUM_PARTICLES + i]){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    double distance  = get_distance(p_j, p_i);
    double potential = newton_potential(distance);
    double *vector   = get_normalized(get_vector(p_j, p_i));

    potential *=GLOBAL_RATIO_MOMENTUM;
   
    PARTICLES[i].x_momentum = PARTICLES[i].x_momentum*GLOBAL_DECAY_RATIO + -(vector[0] * potential)*PARTICLES[j].mass;
    PARTICLES[i].y_momentum = PARTICLES[i].y_momentum*GLOBAL_DECAY_RATIO + -(vector[1] * potential)*PARTICLES[j].mass;
    PARTICLES[i].z_momentum = PARTICLES[i].z_momentum*GLOBAL_DECAY_RATIO + -(vector[2] * potential)*PARTICLES[j].mass;

    PARTICLES[j].x_momentum = PARTICLES[j].x_momentum*GLOBAL_DECAY_RATIO + (vector[0] * potential)*PARTICLES[i].mass;
    PARTICLES[j].y_momentum = PARTICLES[j].y_momentum*GLOBAL_DECAY_RATIO + (vector[1] * potential)*PARTICLES[i].mass;
    PARTICLES[j].z_momentum = PARTICLES[j].z_momentum*GLOBAL_DECAY_RATIO + (vector[2] * potential)*PARTICLES[i].mass;

  
    

    NEIGHBOR_MATRIX_COPY[i + NUM_PARTICLES + j]=false;
    NEIGHBOR_MATRIX_COPY[j + NUM_PARTICLES + i]=false;
    free(vector);
    }
  }

}


void rotate_particles_by_group(double axis_x, double axis_y, double axis_z, double theta, int group) {
    
    double axis[3] = {axis_x, axis_y, axis_z};
    double axis_magnitude = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    
    if (axis_magnitude == 0) {
        return;
    }
    
    
    axis[0] /= axis_magnitude;
    axis[1] /= axis_magnitude;
    axis[2] /= axis_magnitude;

    double u_x = axis[0];
    double u_y = axis[1];
    double u_z = axis[2];

    
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);
    double one_minus_cos_theta = 1.0 - cos_theta;

    
    for (int j = 0; j < NUM_PARTICLES; j++) {

      if(PARTICLES[j].group_id==group){


	
        particle *p = &PARTICLES[j];

        
        double x = p->x;
        double y = p->y;
        double z = p->z;

        
        double new_x = (cos_theta + u_x * u_x * one_minus_cos_theta) * x +
                       (u_x * u_y * one_minus_cos_theta - u_z * sin_theta) * y +
                       (u_x * u_z * one_minus_cos_theta + u_y * sin_theta) * z;

        double new_y = (u_y * u_x * one_minus_cos_theta + u_z * sin_theta) * x +
                       (cos_theta + u_y * u_y * one_minus_cos_theta) * y +
                       (u_y * u_z * one_minus_cos_theta - u_x * sin_theta) * z;

        double new_z = (u_z * u_x * one_minus_cos_theta - u_y * sin_theta) * x +
                       (u_z * u_y * one_minus_cos_theta + u_x * sin_theta) * y +
                       (cos_theta + u_z * u_z * one_minus_cos_theta) * z;

       
        p->x = new_x;
        p->y = new_y;
        p->z = new_z;

        
        double p_x = p->x_momentum;
        double p_y = p->y_momentum;
        double p_z = p->z_momentum;

        double new_p_x = (cos_theta + u_x * u_x * one_minus_cos_theta) * p_x +
                         (u_x * u_y * one_minus_cos_theta - u_z * sin_theta) * p_y +
                         (u_x * u_z * one_minus_cos_theta + u_y * sin_theta) * p_z;

        double new_p_y = (u_y * u_x * one_minus_cos_theta + u_z * sin_theta) * p_x +
                         (cos_theta + u_y * u_y * one_minus_cos_theta) * p_y +
                         (u_y * u_z * one_minus_cos_theta - u_x * sin_theta) * p_z;

        double new_p_z = (u_z * u_x * one_minus_cos_theta - u_y * sin_theta) * p_x +
                         (u_z * u_y * one_minus_cos_theta + u_x * sin_theta) * p_y +
                         (cos_theta + u_z * u_z * one_minus_cos_theta) * p_z;

       
        p->x_momentum = 0.0;
        p->y_momentum = 0.0;
        p->z_momentum = 0.0;

	p->x_momentum_bond = 0.0;
        p->y_momentum_bond = 0.0;
        p->z_momentum_bond = 0.0;
      }
    }
}


void relax_force(int j){

  for(int i = j+1; i < NUM_PARTICLES; i++ ){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    if(BOND_MATRIX_COPY[j * NUM_PARTICLES + i]){
       double distance  = get_distance(p_j, p_i);

       if(distance < 19.0){

	 
	 PARTICLES[i].x_momentum_bond = 0.0;
         PARTICLES[i].y_momentum_bond = 0.0;
         PARTICLES[i].z_momentum_bond = 0.0;

         PARTICLES[j].x_momentum_bond = 0.0;
         PARTICLES[j].y_momentum_bond = 0.0;
         PARTICLES[j].z_momentum_bond = 0.0;

       }

          BOND_MATRIX_COPY[i + NUM_PARTICLES + j]=false;
          BOND_MATRIX_COPY[j + NUM_PARTICLES + i]=false;
       
    }

  }
  

}
void push_force(int j){

 
  for(int i = j+1; i < NUM_PARTICLES; i++ ){



    
    if(BOND_MATRIX_COPY[j * NUM_PARTICLES + i]){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    double distance  = get_distance(p_j, p_i);
    double potential = push_potential(distance);
    double *vector   = get_normalized(get_vector(p_j, p_i));
    
    potential *= BOND_RATIO_MOMENTUM;
    
    PARTICLES[i].x_momentum_bond = PARTICLES[i].x_momentum_bond*BOND_DECAY_RATIO + -(vector[0] * potential)*PARTICLES[j].mass;
    PARTICLES[i].y_momentum_bond = PARTICLES[i].y_momentum_bond*BOND_DECAY_RATIO + -(vector[1] * potential)*PARTICLES[j].mass;
    PARTICLES[i].z_momentum_bond = PARTICLES[i].z_momentum_bond*BOND_DECAY_RATIO + -(vector[2] * potential)*PARTICLES[j].mass;

    PARTICLES[j].x_momentum_bond = PARTICLES[j].x_momentum_bond*BOND_DECAY_RATIO + (vector[0] * potential)*PARTICLES[i].mass;
    PARTICLES[j].y_momentum_bond = PARTICLES[j].y_momentum_bond*BOND_DECAY_RATIO + (vector[1] * potential)*PARTICLES[i].mass;
    PARTICLES[j].z_momentum_bond = PARTICLES[j].z_momentum_bond*BOND_DECAY_RATIO + (vector[2] * potential)*PARTICLES[i].mass;

  
    

    BOND_MATRIX_COPY[i + NUM_PARTICLES + j]=false;
    BOND_MATRIX_COPY[j + NUM_PARTICLES + i]=false;
    free(vector);
    }
  }

}

void bond_force(int j){

  
 
  for(int i =j+1; i < NUM_PARTICLES; i++ ){


      
  if(BOND_MATRIX_COPY[ j * NUM_PARTICLES + i]){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    double distance  = get_distance(p_j, p_i);
    double potential = spring_potential(distance);
    double *vector   = get_normalized(get_vector(p_j, p_i));

    
    
    potential *= BOND_RATIO_MOMENTUM;
    
    PARTICLES[i].x_momentum_bond = PARTICLES[i].x_momentum_bond * BOND_DECAY_RATIO + -(vector[0] * potential);
    PARTICLES[i].y_momentum_bond = PARTICLES[i].y_momentum_bond * BOND_DECAY_RATIO + -(vector[1] * potential);
    PARTICLES[i].z_momentum_bond = PARTICLES[i].z_momentum_bond * BOND_DECAY_RATIO + -(vector[2] * potential);

    PARTICLES[j].x_momentum_bond = PARTICLES[j].x_momentum_bond * BOND_DECAY_RATIO + (vector[0] * potential);
    PARTICLES[j].y_momentum_bond = PARTICLES[j].y_momentum_bond * BOND_DECAY_RATIO + (vector[1] * potential);
    PARTICLES[j].z_momentum_bond = PARTICLES[j].z_momentum_bond * BOND_DECAY_RATIO + (vector[2] * potential);


   

      BOND_MATRIX_COPY[i * NUM_PARTICLES + j]=false;
	BOND_MATRIX_COPY[j * NUM_PARTICLES + i]=false;
    free(vector);
    }
  }

  
}




void lennard_force(int j){

  for(int i = j+1; i < NUM_PARTICLES; i++ ){


  //if there is an edge connecting the particles;
    
  if(NEIGHBOR_MATRIX_COPY[j * NUM_PARTICLES + i]){

    particle p_j = PARTICLES[j];
    particle p_i = PARTICLES[i];

    double distance  = get_distance(p_j, p_i);
    double potential = lennard_jones_potential(distance);
    double *vector   = get_normalized(get_vector(p_j, p_i));

    potential *=GLOBAL_RATIO_MOMENTUM;
 
    PARTICLES[i].x_momentum = PARTICLES[i].x_momentum + -(vector[0] * potential)*PARTICLES[j].mass;
    PARTICLES[i].y_momentum = PARTICLES[i].y_momentum + -(vector[1] * potential)*PARTICLES[j].mass;
    PARTICLES[i].z_momentum = PARTICLES[i].z_momentum + -(vector[2] * potential)*PARTICLES[j].mass;

    PARTICLES[j].x_momentum = PARTICLES[j].x_momentum + (vector[0] * potential)*PARTICLES[i].mass;
    PARTICLES[j].y_momentum = PARTICLES[j].y_momentum + (vector[1] * potential)*PARTICLES[i].mass;
    PARTICLES[j].z_momentum = PARTICLES[j].z_momentum + (vector[2] * potential)*PARTICLES[i].mass;


    NEIGHBOR_MATRIX_COPY[i * NUM_PARTICLES + j]=false;
    NEIGHBOR_MATRIX_COPY[j * NUM_PARTICLES + i]=false;
    free(vector); 
    
    }
  }

}


void move(int j){

  PARTICLES[j].x =PARTICLES[j].x + GLOBAL_RATIO_MOVE * PARTICLES[j].x_momentum;
  PARTICLES[j].y =PARTICLES[j].y + GLOBAL_RATIO_MOVE * PARTICLES[j].y_momentum;
  PARTICLES[j].z =PARTICLES[j].z + GLOBAL_RATIO_MOVE * PARTICLES[j].z_momentum;


}

void bond(int j){
  PARTICLES[j].x =PARTICLES[j].x + BOND_RATIO_MOVE * PARTICLES[j].x_momentum_bond;
  PARTICLES[j].y =PARTICLES[j].y + BOND_RATIO_MOVE * PARTICLES[j].y_momentum_bond;
  PARTICLES[j].z =PARTICLES[j].z + BOND_RATIO_MOVE * PARTICLES[j].z_momentum_bond;
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




void insert_particles(particle *new_particle) {
   
    NUM_PARTICLES++;

  
    PARTICLES = realloc(PARTICLES, NUM_PARTICLES * sizeof(particle));
    if (PARTICLES == NULL) {

     
        perror("Failed to reallocate memory for PARTICLES");
        exit(EXIT_FAILURE);
    }


   
    PARTICLES[NUM_PARTICLES - 1] = *new_particle;

   
    PARTICLES[NUM_PARTICLES - 1].color = random_color();

   
    free_matrices();
  
    init_matrices();
 
}

void pop_particle(int index) {
   
    if (index < 0 || index >= NUM_PARTICLES) {
        fprintf(stderr, "Index out of bounds\n");
        return;
    }

   
    for (int i = index; i < NUM_PARTICLES - 1; i++) {
        PARTICLES[i] = PARTICLES[i + 1];
    }

  
    NUM_PARTICLES--;

  
    PARTICLES = realloc(PARTICLES, NUM_PARTICLES * sizeof(particle));
    if (NUM_PARTICLES > 0 && PARTICLES == NULL) {
        perror("Failed to reallocate memory for PARTICLES");
        exit(EXIT_FAILURE);
    }
}

