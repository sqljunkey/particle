#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>


double rand_zero_to_one() {
    return rand() / (RAND_MAX + 1.0);
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


double get_value(int argc, char *argv[], const char * command , char * error_message ) {

  double result = -1.0;
  char *ptr;
  
  for(int i = 0; i < argc; i++){

    if( strstr(argv[i], command)){
      if(i+1<argc){
	result = strtod(argv[i+1], &ptr);
	break;
      }
    } 
  }

  if(result <0){

    printf(error_message);
    exit(1);

  }

  return result;

}


double get_mean_value( int argc, char *argv[]){


  return get_value(argc, argv, "-m", "please insert -m parameter");

}

double get_standard_deviation( int argc, char *argv[]){


  return get_value(argc, argv, "-std", "please insert -std parameter");
  
}

double get_number_of_particle( int argc , char * argv[]){

  return get_value(argc, argv, "-n", "please insert -n parameter");

}

double get_number_of_dimension(int argc, char * argv[]){

  return get_value(argc, argv, "-d", "please insert -d parameter");
    
}

int main(int argc, char * argv[] ){

 srand(time(NULL));

 double mean = get_mean_value(argc, argv);
 double std  = get_standard_deviation(argc, argv);
 double num_of_particles  = get_number_of_particle(argc, argv);
 double num_of_dimensions=  get_number_of_dimension(argc, argv);

 printf("<particle_system>");

 for(int i = 0; i < num_of_particles; i++){


   printf("<particle><id>%d</id>",i+1);

   double * values = malloc(sizeof(double) * num_of_dimensions);
   
   for(int d = 0; d < num_of_dimensions; d++){
     values[d] =  truncated_normal(mean, std);

   }
   
   for(int d = 0; d < num_of_dimensions; d++){

     printf("<dimension>%lf</dimension>", values[d]);
    
   }

   for(int d = 0; d < num_of_dimensions; d++){

     printf("<previous_dimension>%lf</previous_dimension>", values[d]);
    
   }

   for(int d =0; d < num_of_dimensions; d++){

     printf("<velocity>0.0</velocity>");
   }

   
   printf("</particle>");

   free(values);
 }

 printf("</particle_system>");



  return 0; 
}
