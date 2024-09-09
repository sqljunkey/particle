#include "particle_system.h"


int main(int argc, char *argv[] ){



 
  NUM_PARTICLES =1000;
  RATIO = .01f;
  NUM_FRAMES = 1000;
  const char *filename = "file.xyz";

  /* init particles*/
  
  init_particles();
  delete_file(filename); 


  /* loop through frames */

  double start_time = omp_get_wtime();
 
  for(int i =0; i < NUM_FRAMES; i++ ){

    // printf("Frame: %d\n",i);
    /* all even frames create new neighbor list */

    if(i % 2) {

      create_naive_neighbor_list();
          

    }

    /* delete old edge-matrix copy, create new edge matrix copy */ 
      
   
    make_matrix_copy(); 
    
    /* apply force to all particles in system. */
   
    traverse(lennard_force);
    traverse(move); 
    
    write_particle_positions(filename);

    
  }

    double end_time = omp_get_wtime();
    printf("Execution time: %f seconds\n", end_time - start_time);
  
  free_resources();

  printf("Success\n");
  
  return 0;

}
