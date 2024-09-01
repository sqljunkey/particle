#include "particle_system.h"


int main(int argc, char *argv[] ){



  NUM_PARTICLES =100;
  RATIO = 1.0f;
  NUM_FRAMES = 1000;
  const char *filename = "file.xyz";

  /* init particles*/
  
  init_particles();
  delete_file(filename); 


  /* loop through frames */

 
  for(int i =0; i < NUM_FRAMES; i++ ){

    printf("Frame: %d\n",i);
    /* all even frames create new neighbor list */

    if(i % 2) {

      create_naive_neighbor_list();
          

    }

    /* delete old edge-matrix copy, create new edge matrix copy */ 
      
   
    make_matrix_copy(); 
    
    /* apply force to all particles in system. */
   
    traverse(newton_force);
    traverse(move); 
    
    write_particle_positions(filename);

    
  }
  
  free_matrix_copy();
  free_matrix();
  free_particles();

  printf("Success\n");
  
  return 0;

}
