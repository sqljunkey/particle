#include "particle_system.h"


int main(int argc, char *argv[] ){



  NUM_PARTICLES =10;
  RATIO = 0.0001f;
  NUM_FRAMES = 2;
  init_particles();



  /* loop through frames */

  
  for(int i =0; i < NUM_FRAMES; i ){

    /* all even frames create new neighbor list */

    if(i % 2) {

      create_naive_neighbor_list();
          

    }

    /* delete old edge-matrix copy, create new edge matrix copy */ 
      
    free_matrix_copy();
    create_matrix_copy(); 
    
    /* apply force to all particles in system. */
    
    traverse(newton_force);
    write_particles("file.xyz");

    
  }
  

  
  return 0;

}
