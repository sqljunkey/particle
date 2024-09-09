#include "particle_system.h"


int main(int argc, char *argv[] ){



  SPREAD = 22.0;
  NUM_PARTICLES =200;
  GLOBAL_RATIO_MOVE     = 1e-2;
  GLOBAL_RATIO_MOMENTUM = 1e-3;
  BOND_RATIO_MOVE       = 1e-1;
  BOND_RATIO_MOMENTUM   = 1e-3;
  DECAY_RATIO=0.9999f;
  NUM_FRAMES = 10000;
  const char *filename = "file.xyz";
  int counter = 0;

  

  
   init_particles();

  //init_debug_particles(); 
  delete_file(filename); 



  double start_time = omp_get_wtime();

  create_naive_bond_list();
 
  for(int i =0; i < NUM_FRAMES; i++ ){


    //print_particle_locations(i);

    if(i % 2) {

      create_naive_neighbor_list();
          

    }

    /* delete old edge-matrix copy, create new edge matrix copy */ 
      
   
    make_bond_matrix_copy();
    make_matrix_copy();
    
    /* apply force to all particles in system. */

    counter++;
    traverse(newton_force);
    traverse(bond_force);
    traverse(move_bond);


    /* change color */
    if(i==500){

      printf("changing colors %d\n", i);
      change_colors();
    }
    
   
    // system("PAUSE");
    if(counter > 10){
    write_particle_positions(filename);

    counter =0;
    }
    
  }

    double end_time = omp_get_wtime();
    printf("Execution time: %f seconds\n", end_time - start_time);
  
  free_resources();

  printf("Success\n");

  system("java -jar SimpleViewer.jar file.xyz");
  return 0;

}
