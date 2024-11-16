#include "particle_system.h"


int main(int argc, char *argv[] ){




  init();

  
  SPREAD = 200.0;
  NUM_PARTICLES =400;

  GLOBAL_RATIO_MOVE     = 1e1;
  GLOBAL_RATIO_MOMENTUM = 1e1;

  BOND_RATIO_MOVE = 1e-1;
  BOND_RATIO_MOMENTUM = 1e-2;

  BOND_DECAY_RATIO= .99999;
  GLOBAL_DECAY_RATIO=1.0;
  
  NEWTON_GAP = 50.0;
  NUM_FRAMES = 100000;

  const char *filename = "file.xyz";
  int counter = 0;

  

  
  

  // init_disc_particles(0.0,0.0,0.0, 10000.0); 

  init_spherical_particles(0.0, 0.0, 0.0, 10000.0);
    
  set_group_id(-1, 0);

  //scale_by_group(1,0,0,20.0,0.0,0); 

  set_random_mass_by_group(100000.0,5000.0,0); 

  double theta = 90.0 * (M_PI / 180.0);
  
  rotate_particles_by_group(0.0,1.0,0.0,theta,0);


  SPREAD=10000.0;

 
  add_particles(1000, 50000.0,150000.0,0.0);
  

  set_group_id(-1,1);


  set_random_mass_by_group(2050.0, 100.0,1);

  set_momentum_by_group(-12.5, 2.0,0.0, 1);



  add_particles(1000, -50000.0,200000.0,0.0);
  

  set_group_id(-1,2);


  set_random_mass_by_group(2050.0, 100.0,2);

  set_momentum_by_group(-12.5, 2.0,0.0, 2);





  add_particles(1000, -50000.0,-200000.0,0.0);
  

  set_group_id(-1,3);


  set_random_mass_by_group(2050.0, 100.0,3);

  set_momentum_by_group(12.5, 2.0,0.0, 3);





  assign_group_colors(); 

  
  

  create_naive_neighbor_list();

  display_groups();
  

  delete_file(filename); 



  double start_time = omp_get_wtime();



 
  for(int i =0; i < NUM_FRAMES; i++ ){

     

      counter++;
     
     
      make_matrix_copy();
   
   
      traverse(newton_force);


      double theta = 90.0 * (M_PI / 180.0)*.02;
  

      rotate_particles_by_group(.0,.0,1.0,theta,0);
      

      set_momentum_by_group(0.0,0.0,0.0,0);
      
         
      traverse(move_bond);



      

    if(counter > 10){

      printf("Frame: %d\n", i);
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
