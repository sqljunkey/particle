
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


void init_particles(particles *particles, int num_particles){

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
