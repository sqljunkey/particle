#include<math.h>

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
