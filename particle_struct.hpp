#include<vector>

struct particle{
  long id;
  double x;
  double y;
  double z;

  double x_momentum;
  double y_momentum;
  double z_momentum;

  std::vector<particle> neighbors;

};
