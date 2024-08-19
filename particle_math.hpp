#include<cmath>


 
vector<double> normalize(vector<double> v){

  vector<double> normal;

  double power = 0.0;
  for (const auto &x : v) {

         power += pow(x, 2);
   }

   double magnitude = sqrt(power);

   for (const auto &x : v) {

     normal.push_back(x / magnitude);

   }

  return normal;

}

double newton_potential(double distance, double cutoff) {
  double newton = 0.0;

    if (distance > cutoff) {

	newton = 1.0 / pow(distance,2.0);

	}


  return newton;
}

vector<double> create_vector(){

  vector<double> result;
  for(int i =0; i<3; i++){
    result.push_back(0.0);
    
  }
  return result;
}

void swap_vector(particle &p, vector<double> v){

  p.x = v.at(0);
  p.y = v.at(1);
  p.z = v.at(2);

}

vector<double> get_position_as_vector(const particle &p){
  vector<double> result;

  result.push_back(p.x);
  result.push_back(p.y);
  result.push_back(p.z);

  return result; 

}

vector<double> add(const vector<double> &v, const vector<double> &k){

  vector<double> result;

  if(v.size() == 3 && k.size()==3){
    result.push_back(v.at(0) + k.at(0));
    result.push_back(v.at(1) + k.at(1));
    result.push_back(v.at(2) + k.at(2));
  }
  return result;

  
}

vector<double> scale(const vector<double> &v, double scalar){

  vector<double> result;
  
  for(const auto &a: v){

    result.push_back(a * scalar);
  }

  return result;
}



vector<double> inverted_origin_vector(particle &p){

  vector<double> v;

  v.push_back(-p.x);
  v.push_back(-p.y);
  v.push_back(-p.z);
  return v;
}
double distance_from_center(particle &p){

  return sqrt(p.x*p.x
	      +p.y*p.y
	      +p.z*p.z);



}


double lennard_jones_potential(double distance, double cutoff) {

		double lj_pot_1 = 0.0;
		double lj_pot_2 = 0.0;
		double diff = 1E1;

		double epsilon = 1.1;
		double sigma = 10.0;

		lj_pot_1 = epsilon * (pow((sigma / (distance - diff)), 2.0) - 2.0 * pow((sigma / (distance - diff)), 1.0));
		lj_pot_2 = epsilon * (pow((sigma / (distance)), 2.0) - 2.0 * pow((sigma / (distance)), 1.0));

		if(distance >cutoff) {
		  
		        return (lj_pot_2 - lj_pot_1);

		} else {
			return 0.0;
		}
}

vector<double> get_vector(particle p, particle n) {
    
   
  vector<double> v;
   
    double dX = -(p.x) + (n.x);
    double dY = -(p.y) + (n.y);
    double dZ = -(p.z) + (n.z);

    
    v.push_back(dX);
    v.push_back(dY);
    v.push_back(dZ);

    return v;
}


  
double get_distance(const particle p, const particle n){
  double result=0.0;


  double distance_x = -p.x + n.x;
  double distance_y = -p.y + n.y;
  double distance_z = -p.z + n.z;

  
  
     double distance = sqrt( (
			        pow(distance_x, 2.0)
			      + pow(distance_y, 2.0)
			      + pow(distance_z, 2.0)));
     return distance;
}


