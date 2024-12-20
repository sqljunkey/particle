#include<iostream>
#include<unordered_map>
#include<unordered_set>
#include<string>
#include<random>
#include<cmath>
#include<vector>


using namespace std;



double grid_size = 5.0;
double grid_multiple = 10.0;

void set_grid_multiple(int value){

  grid_multiple = value;

}

void set_gridsize(int value){
  grid_size = value; 
}


double neg_zero(double val) {
    return val == 0.0 ? 0.0 : val;
}


struct Key {
    double x, y, z;

    bool operator==(const Key &other) const {
        const double epsilon = 1e-9;
        return std::fabs(x - other.x) < epsilon &&
               std::fabs(y - other.y) < epsilon &&
               std::fabs(z - other.z) < epsilon;
    }
};

struct KeyHash {
    std::size_t operator()(const Key &k) const {
        std::size_t hx = std::hash<double>{}(k.x);
        std::size_t hy = std::hash<double>{}(k.y);
        std::size_t hz = std::hash<double>{}(k.z);

        // Use a better combination function to reduce collisions
        std::size_t seed = 0;
        seed ^= hx + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hy + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hz + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};


/*
struct Key {
    double x, y, z;

    bool operator==(const Key &other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct KeyHash {
    std::size_t operator()(const Key &k) const {
        std::size_t hx = std::hash<double>{}(k.x);
        std::size_t hy = std::hash<double>{}(k.y);
        std::size_t hz = std::hash<double>{}(k.z);
        return hx ^ (hy << 1) ^ (hz << 2); 
    }
};
*/

unordered_map<Key, vector<particle>, KeyHash> m;


void clear_grid(){
  m.clear();
}


Key generate_key(double x, double y, double z) {
    double round_x = neg_zero(round(x * grid_size) / grid_size);
    double round_y = neg_zero(round(y * grid_size) / grid_size);
    double round_z = neg_zero(round(z * grid_size) / grid_size);
    return Key{round_x, round_y, round_z};
}




void insert_location(particle p) {

  
       Key key = generate_key(p.x, p.y, p.z);
       m[key].push_back(p);

 
 }


vector<particle> retrieve_neighbors(particle &p){
  double i,k,j;

  vector<particle> neighbors;
  
  double bound = 1/grid_size* grid_multiple; 

  for(i=p.x-bound; i < p.x+bound; i+=1/grid_size){

    for(k=p.y-bound; k < p.y+bound ; k+=1/grid_size){

      for(j=p.z-bound; j < p.z+bound; j+=1/grid_size){

	Key key = generate_key(i,k,j);

	
	
	if(m.find(key)!= m.end()){
	 

	  neighbors.insert(neighbors.end(), m[key].begin(), m[key].end());

	}

      }
    }


  }
  return neighbors; 
}



