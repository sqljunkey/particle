#include<iostream>
#include<stdio.h>
#include<unordered_map>
#include<string>
#include<cmath>
#include<vector>
#include <iomanip>

using namespace std;

struct particle{

  int id;
  double x;
  double y;
  double z;

  double prev_x;
  double prev_y;
  double prev_z;

  double vel_x;
  double vel_y;
  double vel_z;

};


struct middle_and_left_over{

  string middle;
  string left_over; 
 
};

middle_and_left_over get_middle(string source, string begin, string end){

  middle_and_left_over result;
  result.middle="";
  result.left_over = "";
  

  size_t end_end_pos = source.find(end) + end.length();
  size_t begin_pos = source.find(begin) + begin.length();
  size_t end_pos = source.find(end)-end.length()+1;


  if(begin_pos!=string::npos && end_pos!=string::npos){

   
    result.middle = source.substr(begin_pos, end_pos);
    result.left_over = source.substr(end_end_pos);
  
  }

  return result; 

}

vector <particle> particles; 

double grid_size = 5.0;
double grid_multiple = 10.0;
double radius = 100;

void set_grid_multiple(int value){

  grid_multiple = value;

}

void set_gridsize(int value){
  grid_size = value; 
}


double neg_zero(double val) {
    return val == 0.0 ? 0.0 : val;
}

double get_distance(particle p, particle l){

 double dx =  p.x - l.x;
 double dy =  p.y - l.y;
 double dz =  p.z - l.z;

 return sqrt( dx * dx + dy * dy + dz * dz);
    

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

      
        std::size_t seed = 0;
        seed ^= hx + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hy + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hz + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};


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




void insert_location(const particle p) {

  
       Key key = generate_key(p.x, p.y, p.z);
       m[key].push_back(p);

 
 }


vector<particle> retrieve_neighbors(const particle p){
  double i,k,j;

  vector<particle> neighbors;
  
  double bound = 1/grid_size* grid_multiple; 

  for(i=p.x-bound; i < p.x+bound; i+=1/grid_size){

    for(k=p.y-bound; k < p.y+bound ; k+=1/grid_size){

      for(j=p.z-bound; j < p.z+bound; j+=1/grid_size){

	Key key = generate_key(i,k,j);


	if(m.find(key)!= m.end()){
	 
	  for(auto &n : m[key]){

	    if(get_distance(n,p)<=radius){

	      neighbors.push_back(n);
	    }

	  }
	  //neighbors.insert(neighbors.end(), m[key].begin(), m[key].end());

	}

      }
    }


  }
  return neighbors;

}

double get_value(int argc, char * argv[], const string c, const string error_message ){

  double result = -1;

  for(int i = 0; i < argc; i++){


    if(argv[i]==c){

      if(i+1<argc)
	{
	  result = stod(argv[i+1]);

	}

    }


  }
 

  if(result < 0){

    cout<<error_message<<endl;

    exit(EXIT_FAILURE);
  }


  return result; 

}


double get_grid_size(int argc, char *argv[]){


  return get_value(argc, argv, "-g","please insert -g parameter");

}

double get_radius(int argc, char *argv[]){


  return get_value(argc, argv, "-r","please insert -r parameter");

}




void parse_string(string line){

 

  while(line.find("</particle>")!=string::npos){
   

    particle particle ;
    middle_and_left_over data = get_middle(line,"<particle>", "</particle>");
    
    line = data.left_over; 
    
  
    string id = get_middle(data.middle, "<id>", "</id>").middle;
    if(!id.empty())
    particle.id = stoi(id);

    vector<double> position;
    vector<double> prev_position;
    vector<double> velocities;
    
    while(data.middle.find("<dimension>")!=string::npos){


      string value = get_middle(data.middle,"<dimension>","</dimension>").middle;

      if(!value.empty()){ 
        position.push_back(stod(value));
        data.middle = get_middle(data.middle,"<dimension>","</dimension>").left_over; 
       }
    }

     while(data.middle.find("<previous_dimension>")!=string::npos){


      string value = get_middle(data.middle,"<previous_dimension>","</previous_dimension>").middle;

      if(!value.empty()){ 
        prev_position.push_back(stod(value));
        data.middle = get_middle(data.middle,"<previous_dimension>","</previous_dimension>").left_over; 
       }
    }

    while(data.middle.find("<velocity>")!=string::npos){


      string value = get_middle(data.middle,"<velocity>","</velocity>").middle;

      if(!value.empty()){ 
        velocities.push_back(stod(value));
        data.middle = get_middle(data.middle,"<velocity>","</velocity>").left_over; 
       }
    }

     
    if(position.size() >= 3){
      particle.x = position.at(0);
      particle.y = position.at(1);
      particle.z = position.at(2); 

    }

    if(prev_position.size() >= 3){
      particle.prev_x = prev_position.at(0);
      particle.prev_y = prev_position.at(1);
      particle.prev_z = prev_position.at(2); 

    }

    if(velocities.size() >= 3){
      particle.vel_x = velocities.at(0);
      particle.vel_y = velocities.at(1);
      particle.vel_z = velocities.at(2); 

    }
    
    particles.push_back(particle);

  }
  

}



void read_data_from_stdin(){


  string line ;

   while (std::getline(std::cin, line)) {

     parse_string(line);

   }

}


void process(){


  for(const auto& p: particles){

    insert_location(p); 


  }

  printf("<particle_system>");
  for(const auto& p: particles){

    vector<particle>neighbors = retrieve_neighbors(p); 

    printf("<particle><id>%d</id><dimension>%lf</dimension><dimension>%lf</dimension><dimension>%lf</dimension>",
	   p.id, p.x, p.y, p.z);
    printf("<previous_dimension>%lf</previous_dimension><previous_dimension>%lf</previous_dimension><previous_dimension>%lf</previous_dimension>",
	   p.prev_x, p.prev_y, p.prev_z);
    printf("<velocity>%lf</velocity><velocity>%lf</velocity><velocity>%lf</velocity>",
	   p.vel_x, p.vel_y, p.vel_z);
    for(const auto &l: neighbors){
      printf("<neighbor_id>%d</neighbor_id>",l.id);
    }
    printf("</particle>");
  }

  printf("</particle_system>");


}


int main(int argc, char * argv[]){


  grid_size = get_grid_size(argc, argv);
  radius  = get_radius(argc, argv);
  
  read_data_from_stdin();
  process();
  

    
}
