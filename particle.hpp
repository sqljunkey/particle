#include<fstream>
#include<vector>
#include<random>
#include<filesystem>



using namespace std;



int interval=0;
vector<particle> particle_list;
double radius;
double particle_decay;
double speed;
double total_density;

double k = 0.012;
double k_near = 3.0;


void set_speed(double value){

  speed = value;
}


void set_decay(double value){
  particle_decay = value; 
}


void set_radius(double value){

  radius = value; 
}

void create_particles(int particle_num){

  random_device rd;
  mt19937 gen(rd());
  normal_distribution<float> d(0.0,1.0);

  
  for(int i=0; i<particle_num; i++){

    particle p;
    p.id = i;
    p.x=d(gen);
    p.y=d(gen);
    p.z=d(gen);


    p.x_momentum=0.0;
    p.y_momentum=0.0;
    p.z_momentum=0.0;

    particle_list.push_back(p);
  }
  
}


void grid_insert(){

  for(auto &p: particle_list ){

    insert_location(p);
     }
 

}

void clear_neighbors(){

  for(auto & p: particle_list){
    p.neighbors.clear();


  }

}

void update_neighbors(){

  if(interval>  10){
    clear_grid();  
    clear_neighbors();
    
    grid_insert();
    

    for(auto &p: particle_list){
      p.neighbors = retrieve_neighbors(p);
      
    }
        
   interval = 0;
 }
  
 interval++;

}

void apply_potential_layer(particle &p ,double pot, vector<double> v, double opacity){

  	p.x_momentum=p.x_momentum*particle_decay + pot*v.at(0)*0.01;
	p.y_momentum=p.y_momentum*particle_decay + pot*v.at(1)*0.01;
	p.z_momentum=p.z_momentum*particle_decay + pot*v.at(2)*0.01;

}


void centeral_velocity(double opacity){


  for(auto &p : particle_list){

    double pot = newton_potential(distance_from_center(p),1.0);
    
    vector<double> v = normalize(inverted_origin_vector(p));
    
    apply_potential_layer(p, pot, v, opacity);
    
  }

}

void double_easement_velocities(){

  for(auto &p : particle_list){
    double density = 0;
    double density_near = 0;
   
    for(auto &n : p.neighbors){

       double distance = get_distance(p,n);

       
      if(n.id !=p.id && distance <= radius ){
	
        double q = distance/radius; 
        double q_minus_one = (1-q);
	density+= q_minus_one*q_minus_one;
	density+= q_minus_one*q_minus_one*q_minus_one;
	
      }

    }

     double  pressure = k * (density - rest_density);
     double  pressure_near = k_near * density_near;
     double  particle_p_displacement = Vector2.Zero();


    for(auto &n : p.neighbors){

       double distance = get_distance(p,n);
       
      if(n.id !=p.id && distance <= radius ){

	double q = distance/radius; 
	vector<double> pn = normalize(get_vector(p, n));
        
	double displacement_term = pow (dt, 2) * (pressure * (1-q)
						   +pressure_near*pow(1-q,2));
	////scale position with displace term or whatever
	vector<double> d = scale(pn, displacement_term);

	swap_vector(particle[n.id], add(get_position_as_vector(particle[n.id]) ,scale(d,0.5)));  
	
      }

    }

    //add ppp

    
    
  }


}


void move_particles(){

  for(auto &p: particle_list){

    p.x= p.x + p.x_momentum * speed;
    p.y= p.y + p.y_momentum * speed;
    p.z= p.z + p.z_momentum * speed;

  }

}
void delete_file_before_write( string filename) {

 try {
   if (filesystem::remove(filename))
       cout << "file " << filename << " deleted.\n";
    else
       cout << "file " << filename << " not found.\n";
  }
 catch(const filesystem::filesystem_error& err) {
       cout << "filesystem error: " << err.what() << '\n';
  }



}

void write_particle_positions( string filename) {

  
  ofstream file(filename,std::ios_base::app);
  
    if (!file.is_open()) {
        cerr << "Error opening file " << filename << std::endl;
        return;
    }

  
        file << particle_list.size() << "\n\n";  

        for (const auto &particle : particle_list) {
            file << "BLUE" << " " 
                 << particle.x << " " 
                 << particle.y << " " 
                 << particle.z << "\n";
        }
    

    file.close();
}


void average_neighbors(){

  double total =0;
  for(const auto &p: particle_list){

   total+= p.neighbors.size();

  }

  cout<<"average particles: "<<total/particle_list.size()<<endl;

}

void particle_system(int iterations, int particle_num, string filename){

  delete_file_before_write(filename);
  create_particles(particle_num);
  

  for(int i=0; i < iterations; i++){

    cout<<"Iterations: "<<i<< endl;
    
   
    
    update_neighbors();
    average_neighbors();
    
    central_velocities();
    move_particles();
    write_particle_positions(filename);
    
  }

}

