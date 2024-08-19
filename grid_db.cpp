#include<iostream>
#include "particle_struct.hpp"
#include "grid_db.hpp"

#include "particle_math.hpp"
#include "particle.hpp"



int main(int argc, char * argv[] ) {

  double animation = 0;
  double particle =0;

  for(int i = 0;i < argc; i++){

  
  string c =   argv[i];  

  double value =0.0;
  

  if(c == "-g"){
            cout<<c<<": "<<argv[i+1]<<endl;
            value = stod(argv[i+1]);
            set_gridsize(value);
	    i++;
  }
  else if (c== "-gm"){
            cout<<c<<": "<<argv[i+1]<<endl;
	    value = stod(argv[i+1]);
            set_grid_multiple(value);
	    i++;
  }
  else if (c== "-s"){
            cout<<c<<": "<<argv[i+1]<<endl;
            value = stod(argv[i+1]);
            set_speed(value);
	    i++;
  }
  else if (c == "-r"){
            cout<<c<<": "<<argv[i+1]<<endl;
            value = stod(argv[i+1]);
            set_radius(value);
	    i++;
  }
  else if (c== "-d"){
            cout<<c<<": "<<argv[i+1]<<endl;
            value = stod(argv[i+1]);
            set_decay(value);
	    i++;
  }
  else if (c== "-p"){
            cout<<c<<": "<<argv[i+1]<<endl;
            particle = stod(argv[i+1]);
            set_decay(value);
	    i++;
  }

  else if (c== "-a"){
            cout<<c<<": "<<argv[i+1]<<endl;
            animation = stod(argv[i+1]);
            set_decay(value);
	    i++;
  }
  
  else {
    cout<< "-a : animation"<< endl;
    cout<< "-p : particle"<<endl;
    cout<< "-d : decay"<<endl;
    cout<< "-r : radius"<<endl;
    cout<< "-gm : grid multiple"<<endl;
    cout<< "-g : grid size"<<endl;
    cout<< "-s : speed" << endl;
    cout<< "Usage: -d 12 -r 23 -s 12 -gm 12 -g 12 "<<endl;
    
  }
}

  cout<<"multiple "<<grid_multiple<<endl;
  
  particle_system(animation,particle,"file.xyz");

  return 0;
}
