#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libxml/parser.h>
#include <libxml/tree.h>

typedef struct {
    int *neighbor_ids;
    size_t num_neighbors;
    int id;
    double *dimensions;
    size_t num_dimensions;
    double *previous_dimensions;
    size_t num_previous_dimensions;
    double *velocities;
    size_t num_velocities;
} Particle;

typedef struct {
    Particle *particles;
    size_t num_particles;
} ParticleSystem;

 ParticleSystem particleSystem;
 double x_center, y_center , z_center;
double speed,cutoff; 

void parseParticle(xmlNode *node, Particle *particle) {
    xmlNode *cur_node = NULL;
    particle->num_dimensions = 0;
    particle->num_previous_dimensions = 0;
    particle->num_velocities = 0;
    particle->num_neighbors = 0;
    particle->dimensions = NULL;
    particle->previous_dimensions = NULL;
    particle->velocities = NULL;
    particle->neighbor_ids = NULL;

    for (cur_node = node->children; cur_node; cur_node = cur_node->next) {
        if (cur_node->type == XML_ELEMENT_NODE) {
            if (xmlStrcmp(cur_node->name, (const xmlChar *)"id") == 0) {
                particle->id = atoi((const char *)xmlNodeGetContent(cur_node));
            } else if (xmlStrcmp(cur_node->name, (const xmlChar *)"dimension") == 0) {
                particle->num_dimensions++;
                particle->dimensions = realloc(particle->dimensions, particle->num_dimensions * sizeof(double));
                particle->dimensions[particle->num_dimensions - 1] = atof((const char *)xmlNodeGetContent(cur_node));
            } else if (xmlStrcmp(cur_node->name, (const xmlChar *)"previous_dimension") == 0) {
                particle->num_previous_dimensions++;
                particle->previous_dimensions = realloc(particle->previous_dimensions, particle->num_previous_dimensions * sizeof(double));
                particle->previous_dimensions[particle->num_previous_dimensions - 1] = atof((const char *)xmlNodeGetContent(cur_node));
            } else if (xmlStrcmp(cur_node->name, (const xmlChar *)"velocity") == 0) {
                particle->num_velocities++;
                particle->velocities = realloc(particle->velocities, particle->num_velocities * sizeof(double));
                particle->velocities[particle->num_velocities - 1] = atof((const char *)xmlNodeGetContent(cur_node));
            } else if (xmlStrcmp(cur_node->name, (const xmlChar *)"neighbor_id") == 0) {
                particle->num_neighbors++;
                particle->neighbor_ids = realloc(particle->neighbor_ids, particle->num_neighbors * sizeof(int));
                particle->neighbor_ids[particle->num_neighbors - 1] = atoi((const char *)xmlNodeGetContent(cur_node));
            }
        }
    }
}

void parseParticleSystem(xmlNode *node, ParticleSystem *particleSystem) {
    xmlNode *cur_node = NULL;
    particleSystem->num_particles = 0;
    particleSystem->particles = NULL;

    for (cur_node = node; cur_node; cur_node = cur_node->next) {
        if (cur_node->type == XML_ELEMENT_NODE && xmlStrcmp(cur_node->name, (const xmlChar *)"particle") == 0) {
            particleSystem->num_particles++;
            particleSystem->particles = realloc(particleSystem->particles, particleSystem->num_particles * sizeof(Particle));
            parseParticle(cur_node, &particleSystem->particles[particleSystem->num_particles - 1]);
        }
    }
}

double newton_potential(double distance, double cutoff) {
  double newton = 0.0;

    if (distance > cutoff) {

	newton = 1.0 / pow(distance,2.0);

	}


  return newton;
}

double get_distance(double x, double y , double z, double xb, double yb, double zb){

  double dx = x-xb;
  double dy = y-yb;
  double dz = z-zb;
  
  return sqrt(dx*dx + dy*dy + dz*dz );
}

double* normalize(double x, double y, double z, double xb, double yb, double zb){

  double * normal = malloc(sizeof(double) * 3);



  double magnitude = get_distance(x, y, z, xb, yb, zb);

  
 normal[0] = x/magnitude;
 normal[1] = y/magnitude;
 normal[2] = z/magnitude;


  return normal;

}

double * get_vector(double x, double y, double z, double xb, double yb, double zb){

  double * v = malloc(sizeof(double)*3);

  v[0]=(-x+xb);
  v[1]=(-y+yb);
  v[2]=(-z+zb);

  return v;
}

void move_particle(){

    for (size_t i = 0; i < particleSystem.num_particles; i++) {

        Particle p = particleSystem.particles[i];

	double distance = get_distance(p.dimensions[0], p.dimensions[1], p.dimensions[2], x_center, y_center, z_center);
	double newton = newton_potential(distance, cutoff);
	double *v  = get_vector(p.dimensions[0], p.dimensions[1], p.dimensions[2], x_center, y_center, z_center); 
	double *normal = normalize(v[0],v[1],v[2],x_center, y_center, z_center);


       
	particleSystem.particles[i].velocities[0] += normal[0]*newton * speed;
        particleSystem.particles[i].velocities[1] += normal[1]*newton * speed;
	particleSystem.particles[i].velocities[2] += normal[2]*newton * speed;

	particleSystem.particles[i].previous_dimensions[0] = 	particleSystem.particles[i].dimensions[0];
	particleSystem.particles[i].previous_dimensions[1] = 	particleSystem.particles[i].dimensions[1];
	particleSystem.particles[i].previous_dimensions[2] = 	particleSystem.particles[i].dimensions[2];
	
	particleSystem.particles[i].dimensions[0] += particleSystem.particles[i].velocities[0];
	particleSystem.particles[i].dimensions[1] += particleSystem.particles[i].velocities[1];
	particleSystem.particles[i].dimensions[2] += particleSystem.particles[i].velocities[2];
         
	
    }

}


double get_value(int argc, char *argv[], const char * command , char * error_message ) {

  double result = -1.0;
  char *ptr;
  
  for(int i = 0; i < argc; i++){

    if( strstr(argv[i], command)){
      if(i+1<argc){
	result = strtod(argv[i+1], &ptr);
	break;
      }
    } 
  }

  if(result <0){

    printf(error_message);
    exit(1);

  }

  return result;

}


double get_center_x(int argc, char * argv[]){

  return get_value(argc, argv, "-x", "insert -x ( x center )into parameter");
  
}

double get_center_y(int argc, char * argv[]){

  return get_value(argc, argv, "-y", "insert -y ( y center )into parameter");
  
}
double get_center_z(int argc, char * argv[]){

  return get_value(argc, argv, "-z", "insert -z ( z center )into parameter");
  
}

double get_cutoff(int argc, char * argv[]){

  return get_value(argc, argv, "-c", "insert -c (newton center cuttoff ) into parameter");

}

double get_speed(int argc, char * argv[]){

  return get_value(argc, argv, "-s", "insert -s(speed) into parameter");

}

int main(int argc, char * argv[]) {


  x_center = get_center_x(argc, argv);
  y_center = get_center_y(argc, argv);
  z_center = get_center_z(argc, argv);
  cutoff = get_cutoff(argc, argv);
  speed = get_speed(argc, argv); 
   
    LIBXML_TEST_VERSION

  
    char *buffer = NULL;
    size_t size = 0;
    size_t capacity = 1024;
    buffer = (char *)malloc(capacity);
    if (!buffer) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    size_t read_size;
    while ((read_size = fread(buffer + size, 1, capacity - size, stdin)) > 0) {
        size += read_size;
        if (size >= capacity) {
            capacity *= 2;
            buffer = (char *)realloc(buffer, capacity);
            if (!buffer) {
                fprintf(stderr, "Memory reallocation failed\n");
                return 1;
            }
        }
    }
    buffer[size] = '\0'; 
    
    xmlDoc *doc = xmlReadMemory(buffer, size, "noname.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Failed to parse XML\n");
        free(buffer);
        return 1;
    }
    free(buffer);

    xmlNode *root_element = xmlDocGetRootElement(doc);

    
 
    parseParticleSystem(root_element->children, &particleSystem);


     
    move_particle();

    printf("<particle_system>");
    for (size_t i = 0; i < particleSystem.num_particles; i++) {
        Particle p = particleSystem.particles[i];
        printf("<particle><id>%d</id>", p.id);
        
        for (size_t j = 0; j < p.num_dimensions; j++) {
            printf("<dimensions>%lf</dimensions>", p.dimensions[j]);
        }
        
        for (size_t j = 0; j < p.num_previous_dimensions; j++) {
            printf("<previous_dimension>%lf</previous_dimension>", p.previous_dimensions[j]);
        }
      
        for (size_t j = 0; j < p.num_velocities; j++) {
            printf("<velocity>%lf</velocity>", p.velocities[j]);
        }
       
        for (size_t j = 0; j < p.num_neighbors; j++) {
            printf("<neighbor_id>%d</neighbor_id>", p.neighbor_ids[j]);
        }
        printf("</particle>");
    }
     printf("</particle_system>");
   
    for (size_t i = 0; i < particleSystem.num_particles; i++) {
        free(particleSystem.particles[i].dimensions);
        free(particleSystem.particles[i].previous_dimensions);
        free(particleSystem.particles[i].velocities);
        free(particleSystem.particles[i].neighbor_ids);
    }
    free(particleSystem.particles);

  
    xmlFreeDoc(doc);

   
    xmlCleanupParser();

    return 0;
}
