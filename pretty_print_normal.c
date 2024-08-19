#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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




double distance(double x, double y, double z){



}

void move_particles(){

}

int main(void) {
    // Initialize the library and check potential ABI mismatches
    LIBXML_TEST_VERSION

    // Read the entire stdin into a buffer
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
    buffer[size] = '\0'; // Null-terminate the buffer

    
    xmlDoc *doc = xmlReadMemory(buffer, size, "noname.xml", NULL, 0);
    if (doc == NULL) {
        fprintf(stderr, "Failed to parse XML\n");
        free(buffer);
        return 1;
    }
    free(buffer);

    xmlNode *root_element = xmlDocGetRootElement(doc);

    
    ParticleSystem particleSystem;
    parseParticleSystem(root_element->children, &particleSystem);

    // Print the parsed data
    for (size_t i = 0; i < particleSystem.num_particles; i++) {
        Particle p = particleSystem.particles[i];
        printf("Particle ID: %d\n", p.id);
        printf("Dimensions: ");
        for (size_t j = 0; j < p.num_dimensions; j++) {
            printf("%f ", p.dimensions[j]);
        }
        printf("\nPrevious Dimensions: ");
        for (size_t j = 0; j < p.num_previous_dimensions; j++) {
            printf("%f ", p.previous_dimensions[j]);
        }
        printf("\nVelocities: ");
        for (size_t j = 0; j < p.num_velocities; j++) {
            printf("%f ", p.velocities[j]);
        }
        printf("\nNeighbors: ");
        for (size_t j = 0; j < p.num_neighbors; j++) {
            printf("%d ", p.neighbor_ids[j]);
        }
        printf("\n");
    }

   
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
