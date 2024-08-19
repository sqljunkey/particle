#include <GL/glew.h>
#include <GL/freeglut.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>

#define WIDTH 800
#define HEIGHT 600



int animation_position=0;
int number_of_frames =0;
int number_of_particles = 0;


GLfloat cameraX = 0.0f;
GLfloat cameraY = 0.0f;
GLfloat cameraZ = 5.0f;
GLfloat cameraAngleX = 0.0f;
GLfloat cameraAngleY = 0.0f;
GLfloat zoomStep = 0.1f; 


typedef struct{
  double x;
  double y;
  double z;
  double x_momentum;
  double y_momentum;
  double z_momentum;
  double mass;
  char *color;
}particle;


typedef struct{

  particle *particles;
   

}frame;          


void mouse_wheel_callback(int wheel, int direction, int x, int y) {
    
        if (direction > 0) {
            cameraZ -= zoomStep;
        } else {
            cameraZ += zoomStep;
        }
   
}




void key_callback(unsigned char key, int x, int y) {
    GLfloat moveSpeed = 0.1f; 

    switch (key) {
        case 'w': 
            cameraX -= moveSpeed * sin(cameraAngleY);
            cameraY += moveSpeed * sin(cameraAngleX);
            cameraZ -= moveSpeed * cos(cameraAngleY);
            break;
        case 's': 
            cameraX += moveSpeed * sin(cameraAngleY);
            cameraY -= moveSpeed * sin(cameraAngleX);
            cameraZ += moveSpeed * cos(cameraAngleY);
            break;
        case 'a': // Move left
            cameraX -= moveSpeed * cos(cameraAngleY);
            cameraZ -= moveSpeed * sin(cameraAngleY);
            break;
        case 'd': // Move right
            cameraX += moveSpeed * cos(cameraAngleY);
            cameraZ += moveSpeed * sin(cameraAngleY);
            break;
    }

    glutPostRedisplay();
}


void loadXYZ(const char* filename) {


  
  for(int i = 0; i < 9000; i++){
    

     points[i][0]=normal_distribution(0.0f,0.001f);
     points[i][1]=normal_distribution(0.0f,0.001f);
     points[i][2]=normal_distribution(0.0f,0.001f);

     
     numPoints++;
  }
   
    
    

    
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Failed to open XYZ file: %s\n", filename);
	 return;
    }

   
    GLfloat x, y, z;
    while (fscanf(file, "%f %f %f", &x, &y, &z) != EOF && numPoints < MAX_POINTS) {
        points[numPoints][0] = x;
        points[numPoints][1] = y;
        points[numPoints][2] = z;
        numPoints++;
    }

    fclose(file);

 
}


void drawPoints() {
    glBegin(GL_POINTS);
    for (int i = 0; i < numPoints; i++) {
        glVertex3f(points[i][0], points[i][1], points[i][2]);
    }
    glEnd();
}


void animatePoints() {

  Sleep((int)(1000.0f/66.66f));
 
  //fix
 
  glutPostRedisplay();   
}


void display() {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (GLfloat)WIDTH / (GLfloat)HEIGHT, 0.1f, 10000.0f);

   
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(cameraX, cameraY, cameraZ, 0.0f, 0.0f, 0.0f, 0.0f, 20.0f, 0.0f);

 
    glColor3f(1.0f, 1.0f, 1.0f);  
    drawPoints();

    
    animatePoints();

    
    glutSwapBuffers();

       
}


void init() {
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  
    glEnable(GL_DEPTH_TEST);               
    glEnable(GL_POINT_SMOOTH);            
}

int main(int argc, char** argv) {

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(WIDTH, HEIGHT);
    glutCreateWindow("OpenGL XYZ Loader");

   
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        fprintf(stderr, "Error initializing GLEW: %s\n", glewGetErrorString(err));
        return EXIT_FAILURE;
    }

 
    glutDisplayFunc(display);
    glutKeyboardFunc(key_callback);
    glutMouseWheelFunc(mouse_wheel_callback);

    loadXYZ("file.xyz");

   
    init();
    glutMainLoop();

    return EXIT_SUCCESS;
}
