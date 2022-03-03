#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "graphics/graphics.h"
#include <omp.h>

typedef struct {
  double x, y, mass, vel_x, vel_y, brightness, f_x, f_y; 
} particle;


void runSimulation(int N, particle* particles, double delta_t, int nsteps, int useGraphics, char* argv[]) {
  // define constants
  const double smoothing = 1e-3;
  double G = 100.0 / (double)N;

  const int WIDTH = 800;
  const int HEIGHT = 800;

  if(useGraphics) {
    InitializeGraphics(argv[0],WIDTH,HEIGHT);
    SetCAxes(0,1);
  }

  // Run simulation for nsteps
  for (int step = 0; step < nsteps; step++) {
    // Reset force accumulators to zero
    for (int i = 0; i < N; i++) {
      particles[i].f_x = 0;
      particles[i].f_y = 0;
    }
    // Calculate the forces on all particles
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
      particle p_i = particles[i];
      for (int j = 0; j < N; j++) {
        particle p_j = particles[j];
        // Force calculations
        double diff_x = p_i.x - p_j.x;
        double diff_y = p_i.y - p_j.y;
        double denom = (sqrt(diff_x * diff_x + diff_y * diff_y) + smoothing);
        denom = denom*denom*denom;
        
        double k = -G * p_i.mass * p_j.mass / denom;

        double f_x = k*diff_x;
        double f_y = k*diff_y;

        particles[i].f_x += f_x;
        particles[i].f_y += f_y;
        
        //particles[j].f_x -= f_x;
        //particles[j].f_y -= f_y;
      }
    }

    // Perform Symplectic Euler step
    for (int i = 0; i < N; i++) {
      double acc_x = particles[i].f_x / particles[i].mass;
      double acc_y = particles[i].f_y / particles[i].mass;
      particles[i].vel_x += delta_t * acc_x; 
      particles[i].vel_y += delta_t * acc_y; 
      particles[i].x += delta_t * particles[i].vel_x;
      particles[i].y += delta_t * particles[i].vel_y;
    }

    // Show graphics
    if (useGraphics) {
      ClearScreen();
      for (int i = 0; i < N; i++) {
        DrawCircle(particles[i].x, particles[i].y, 1, 1, 0.005*particles[i].mass, particles[i].brightness/5);
      }
      Refresh();
      if(CheckForQuit())
        break;
    }

  }
  
  // Keep the last frame active until user hits q
  if(useGraphics) {
    printf("Done.");
    while(!CheckForQuit()) {
      Refresh();
      usleep(300);
    }
  }
}


int main(int argc, char* argv[]) {
  if (argc != 7) {
    printf("Exactly 5 parameters should be given like this: ./galsim N filename nsteps delta_t graphics\n");
    return 1;
  }

  // Read in the cmdline arguments
  int N, nsteps, useGraphics, nThreads;
  char* filename;
  double delta_t;
  N = atoi(argv[1]);
  filename = argv[2];
  nsteps = atoi(argv[3]);
  delta_t = atof(argv[4]);
  useGraphics = atoi(argv[5]);
  nThreads = atoi(argv[6]);

  omp_set_num_threads(nThreads);
  // Allocate a list to store all particles
  particle* particles = (particle*) malloc(N*sizeof(particle));

  // Open file and raise error if it doesn't exists 
  FILE* f = fopen(filename, "r"); 
  if (f == NULL) {
    printf("file %s doesn't exists!", filename);
    return 1;
  } 

  // Read in all particles
  for (int i = 0; i < N; i++) {
    fread(&(particles[i].x), sizeof(double), 1, f);     
    fread(&(particles[i].y), sizeof(double), 1, f);     
    fread(&(particles[i].mass), sizeof(double), 1, f);     
    fread(&(particles[i].vel_x), sizeof(double), 1, f);     
    fread(&(particles[i].vel_y), sizeof(double), 1, f);     
    fread(&(particles[i].brightness), sizeof(double), 1, f);
  }
  fclose(f);

  runSimulation(N, particles, delta_t, nsteps, useGraphics, argv);

  // Open the result file 
  f = fopen("result.gal", "w");
  
  // Write the list of particles to file
  for (int i = 0; i < N; i++) {
    fwrite(&(particles[i].x), sizeof(double), 1, f);     
    fwrite(&(particles[i].y), sizeof(double), 1, f);     
    fwrite(&(particles[i].mass), sizeof(double), 1, f);     
    fwrite(&(particles[i].vel_x), sizeof(double), 1, f);     
    fwrite(&(particles[i].vel_y), sizeof(double), 1, f);     
    fwrite(&(particles[i].brightness), sizeof(double), 1, f);     
    printf("x, y, mass: %f, %f, %f\n", particles[i].x, particles[i].y, particles[i].mass);
  } 
  fclose(f);

  free(particles);
  return 0;
}
