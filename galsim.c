#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "graphics/graphics.h"

typedef struct {
  double x, y, mass, vel_x, vel_y, brightness; 
} particle;


void runSimulation(int N, particle* particles, double delta_t, int nsteps, int useGraphics, char* argv[]) {
  const double smoothing = 1e-3;
  double* force = (double*)malloc(sizeof(double)*2*N);
  double G = 100.0 / (double)N;

  const int WIDTH = 800;
  const int HEIGHT = 800;

  if(useGraphics) {
    InitializeGraphics(argv[0],WIDTH,HEIGHT);
    SetCAxes(0,1);
  }

  for (int step = 0; step < nsteps; step++) {
    for (int i = 0; i < N; i++) {
      force[2*i] = 0.0;
      force[2*i+1] = 0.0;
      for (int j = 0; j < N; j++) {
        if (i != j) {
          double diff_x = particles[i].x - particles[j].x;
          double diff_y = particles[i].y - particles[j].y;
          double r = sqrt(diff_x * diff_x + diff_y * diff_y);
          double denom = (r + smoothing) * (r + smoothing) * (r + smoothing);
          
          force[i*2] += -G * particles[i].mass * particles[j].mass * diff_x / denom;   
          force[i*2+1] += -G * particles[i].mass * particles[j].mass * diff_y / denom;   
        }
      }
    }

    for (int i = 0; i < N; i++) {
      double acc_x = force[2*i] / particles[i].mass;
      double acc_y = force[2*i+1] / particles[i].mass;
      particles[i].vel_x += delta_t * acc_x; 
      particles[i].vel_y += delta_t * acc_y; 
      particles[i].x += delta_t * particles[i].vel_x;
      particles[i].y += delta_t * particles[i].vel_y;
    }

    if (useGraphics) {
      ClearScreen();
      for (int i = 0; i < N; i++) {
        DrawCircle(particles[i].x, particles[i].y, 1, 1, 0.003, 0.5);
      }
      Refresh();
    }

  } 

  if(useGraphics) {
    printf("Done.");
    while(!CheckForQuit()) {
      Refresh();
      usleep(300);
    }
  }
}


int main(int argc, char* argv[]) {
  if (argc != 6) {
    printf("Exactly 5 parameters should be given!\n");
    return 1;
  }

  int N, nsteps, useGraphics;
  char* filename;
  double delta_t;
  N = atoi(argv[1]);
  filename = argv[2];
  nsteps = atoi(argv[3]);
  delta_t = atof(argv[4]);
  useGraphics = atoi(argv[5]);
  printf("N = %d, filename = %s, nsteps = %d, delta_t = %f, useGraphics = %d\n", N, filename, nsteps, delta_t, useGraphics);

  particle* particles = (particle*) malloc(N*sizeof(particle));

  FILE* f = fopen(filename, "r"); 
  if (f == NULL) {
    printf("file %s doesn't exists!", filename);
    return 1;
  } 
  for (int i = 0; i < N; i++) {
    fread(&(particles[i].x), sizeof(double), 1, f);     
    fread(&(particles[i].y), sizeof(double), 1, f);     
    fread(&(particles[i].mass), sizeof(double), 1, f);     
    fread(&(particles[i].vel_x), sizeof(double), 1, f);     
    fread(&(particles[i].vel_y), sizeof(double), 1, f);     
    fread(&(particles[i].brightness), sizeof(double), 1, f);     
    printf("x, y, brightness: %f, %f, %f\n", particles[i].x, particles[i].y, particles[i].brightness);
  }
  fclose(f);

  runSimulation(N, particles, delta_t, nsteps, useGraphics, argv);

  f = fopen("result.gal", "w");
  
  for (int i = 0; i < N; i++) {
    fwrite(&(particles[i].x), sizeof(double), 1, f);     
    fwrite(&(particles[i].y), sizeof(double), 1, f);     
    fwrite(&(particles[i].mass), sizeof(double), 1, f);     
    fwrite(&(particles[i].vel_x), sizeof(double), 1, f);     
    fwrite(&(particles[i].vel_y), sizeof(double), 1, f);     
    fwrite(&(particles[i].brightness), sizeof(double), 1, f);     
    printf("x, y, brightness: %f, %f, %f\n", particles[i].x, particles[i].y, particles[i].brightness);
  } 
  fclose(f);

  free(particles);
  return 0;
}
