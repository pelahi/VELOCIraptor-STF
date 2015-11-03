/* Tipsy header file which defines Tispy structures */

#ifndef TIPSY_STRUCTS_H
#define TIPSY_STRUCTS_H

struct dark_particle {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
} ;


struct star_particle {
    float mass;
    float pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
} ;

struct gas_particle {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
} ;

struct dump {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
} ;

struct simple_particle {
    float mass;
    float pos[3];
    float vel[3];
} ;

#endif
