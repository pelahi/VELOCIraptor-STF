/* Tipsy header file which defines Tispy structures */

#ifndef TIPSY_STRUCTS_H
#define TIPSY_STRUCTS_H

#include "endianutils.h"

struct dark_particle {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
    //assume that want big endian tipys format
    void SwitchtoBigEndian(void){
        mass=BigFloat(mass);
        eps=BigFloat(eps);
        phi=BigFloat(phi);
        for (int i=0;i<3;i++) {pos[i]=BigFloat(pos[i]);vel[i]=BigFloat(vel[i]);}
    }
} ;


struct star_particle {
    float mass;
    float pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
    void SwitchtoBigEndian(void){
        mass=BigFloat(mass);
        eps=BigFloat(eps);
        phi=BigFloat(phi);
        metals=BigFloat(metals);
        tform=BigFloat(tform);
        for (int i=0;i<3;i++) {pos[i]=BigFloat(pos[i]);vel[i]=BigFloat(vel[i]);}
    }
} ;

struct gas_particle {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float eps;
    float metals;
    float phi;
    void SwitchtoBigEndian(void){
        mass=BigFloat(mass);
        rho=BigFloat(rho);
        temp=BigFloat(temp);
        eps=BigFloat(eps);
        phi=BigFloat(phi);
        metals=BigFloat(metals);
        for (int i=0;i<3;i++) {pos[i]=BigFloat(pos[i]);vel[i]=BigFloat(vel[i]);}
    }
} ;

struct dump {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
    int pad;
    void SwitchtoBigEndian(void){
        time=BigDouble(time);
        nbodies=BigInt(nbodies);
        ndim=BigInt(ndim);
        nsph=BigInt(nsph);
        ndark=BigInt(ndark);
        nstar=BigInt(nstar);
        pad=BigInt(pad);
    }
    
} ;

struct simple_particle {
    float mass;
    float pos[3];
    float vel[3];
} ;

#endif
