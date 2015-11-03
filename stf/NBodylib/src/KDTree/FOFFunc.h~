/*! \file FOFFunc.h
 *  \brief header file for inline fof comparison functions
    
    A \b FOF comparison function is of a particular form:\n
    NOTE: parameters for tree search and comparison are contained in the \b \e params variable
    \arg The first 5 indices of the params variable indicating how tree is searched. 
    \arg The [0] address is used to store treetype. The second value and third value are the 
    \arg associated distance^2 measures for physical and velocity space, the next 3 spaces are 
    left for possible changes. \n 
    The remaining parameters used in comparison function start at params[6] & the use depends on the comparison function.\n

    NOTE: comparison function is used only for leaf nodes. 
 */

#ifndef FOFCOMPFUNC_H
#define FOFCOMPFUNC_H

namespace NBody
{

    //-- FOF Comparison Functions

    ///FOFcompfunc type
    typedef int (*FOFcompfunc) (Particle&, Particle&, Double_t *);
    ///FOFcheckfunc type where return -1 if no need to check the particle in fof search
    typedef int (*FOFcheckfunc) (Particle&, Double_t *);

    // inline FOF algorithms

    inline int FOF3d(Particle &a, Particle &b, Double_t *params){
        Double_t total=0;
        for (int j=0;j<3;j++)
            total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
        return (total<1);
    }

    //velocity fof
    //here params 6 velocity linking lengths
    inline int FOFVel(Particle &a, Particle &b, Double_t *params){
        Double_t total=0;
        for (int j=0;j<3;j++)
            total+=(a.GetVelocity(j)-b.GetVelocity(j))*(a.GetVelocity(j)-b.GetVelocity(j))/params[6];
        return (total<1);
    }

    //simple 6D fof
    //here params 6 and 7 are physical and velocity linking lengths
    inline int FOF6d(Particle &a, Particle &b, Double_t *params){
        Double_t total=0;
        for (int j=0;j<3;j++){
            total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
            total+=(a.GetVelocity(j)-b.GetVelocity(j))*(a.GetVelocity(j)-b.GetVelocity(j))/params[7];
        }
        return (total<1);
    }
    //simple no check function just as placeholder
    inline int Pnocheck(Particle &a, Double_t *params){return 0;}
}
#endif 
