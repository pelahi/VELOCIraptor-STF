/*! \file fofalgo.cxx
 *  \brief this file contains function declarations of fof algoritms
 */

#include "stf.h"

int FOFStream(Particle &a, Particle &b, Double_t *params){
    Double_t total=0,v1=0,v2=0,vdot=0,vnorm;
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    vdot*=vnorm;
    return (total<1.0&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOFStreamwithprob(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t total=0,v1=0,v2=0,vdot=0,vnorm;
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    vdot*=vnorm;
    return (total<1.0&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOFStreamwithprobIterative(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]&&b.GetPotential()<params[9]) return 0;
    Double_t total=0,v1=0,v2=0,vdot=0,vnorm;
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    vdot*=vnorm;
    return (total<1&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}


int FOFStreamwithprobNN(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t total=0,v1=0,v2=0,vdot=0,vnorm;
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    vdot*=vnorm;
    return (total<=1.0&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOFStreamwithprobNNNODIST(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t v1=0,v2=0,vdot=0,vnorm;
    for (int j=0;j<3;j++){
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    vdot*=vnorm;
    return (vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOFStreamwithprobLX(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t ds1=0,ds2=0,v1=0,v2=0,vdot=0,vnorm,total=0;
    for (int j=0;j<3;j++){
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    for (int j=0;j<3;j++){
        ds1+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/(params[6]*0.25*(1.0+a.GetVelocity(j)*a.GetVelocity(j)/v1)*(1.0+a.GetVelocity(j)*a.GetVelocity(j)/v1));
        ds2+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/(params[6]*0.25*(1.0+b.GetVelocity(j)*b.GetVelocity(j)/v2)*(1.0+b.GetVelocity(j)*b.GetVelocity(j)/v2));
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/v1/v2;
    total=min(ds1,ds2);
    vdot*=vnorm;
    return (total<=1.0&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOFStreamwithprobNNLX(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t dist2=params[6];
    Double_t ds1=0,ds2=0,v1=0,v2=0,vdot=0,vnorm,total=0;
    for (int j=0;j<3;j++){
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    for (int j=0;j<3;j++){
        ds1+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/(dist2*0.25*(1.0+a.GetVelocity(j)*a.GetVelocity(j)/v1)*(1.0+a.GetVelocity(j)*a.GetVelocity(j)/v1));
        ds2+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/(dist2*0.25*(1.0+b.GetVelocity(j)*b.GetVelocity(j)/v2)*(1.0+b.GetVelocity(j)*b.GetVelocity(j)/v2));
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    total=min(ds1,ds2);
    vdot*=vnorm;
    return (total<=1.0&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOFStreamwithprobscaleell(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t total=0,v1=0,v2=0,vdot=0,vnorm,ellscale;
    //use smaller ellscale where linking length scales with mass relative to smallest mass params[10]
    if (a.GetMass()<=b.GetMass()) ellscale=params[6]*pow(a.GetMass()/params[10],(Double_t)(2./3.0));
    else ellscale=params[6]*pow(b.GetMass()/params[10],(Double_t)(2./3.0));
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/ellscale;
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    vdot*=vnorm;
    return (total<1&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOFStreamwithprobscaleellNN(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t total=0,v1=0,v2=0,vdot=0,vnorm,ellscale;
    //use smaller ellscale where linking length scales with mass relative to smallest mass params[10]
    if (a.GetMass()<=b.GetMass()) ellscale=params[6]*pow(2.1,2.0*log(a.GetMass()/params[10])/log(2.0)/3.0);
    else ellscale=params[6]*pow(2.1,2.0*log(a.GetMass()/params[10])/log(2.0)/3.0);
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/ellscale;
        v1+=a.GetVelocity(j)*a.GetVelocity(j);
        v2+=b.GetVelocity(j)*b.GetVelocity(j);
        vdot+=a.GetVelocity(j)*b.GetVelocity(j);
    }
    v1=sqrt(v1);v2=sqrt(v2);
    vnorm=1.0/(v1*v2);
    vdot*=vnorm;
    return (total<1&&vdot>params[8]&&v1/v2<params[7]&&v1/v2>1.0/params[7]);
}

int FOF6dbg(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()>=params[9]||b.GetPotential()>=params[9]) return 0;
    Double_t total=0;
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
        total+=(a.GetVelocity(j)-b.GetVelocity(j))*(a.GetVelocity(j)-b.GetVelocity(j))/params[7];
    }
    return (total<1);
}
int FOF6dbgup(Particle &a, Particle &b, Double_t *params){
    if (a.GetPotential()<params[9]||b.GetPotential()<params[9]) return 0;
    Double_t total=0;
    for (int j=0;j<3;j++){
        total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
        total+=(a.GetVelocity(j)-b.GetVelocity(j))*(a.GetVelocity(j)-b.GetVelocity(j))/params[7];
    }
    return (total<1);
}

/// Optimised version of 6DFOF that performs no divisions
int FOF6d_opt(Particle &a, Particle &b, Double_t *params){
    Double_t total_x=0;
    Double_t total_v=0;
    for (int j=0;j<3;j++){
        total_x+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j));
        total_v+=(a.GetVelocity(j)-b.GetVelocity(j))*(a.GetVelocity(j)-b.GetVelocity(j));
    }
    total_x*=params[7];
    total_v*=params[6];
    return (total_x + total_v < params[7] * params[6]);
}

///stream FOF algorithm that requires the primary particle to be dark matter for a link to occur
int FOF3dDM(Particle &a, Particle &b, Double_t *params){
    if (a.GetType()!=int(params[7])) return 0;
    Double_t total=0;
    for (int j=0;j<3;j++) total+=(a.GetPosition(j)-b.GetPosition(j))*(a.GetPosition(j)-b.GetPosition(j))/params[6];
    return (total<1);
}

int FOFPositivetypes(Particle &a, Particle &b, Double_t *params){
    return (a.GetType()>=0 && b.GetType()>=0);
}


int FOFchecksub(Particle &a, Double_t *params){
    if (a.GetPotential()>=params[9]) return 0;
    else return -1;
}
int FOFcheckbg(Particle &a, Double_t *params){
    if (a.GetPotential()<params[9]) return 0;
    else return -1;
}
int FOFchecktype(Particle &a, Double_t *params){
    if (a.GetType()==int(params[7])) return 0;
    else return -1;
}
int FOFcheckpositivetype(Particle &a, Double_t *params){
    if (a.GetType()>=0) return 0;
    else return -1;
}

//@}
