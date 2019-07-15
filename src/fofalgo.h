/*! \file fofalgo.h
 *  \brief this file contains function declarations of fof algoritms
 */


#ifndef STFFOF_H
#define STFFOF_H

/// \name FOF algorithms
//@{
///stream FOF algorithm.
///Checks to see if particles have similar velocities and are linked in physical space.
///here param 6 is physical linking length, 7 is speed ratio  8 is min cosine of the velocity angle
int FOFStream(Particle &a, Particle &b, Double_t *params);
///Similar to \ref FOFStream but also requires potential to be outlier value and both particles must be above threshold param 9
int FOFStreamwithprob(Particle &a, Particle &b, Double_t *params);
///Similar to \ref FOFStreamwithprob but only one particle needs to be above threshold
int FOFStreamwithprobIterative(Particle &a, Particle &b, Double_t *params);
///similar to FOFStreamwithprob but for NN
int FOFStreamwithprobNN(Particle &a, Particle &b, Double_t *params);
///again similar to above but distance not checked
int FOFStreamwithprobNNNODIST(Particle &a, Particle &b, Double_t *params);
///similar to FOFStreamwithprob but the length scale is adjusted by the offset in velocities
int FOFStreamwithprobLX(Particle &a, Particle &b, Double_t *params);
///similar to above but for NN search
int FOFStreamwithprobNNLX(Particle &a, Particle &b, Double_t *params);
/// here similar but useful for AMR with linking length scaled
int FOFStreamwithprobscaleell(Particle &a, Particle &b, Double_t *params);
///similar to \ref FOFStreamwithprobscaleell but valid for nearest neighbour FOF search
int FOFStreamwithprobscaleellNN(Particle &a, Particle &b, Double_t *params);
///6dfof linking but with a potential checks to see if potential is below threshold value param 9
int FOF6dbg(Particle &a, Particle &b, Double_t *params);
/// optimised 6d FOF
int FOF6d_opt(Particle &a, Particle &b, Double_t *params);
///similar to \ref FOF6dbg but here particles have to be above threshold value
int FOF6dbgup(Particle &a, Particle &b, Double_t *params);
///checks to see if particles have positive types (useful for \ref GetVelocityDensity calculation with \ref STRUCDEN flag)
int FOFPositivetypes(Particle &a, Particle &b, Double_t *params);

///stream FOF algorithm that requires one of the particles to be dark matter for a link to occur
int FOF3dDM(Particle &a, Particle &b, Double_t *params);
//@}

/// \name FOF precheck algorithms
//@{
///checks to see if particle should be ingored based on potential being above threshold value stored in param 9
int FOFchecksub(Particle &a, Double_t *params);
///checks to see if particle should be ingored based on potential being below threshold value stored in param 9
int FOFcheckbg(Particle &a, Double_t *params);
///checks to see if particle appropriate type
int FOFchecktype(Particle &a, Double_t *params);
///checks to see if particle has a positive type
int FOFcheckpositivetype(Particle &a, Double_t *params);
//@}

#endif
