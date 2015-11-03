/*! \file DistFunc.h
 *  \brief header file for inline distance calculations
 */

#ifndef DISTFUNC_H
#define DISTFUNC_H
#include <NBodyMath.h>

namespace NBody
{
    /// Inline functions to get distances
    //@{
    /// Metric distance (squared).
    inline Double_t DistanceSqd(const Double_t *p1, const Double_t *p2, int dim)
    {
		Double_t total=0;
		for (int i=0;i<dim;i++) total+=(p1[i] - p2[i]) * (p1[i] - p2[i]);
		return total;
    }
    /// Velocity distance (squared).
    inline Double_t VelDistSqd(const Double_t *v1, const Double_t *v2, int dim)
    {
		Double_t total=0;
		for (int i=0;i<dim;i++) total+=(v1[i] - v2[i]) * (v1[i] - v2[i]);
		return total;
    }
    inline Double_t DistanceSqd(const Double_t *p1, const Double_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }
    /// Velocity distance (squared).
    inline Double_t VelDistSqd(const Double_t *v1, const Double_t *v2)
    {
        return (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    /// phase-space distance (squared).
    inline Double_t PhaseDistSqd(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]) +
            (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    /// Projected distance (squared).
    inline Double_t DistanceProjSqd(const Double_t *p1, const Double_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]);
    }

    /// like \ref DistanceProjSqd but with a metric.
    inline Double_t MetricDistSqd(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *metric)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++) total+=metric[i]*(coord1[i] - coord2[i]) * (coord1[i] - coord2[i]);
        return total;
    }
    /// distance within a given subset of dimensions (squared).
    inline Double_t MetricTensorDistSqd(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const GMatrix gm)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}

        for (int i=0;i<6;i++)
        for (int j=0;j<6;j++)
            total+=gm(i,j)*(coord1[i] - coord2[i]) * (coord1[j] - coord2[j]);
        return total;
    }
    /// distance within a given subset of dimensions using vector m0, m1, and matrix gm (squared). 
    /// Based on Enbid which discusses how to appropriate distance in arbitrary space with anisotropic smoothing (that is optimal local Mahalanobis metric) 
    inline Double_t MetricwithTensorDistSqd(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)
    {
        Double_t total=0.,temp;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++){
            temp=0.;
            for (int j=0;j<6;j++)
                temp+=gm(i,j)*(coord1[j] - coord2[j])*m0[j];
            temp*=m1[i];
            total+=temp*temp;
        }
        return total;
        //note that enbid uses the following to calculate distance, note here metric1 stores initial metric
        //and metric stoes sqrt(maxeigenvalue/(eigenvalues)) beloning to gmatrix
        /*for(k=0,pq->r=0.0; k<ND; k++)
        {
        for(l=0,temp=0; l<ND; ++l)
        {
            temp+=(pq->x[l]-searchcenter[l])*gmatrix[l][k]*metric1[l];
        }
        temp=temp*metric[k];
        pq->r+=temp*temp;
        }*/
    }
    //@}
    ///Inline functions with different types of pointers
    //@{
    inline Double_t DistanceSqd(const Real_t *p1, const Real_t *p2, int dim)
    {
        Double_t total=0;
        for (int i=0;i<dim;i++) total+=(p1[i] - p2[i]) * (p1[i] - p2[i]);
        return total;
    }
    inline Double_t VelDistSqd(const Real_t *v1, const Real_t *v2, int dim)
    {
        Double_t total=0;
        for (int i=0;i<dim;i++) total+=(v1[i] - v2[i]) * (v1[i] - v2[i]);
        return total;
    }
    inline Double_t DistanceSqd(const Real_t *p1, const Real_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }
    inline Double_t VelDistSqd(const Real_t *v1, const Real_t *v2)
    {
        return (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    inline Double_t PhaseDistSqd(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]) +
            (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    inline Double_t DistanceProjSqd(const Real_t *p1, const Real_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]);
    }
    inline Double_t MetricDistSqd(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *metric)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++) total+=metric[i]*(coord1[i] - coord2[i]) * (coord1[i] - coord2[i]);
        return total;
    }
    inline Double_t MetricTensorDistSqd(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const GMatrix gm)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}

        for (int i=0;i<6;i++)
        for (int j=0;j<6;j++)
            total+=gm(i,j)*(coord1[i] - coord2[i]) * (coord1[j] - coord2[j]);
        return total;
    }
    inline Double_t MetricwithTensorDistSqd(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)
    {
        Double_t total=0.,temp;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++){
            temp=0.;
            for (int j=0;j<6;j++)
                temp+=gm(i,j)*(coord1[j] - coord2[j])*m0[j];
            temp*=m1[i];
            total+=temp*temp;
        }
        return total;
    }
    
    inline Double_t DistanceSqd(const Double_t *p1, const Real_t *p2, int dim)
    {
        Double_t total=0;
        for (int i=0;i<dim;i++) total+=(p1[i] - p2[i]) * (p1[i] - p2[i]);
        return total;
    }
    inline Double_t VelDistSqd(const Double_t *v1, const Real_t *v2, int dim)
    {
        Double_t total=0;
        for (int i=0;i<dim;i++) total+=(v1[i] - v2[i]) * (v1[i] - v2[i]);
        return total;
    }
    inline Double_t DistanceSqd(const Double_t *p1, const Real_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }
    inline Double_t VelDistSqd(const Double_t *v1, const Real_t *v2)
    {
        return (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    inline Double_t PhaseDistSqd(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]) +
            (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    inline Double_t DistanceProjSqd(const Double_t *p1, const Real_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]);
    }
    inline Double_t MetricDistSqd(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *metric)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++) total+=metric[i]*(coord1[i] - coord2[i]) * (coord1[i] - coord2[i]);
        return total;
    }
    inline Double_t MetricTensorDistSqd(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const GMatrix gm)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}

        for (int i=0;i<6;i++)
        for (int j=0;j<6;j++)
            total+=gm(i,j)*(coord1[i] - coord2[i]) * (coord1[j] - coord2[j]);
        return total;
    }
    inline Double_t MetricwithTensorDistSqd(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)
    {
        Double_t total=0.,temp;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++){
            temp=0.;
            for (int j=0;j<6;j++)
                temp+=gm(i,j)*(coord1[j] - coord2[j])*m0[j];
            temp*=m1[i];
            total+=temp*temp;
        }
        return total;
    }
    
    inline Double_t DistanceSqd(const Real_t *p1, const Double_t *p2, int dim)
    {
        Double_t total=0;
        for (int i=0;i<dim;i++) total+=(p1[i] - p2[i]) * (p1[i] - p2[i]);
        return total;
    }
    inline Double_t VelDistSqd(const Real_t *v1, const Double_t *v2, int dim)
    {
        Double_t total=0;
        for (int i=0;i<dim;i++) total+=(v1[i] - v2[i]) * (v1[i] - v2[i]);
        return total;
    }
    inline Double_t DistanceSqd(const Real_t *p1, const Double_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]);
    }
    inline Double_t VelDistSqd(const Real_t *v1, const Double_t *v2)
    {
        return (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    inline Double_t PhaseDistSqd(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]) +
            (p1[2] - p2[2]) * (p1[2] - p2[2]) +
            (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]);
    }
    inline Double_t DistanceProjSqd(const Real_t *p1, const Double_t *p2)
    {
        return (p1[0] - p2[0]) * (p1[0] - p2[0]) +
            (p1[1] - p2[1]) * (p1[1] - p2[1]);
    }
    inline Double_t MetricDistSqd(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *metric)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++) total+=metric[i]*(coord1[i] - coord2[i]) * (coord1[i] - coord2[i]);
        return total;
    }
    inline Double_t MetricTensorDistSqd(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const GMatrix gm)
    {
        Double_t total=0;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}

        for (int i=0;i<6;i++)
        for (int j=0;j<6;j++)
            total+=gm(i,j)*(coord1[i] - coord2[i]) * (coord1[j] - coord2[j]);
        return total;
    }
    inline Double_t MetricwithTensorDistSqd(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)
    {
        Double_t total=0.,temp;
        Double_t coord1[6],coord2[6];
        for (int i=0;i<3;i++){coord1[i]=p1[i];coord1[i+3]=v1[i];coord2[i]=p2[i];coord2[i+3]=v2[i];}
        for (int i=0;i<6;i++){
            temp=0.;
            for (int j=0;j<6;j++)
                temp+=gm(i,j)*(coord1[j] - coord2[j])*m0[j];
            temp*=m1[i];
            total+=temp*temp;
        }
        return total;
    }
    //@}
    /// Periodic based functions
    //@{
    inline Double_t PeriodicReflection1D(const Coordinate &x0, Coordinate &xp, Coordinate p,int k){
        Double_t sval=0;
        xp=x0;
        if (x0[k]<p[k]/2.0) {xp[k]=x0[k]+p[k];sval=x0[k];}
        else {xp[k]=x0[k]-p[k];sval=-xp[k];}
        return sval;
    }
    /*
    inline double PeriodicReflectionDistance1D(Coordinate x0, Coordinate p, int k){
        if (x0[k]<p[k]/2.0) return=x0[k];
        else return p[k]-x0[k];
    }*/
    inline Double_t PeriodicReflection2D(const Coordinate &x0, Coordinate &xp, Coordinate p, int k1, int k2){
        Double_t sval=0;
        xp=x0;
        if (x0[k1]<p[k1]/2.0) {xp[k1]=x0[k1]+p[k1];sval+=x0[k1]*x0[k1];}
        else {xp[k1]=x0[k1]-p[k1];sval+=xp[k1]*xp[k1];}
        if (x0[k2]<p[k2]/2.0) {xp[k2]=x0[k2]+p[k2];sval+=x0[k2]*x0[k2];}
        else {xp[k2]=x0[k2]-p[k2];sval+=xp[k2]*xp[k2];}
        return sqrt(sval);
    }
    inline Double_t PeriodicReflectionND(const Coordinate &x0, Coordinate &xp, Coordinate p, int ndim){
        Double_t sval=0;
        xp=x0;
        for (int k=0;k<ndim;k++) {
            if (x0[k]<p[k]/2.0) {xp[k]=x0[k]+p[k];sval+=x0[k]*x0[k];}
            else {xp[k]=x0[k]-p[k];sval+=xp[k]*xp[k];}
        }
        return sqrt(sval);
    }

    inline Double_t PeriodicReflection1D(const Particle &x0, Particle &xp, Coordinate p,int k){
        Double_t sval=0;
        xp=x0;
        if (x0.GetPosition(k)<p[k]/2.0) {xp.SetPosition(k,x0.GetPosition(k)+p[k]);sval=x0.GetPosition(k);}
        else {xp.SetPosition(k,x0.GetPosition(k)-p[k]);sval=-xp.GetPosition(k);}
        return sval;
    }
    inline Double_t PeriodicReflection2D(const Particle &x0, Particle &xp, Coordinate p, int k1, int k2){
        Double_t sval=0;
        xp=x0;
        if (x0.GetPosition(k1)<p[k1]/2.0) {xp.SetPosition(k1,x0.GetPosition(k1)+p[k1]);sval+=x0.GetPosition(k1)*x0.GetPosition(k1);}
        else {xp.SetPosition(k1,x0.GetPosition(k1)-p[k1]);sval+=xp.GetPosition(k1)*xp.GetPosition(k1);}
        if (x0.GetPosition(k2)<p[k2]/2.0) {xp.SetPosition(k2,x0.GetPosition(k2)+p[k2]);sval+=x0.GetPosition(k2)*x0.GetPosition(k2);}
        else {xp.SetPosition(k2,x0.GetPosition(k2)-p[k2]);sval+=xp.GetPosition(k2)*xp.GetPosition(k2);}
        return sqrt(sval);
    }
    inline Double_t PeriodicReflectionND(const Particle &x0, Particle &xp, Coordinate p, int ndim){
        Double_t sval=0;
        xp=x0;
        for (int k=0;k<ndim;k++) {
            if (x0.GetPosition(k)<p[k]/2.0) {xp.SetPosition(k,x0.GetPosition(k)+p[k]);sval+=x0.GetPosition(k)*x0.GetPosition(k);}
            else {xp.SetPosition(k,x0.GetPosition(k)-p[k]);sval+=xp.GetPosition(k)*xp.GetPosition(k);}
        }
        return sqrt(sval);
    }

    //@}
}
#endif // DISTFUNC_H
