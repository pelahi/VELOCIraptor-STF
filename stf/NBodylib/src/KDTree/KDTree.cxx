/*! \file KDTree.cxx
 *  \brief This file contains subroutines involving tree construction

    NOTE: when possible, openmp parallelism implemented otherwise have MPI or simple serial code
    NOTE: openmp reductions are not implemented for min, max searches as not part of OpenMP API yet.
    and realistically, because one would have to implement a critical stop to check global minimum
    not really worth it. With MPI, obviously could sent bucket list to appropriate threads and min/max
    for that thread and do global broadcast to find global min/max but again overhead means it might not be
    worth it.
    I could also generate a array of size numthreads and then find min/max in that thread
    Must implement check that size of array worthwhile running in parallel
    \todo I can also reduce tree build time related to select the median point. A simple heuristic to avoid coding a complex linear-time median-finding algorithm, or using an O(n log n) sort of all n points, is to use sort to find the median of a fixed number of randomly selected points to serve as the splitting plane. In practice, this technique often results in nicely balanced trees.

*/

#include <KDTree.h>

namespace NBody
{

    // -- Inline functions that get called often when building the tree.

    /// \name Find most spread dimension
    //@{
    inline Double_t KDTree::SpreadestPos(int j, Int_t start, Int_t end, Double_t *bnd)
    {
        Double_t min = bucket[start].GetPosition(j);
        Double_t max = min;
        Int_t i;
#ifndef USEOPENMP
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPosition(j) < min) min = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > max) max = bucket[i].GetPosition(j);
        }
#else
        if (end-start<CRITPARALLELSIZE){
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPosition(j) < min) min = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > max) max = bucket[i].GetPosition(j);
        }
        }
        else {
        int nthreads;
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
        Double_t *mina=new Double_t[nthreads];
        Double_t *maxa=new Double_t[nthreads];
        for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=min;
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic) nowait
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPosition(j) < mina[omp_get_thread_num()]) mina[omp_get_thread_num()] = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > maxa[omp_get_thread_num()]) maxa[omp_get_thread_num()] = bucket[i].GetPosition(j);
        }
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < min) min = mina[i];
            if (maxa[i] > max) max = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
        }
#endif
        bnd[0]=min;bnd[1]=max;
        return max - min;
    }
    inline Double_t KDTree::SpreadestVel(int j, Int_t start, Int_t end, Double_t *bnd)
    {
        Double_t min = bucket[start].GetVelocity(j);
        Double_t max = min;
        Int_t i;
#ifndef USEOPENMP
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetVelocity(j) < min) min = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > max) max = bucket[i].GetVelocity(j);
        }
#else
        if (end-start<CRITPARALLELSIZE){
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetVelocity(j) < min) min = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > max) max = bucket[i].GetVelocity(j);
        }
        }
        else {
        int nthreads;
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
        Double_t *mina=new Double_t[nthreads];
        Double_t *maxa=new Double_t[nthreads];
        for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=min;
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic) nowait
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetVelocity(j) < mina[omp_get_thread_num()]) mina[omp_get_thread_num()] = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > maxa[omp_get_thread_num()]) maxa[omp_get_thread_num()] = bucket[i].GetVelocity(j);
        }
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < min) min = mina[i];
            if (maxa[i] > max) max = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
        }
#endif
        bnd[0]=min;bnd[1]=max;
        return max - min;
    }
    inline Double_t KDTree::SpreadestPhs(int j, Int_t start, Int_t end, Double_t *bnd)
    {
        Double_t min = bucket[start].GetPhase(j);
        Double_t max = min;
        Int_t i;
#ifndef USEOPENMP
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPhase(j) < min) min = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > max) max = bucket[i].GetPhase(j);
        }
#else
        if (end-start<CRITPARALLELSIZE){
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPhase(j) < min) min = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > max) max = bucket[i].GetPhase(j);
        }
        }
        else {
        int nthreads;
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
        Double_t *mina=new Double_t[nthreads];
        Double_t *maxa=new Double_t[nthreads];
        for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=min;
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic) nowait
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPhase(j) < mina[omp_get_thread_num()]) mina[omp_get_thread_num()] = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > maxa[omp_get_thread_num()]) maxa[omp_get_thread_num()] = bucket[i].GetPhase(j);
        }
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < min) min = mina[i];
            if (maxa[i] > max) max = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
        }
#endif
        bnd[0]=min;bnd[1]=max;
        return max - min;
    }
    //@}
    /// \name Find the boundary of the data and return mean
    //@{
    inline Double_t KDTree::BoundaryandMeanPos(int j, Int_t start, Int_t end, Double_t *bnd)
    {
        Double_t mean=bucket[start].GetPosition(j);
        bnd[0] = mean;
        bnd[1] = bnd[0];
        Int_t i;
#ifndef USEOPENMP
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPosition(j) < bnd[0]) bnd[0] = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > bnd[1]) bnd[1] = bucket[i].GetPosition(j);
            mean+=bucket[i].GetPosition(j);
        }
#else
        if (end-start<CRITPARALLELSIZE){
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPosition(j) < bnd[0]) bnd[0] = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > bnd[1]) bnd[1] = bucket[i].GetPosition(j);
            mean+=bucket[i].GetPosition(j);
        }
        }
        else {
        int nthreads;
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
        Double_t *mina=new Double_t[nthreads];
        Double_t *maxa=new Double_t[nthreads];
        for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=bnd[0];
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic) nowait
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPosition(j) < mina[omp_get_thread_num()]) mina[omp_get_thread_num()] = bucket[i].GetPosition(j);
            if (bucket[i].GetPosition(j) > maxa[omp_get_thread_num()]) maxa[omp_get_thread_num()] = bucket[i].GetPosition(j);
        }
    #pragma omp for reduction(+:mean)
        for (i = start+1; i < end; i++) mean+=bucket[i].GetPosition(j);
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < bnd[0]) bnd[0] = mina[i];
            if (maxa[i] > bnd[1]) bnd[1] = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
        }
#endif
        mean/=(Double_t)(end-start);
        return mean;
    }
    inline Double_t KDTree::BoundaryandMeanVel(int j, Int_t start, Int_t end, Double_t *bnd)
    {
        Double_t mean=bucket[start].GetVelocity(j);
        bnd[0] = bnd[1] = mean;
        Int_t i;
#ifndef USEOPENMP
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetVelocity(j) < bnd[0]) bnd[0] = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > bnd[1]) bnd[1] = bucket[i].GetVelocity(j);
            mean+=bucket[i].GetVelocity(j);
        }
#else
        if (end-start<CRITPARALLELSIZE){
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetVelocity(j) < bnd[0]) bnd[0] = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > bnd[1]) bnd[1] = bucket[i].GetVelocity(j);
            mean+=bucket[i].GetVelocity(j);
        }
        }
        else {
        int nthreads;
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
        Double_t *mina=new Double_t[nthreads];
        Double_t *maxa=new Double_t[nthreads];
        for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=bnd[0];
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic) nowait
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetVelocity(j) < mina[omp_get_thread_num()]) mina[omp_get_thread_num()] = bucket[i].GetVelocity(j);
            if (bucket[i].GetVelocity(j) > maxa[omp_get_thread_num()]) maxa[omp_get_thread_num()] = bucket[i].GetVelocity(j);
        }
    #pragma omp for reduction(+:mean)
        for (i = start+1; i < end; i++) mean+=bucket[i].GetVelocity(j);
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < bnd[0]) bnd[0] = mina[i];
            if (maxa[i] > bnd[1]) bnd[1] = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
        }
#endif
        mean/=(Double_t)(end-start);
        return mean;
    }
    inline Double_t KDTree::BoundaryandMeanPhs(int j, Int_t start, Int_t end, Double_t *bnd)
    {
        Double_t mean=bucket[start].GetPhase(j);
        bnd[0] = bnd[1] = mean;
        Int_t i;
#ifndef USEOPENMP
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPhase(j) < bnd[0]) bnd[0] = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > bnd[1]) bnd[1] = bucket[i].GetPhase(j);
            mean+=bucket[i].GetPhase(j);
        }
#else
        if (end-start<CRITPARALLELSIZE){
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPhase(j) < bnd[0]) bnd[0] = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > bnd[1]) bnd[1] = bucket[i].GetPhase(j);
            mean+=bucket[i].GetPhase(j);
        }
        }
        else {
        int nthreads;
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
        Double_t *mina=new Double_t[nthreads];
        Double_t *maxa=new Double_t[nthreads];
        for (i = 0; i < nthreads; i++)mina[i]=maxa[i]=bnd[0];
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for schedule(dynamic) nowait
        for (i = start + 1; i < end; i++)
        {
            if (bucket[i].GetPhase(j) < mina[omp_get_thread_num()]) mina[omp_get_thread_num()] = bucket[i].GetPhase(j);
            if (bucket[i].GetPhase(j) > maxa[omp_get_thread_num()]) maxa[omp_get_thread_num()] = bucket[i].GetPhase(j);
        }
    #pragma omp for reduction(+:mean)
        for (i = start+1; i < end; i++) mean+=bucket[i].GetPhase(j);
    }
        for (i = 0; i < nthreads; i++)
        {
            if (mina[i] < bnd[0]) bnd[0] = mina[i];
            if (maxa[i] > bnd[1]) bnd[1] = maxa[i];
        }
        delete[] mina;
        delete[] maxa;
        }
#endif
        mean/=(Double_t)(end-start);
        return mean;
    }
    //@}
    /// \name Find the dispersion in a dimension (biased variance using 1/N as opposed to 1/(N-1) so that if N=2, doesn't crash)
    //@{
    inline Double_t KDTree::DispersionPos(int j, Int_t start, Int_t end, Double_t mean)
    {
        Double_t disp=0;
        Int_t i;
#ifdef USEOPENMP
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for reduction(+:disp)
#endif
        for (i = start; i < end; i++)
            disp+=(bucket[i].GetPosition(j)-mean)*(bucket[i].GetPosition(j)-mean);
#ifdef USEOPENMP
    }
#endif
        disp/=(Double_t)(end-start);
        return disp;
    }
    inline Double_t KDTree::DispersionVel(int j, Int_t start, Int_t end, Double_t mean)
    {
        Double_t disp=0;
        Int_t i;
#ifdef USEOPENMP
    #pragma omp parallel default(shared)     \
    private(i)
    {
    #pragma omp for reduction(+:disp)
#endif
        for (i = start ; i < end; i++)
            disp+=(bucket[i].GetVelocity(j)-mean)*(bucket[i].GetVelocity(j)-mean);
#ifdef USEOPENMP
    }
#endif
        disp/=(Double_t)(end-start);
        return disp;
    }
    inline Double_t KDTree::DispersionPhs(int j, Int_t start, Int_t end, Double_t mean)
    {
        Double_t disp=0;
        Int_t i;
#ifdef USEOPENMP
    #pragma omp parallel default(shared) \
    private(i)
    {
    #pragma omp for reduction(+:disp)
#endif
        for (i = start; i < end; i++)
            disp+=(bucket[i].GetPhase(j)-mean)*(bucket[i].GetPhase(j)-mean);
#ifdef USEOPENMP
    }
#endif
        disp/=(Double_t)(end-start);
        return disp;
    }
    //@}

/*
    //code that applies a correction to the boundary of a node, ibnd is initial boundary range estimate,
    //xbnd is the current parent nodes estimate
    //at the moment the code is not setup to correct for underestimation in outer and inner parts
    //(I'm not exactly sure how this corrects this, must check Enbid paper)
    inline void BoundaryCor(int j, Int_t count, Int_t dim, Int_t numparts, Double_t *ibnd, Double_t *xbnd){
        //factors in the correction (taken from Enbind where don't assume cubic cells)
        Double_t fac1=10.0, fac2=0.2, fac3=2.0;
        Double_t temp1,temp2,temp3;
        temp1=fac1*pow((Double_t)numparts,1.0/(Double_t)dim);
        temp2=fac2/temp1;
        if((ibnd[1]-ibnd[0])/(xbnd[1]-xbnd[0])> 1.0/(Double_t)numparts)
        {
            temp3=(ibnd[1]-ibnd[0])/(count-1);
            if(((xbnd[1]-ibnd[1])>(temp2*temp3))&&((ibnd[0]-xbnd[0])>(temp2*temp3)))
            {
                xbnd[1]=ibnd[1]+fac3*temp3;
                xbnd[0]=ibnd[0]-fac3*temp3;
            }
            else
            {
                if((xbnd[1]-ibnd[1])>(temp1*temp3)) xbnd[1]=ibnd[1]+fac3*temp3;
                if((ibnd[0]-xbnd[0])>(temp1*temp3)) xbnd[0]=ibnd[0]-fac3*temp3;
            }
        }
    }
*/
     /// \name Calculate the entropy in a given dimension. This can be used as a node splitting criterion
     /// This calculates Shannon Entropy, where the region is split into nbins=pow(N,1/3) (N is number of particles) where minimum nbins=1,
     /// and can be used instead of most spread dimension
     //@{
    inline Double_t KDTree::EntropyPos(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *nientropy)
    {
        Int_t ibin, i;
        Double_t mtot=0.,entropy=0.;
        Double_t dx=(up-low)/nbins;
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=(Int_t)((bucket[i].GetPosition(j)-low)/dx);
            nientropy[ibin]+=bucket[i].GetMass();
        }
        mtot=1.0/mtot;
        for (i=0;i<nbins;i++) {
            if (nientropy[i]>0) {
                Double_t temp=nientropy[i]*mtot;
                entropy-=temp*log10(temp);
            }
        }
        return entropy/log10((Double_t)nbins);
    }
    inline Double_t KDTree::EntropyVel(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *nientropy)
    {
        Int_t ibin,i;
        Double_t mtot=0.,entropy=0.;
        Double_t dx=(up-low)/nbins;
        for (i=0;i<nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=(Int_t)((bucket[i].GetVelocity(j)-low)/dx);
            nientropy[ibin]+=bucket[i].GetMass();
        }
        mtot=1.0/mtot;
        for (i=0;i<nbins;i++) {
            if (nientropy[i]>0) {
                Double_t temp=nientropy[i]*mtot;
                entropy-=temp*log10(temp);
            }
        }
        return entropy/log10((Double_t)nbins);
    }
    inline Double_t KDTree::EntropyPhs(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *nientropy)
    {
        Int_t ibin,i;
        Double_t mtot=0.,entropy=0.;
        Double_t dx=(up-low)/nbins;
        for (i=0;i<=nbins;i++) nientropy[i]=0.;
        for (i=start;i<end;i++){
            mtot+=bucket[i].GetMass();
            ibin=(Int_t)((bucket[i].GetPhase(j)-low)/dx);
            nientropy[ibin]+=bucket[i].GetMass();
        }
        mtot=1.0/mtot;
        for (Int_t i=0;i<nbins;i++) {
            if (nientropy[i]>0) {
                Double_t temp=nientropy[i]*mtot;
                entropy-=temp*log10(temp);
            }
        }
        return entropy/log10(nbins);
    }
    //@}

    /// \name Determine the median coordinates in some space
    //@{
    inline Double_t KDTree::MedianPos(int d, Int_t k, Int_t start, Int_t end, bool balanced)
    {
        Int_t left = start;
        Int_t right = end-1;
        Int_t i, j;
        Double_t x;
        Particle w;

        //produced a balanced tree
        if (balanced){
        while (left < right)
        {
            x = bucket[k].GetPosition(d);
            w = bucket[right];
            bucket[right] = bucket[k];
            bucket[k] = w;
            i = left-1;
            j = right;
            while (1) {
                while (i < j) if (bucket[++i].GetPosition(d) >= x) break;
                while (i < j) if (bucket[--j].GetPosition(d) <= x) break;
                w = bucket[i];
                bucket[i] = bucket[j];
                bucket[j] = w;
                if (j <= i) break;
            }
            bucket[j] = bucket[i];
            bucket[i] = bucket[right];
            bucket[right] = w;
            if (i >= k) right = i - 1;
            if (i <= k) left = i + 1;
        }
        return bucket[k].GetPosition(d);
        }
        //much quicker but does not guarantee a balanced tree
        else
        {
            printf("Note yet implemented\n");
            exit(9);
        }
    }
    inline Double_t KDTree::MedianVel(int d, Int_t k, Int_t start, Int_t end, bool balanced)
    {
        Int_t left = start;
        Int_t right = end - 1;
        Int_t i, j;
        Double_t x;
        Particle w;

        if (balanced){
        while (left < right)
        {
            x = bucket[k].GetVelocity(d);
            w = bucket[right];
            bucket[right] = bucket[k];
            bucket[k] = w;
            i = left-1;
            j = right;
            while (1) {
                while (i < j) if (bucket[++i].GetVelocity(d) >= x) break;
                while (i < j) if (bucket[--j].GetVelocity(d) <= x) break;
                w = bucket[i];
                bucket[i] = bucket[j];
                bucket[j] = w;
                if (j <= i) break;
            }
            bucket[j] = bucket[i];
            bucket[i] = bucket[right];
            bucket[right] = w;
            if (i >= k) right = i - 1;
            if (i <= k) left = i + 1;
        }
        return bucket[k].GetVelocity(d);
        }
        //much quicker but does not guarantee a balanced tree
        else
        {
            printf("Note yet implemented\n");
            exit(9);
        }
    }
    inline Double_t KDTree::MedianPhs(int d, Int_t k, Int_t start, Int_t end, bool balanced)
    {
        Int_t left = start;
        Int_t right = end - 1;
        Int_t i, j;
        Double_t x;
        Particle w;

        if (balanced){
        while (left < right)
        {
            //if (d<3) x = bucket[k].GetPosition(d);
            //else x = bucket[k].GetVelocity(d);
            x=bucket[k].GetPhase(d);
            w = bucket[right];
            bucket[right] = bucket[k];
            bucket[k] = w;
            i = left-1;
            j = right;
            while (1) {
                while (i < j) {if (bucket[++i].GetPhase(d) >= x) break;}
                while (i < j) {if (bucket[--j].GetPhase(d) <= x) break;}
                w = bucket[i];
                bucket[i] = bucket[j];
                bucket[j] = w;
                if (j <= i) break;
            }
            bucket[j] = bucket[i];
            bucket[i] = bucket[right];
            bucket[right] = w;
            if (i >= k) right = i - 1;
            if (i <= k) left = i + 1;
        }

        return bucket[k].GetPhase(d);
        }
        //much quicker but does not guarantee a balanced tree
        else
        {
            printf("Note yet implemented\n");
            exit(9);
        }
    }
    //@}
    //-- End of inline functions

    //-- Private functions used to build tree

    /// Recursively build the nodes of the tree.  This works by first finding the dimension under
    /// which the data has the most spread, and then splitting the data about the median
    /// in that dimension.  BuildNodes() is then called on each half. Once the size of the data is
    /// small enough, a leaf node is formed.
    Node *KDTree::BuildNodes(Int_t start, Int_t end)
    {
        Double_t bnd[6][2];
        Int_t size = end - start;
        if (size <= b)
        {
            numleafnodes++;numnodes++;
            for (int j=0;j<ND;j++) (this->*bmfunc)(j, start, end, bnd[j]);
            return new LeafNode(numnodes-1,start, end,  bnd, ND);
        }
        else
        {
            numnodes++;
            int splitdim=0,j, enflag;
            Int_t k = start + (size - 1) / 2;
            Double_t maxspread, minentropy, maxsig,splitvalue;
            Double_t nbins;
            //if using shannon entropy criterion
            if(splittingcriterion==1) if(end-start>8) nbins=ceil(pow((end-start),1./3.));else nbins=2;
#ifdef USEOMP
            //if openmp and region is above the critical size then run otherwise, set omp num threads=1;
            int nthreads;
#pragma omp parallel
    {
            if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
    }
            if (size<CRITPARALLELSIZE)nthreads=1;
            else nthreads=6;
#pragma omp master
    {
            omp_set_num_threads(nthreads);
    }
    #pragma omp parallel default(shared) \
    private(j)
    {
    #pragma omp for schedule(dynamic,1) nowait
#endif
            for (j = 0; j < ND; j++)
            {
                if(splittingcriterion==1) {
                    spreada[j] = (this->*spreadfunc)(j, start, end, bnd[j])+1e-32;//addition incase lattice and no spread
                    Double_t low, up;
                    low=bnd[j][0]-2.0*(spreada[j])/(Double_t)(end-start);
                    up=bnd[j][1]+2.0*(spreada[j])/(Double_t)(end-start);
                    entropya[j] = (this->*entropyfunc)(j, start, end, low, up, nbins, nientropy[j]);
                }
                else if (splittingcriterion==2) {
                    meana[j] = (this->*bmfunc)(j, start, end, bnd[j]);
                    vara[j] = (this->*dispfunc)(j, start, end, meana[j]);
                }
                else {
                    spreada[j] = (this->*spreadfunc)(j, start, end, bnd[j]);
                }
            }
#ifdef USEOMP
    }
#endif


            splitdim=0; maxspread=spreada[0]; minentropy=entropya[0];maxsig=vara[0];
            //splitdim=0; maxspread=0.0; minentropy=1.0;enflag=0;
            //for since entropy can only be used in cases where the subspace is not sparse or does not have lattice structure must check
            //the question is how? At the moment, I do not check for this, though the idea would be only check dimensions that meet the criteria
            //and if non of them meet it, then enflag still zero and perhaps, use most spread dimension
            for (j = 1; j < ND; j++)
            {
                if(splittingcriterion==1) {
                    if (minentropy>=entropya[j]) {
                        minentropy = entropya[j];
                        splitdim = j;
                    }
                }
                else if (splittingcriterion==2) {
                    if (vara[j] > maxsig) {
                        maxsig = vara[j];
                        splitdim = j;
                    }
                }
                else {
                    if (spreada[j] > maxspread) {
                        maxspread = spreada[j];
                        splitdim = j;
                    }
                }
            }

            splitvalue= (this->*medianfunc)(splitdim, k, start, end,true);

            return new SplitNode(numnodes-1, splitdim, splitvalue, size, bnd, start, end, ND, BuildNodes(start, k+1),BuildNodes(k+1, end));
        }
    }

    ///scales the space and calculates the corrected volumes
    ///note here this is not mass weighted which may lead to issues later on.
    void KDTree::ScaleSpace(){

        if (treetype==TPHYS||treetype==TPROJ)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xmean[j]+=bucket[i].GetPosition(j);
        else if (treetype==TVEL)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xmean[j]+=bucket[i].GetVelocity(j);
        else if (treetype==TPHS||treetype==TMETRIC)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xmean[j]+=bucket[i].GetPhase(j);

        for(int j=0;j<ND;j++) xmean[j]/=(Double_t)numparts;

        if (treetype==TPHYS||treetype==TPROJ)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xvar[j]+=(bucket[i].GetPosition(j)-xmean[j])*(bucket[i].GetPosition(j)-xmean[j]);
        else if (treetype==TVEL)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xvar[j]+=(bucket[i].GetVelocity(j)-xmean[j])*(bucket[i].GetVelocity(j)-xmean[j]);
        else if (treetype==TPHS||treetype==TMETRIC)
            for(Int_t i=0; i<numparts; i++)
                for(int j=0;j<ND;j++) xvar[j]+=(bucket[i].GetPhase(j)-xmean[j])*(bucket[i].GetPhase(j)-xmean[j]);

        for(int j=0;j<ND;j++){xvar[j]=sqrt(xvar[j]/(Double_t)numparts);ixvar[j]=1./xvar[j];}
        if (treetype==TPHYS||treetype==TPROJ)
            for (Int_t i=0;i<numparts;i++)
    	        for (int j=0;j<ND;j++) bucket[i].SetPosition(j,bucket[i].GetPosition(j)*ixvar[j]);
        else if (treetype==TVEL)
            for (Int_t i=0;i<numparts;i++)
    	        for (int j=0;j<ND;j++) bucket[i].SetVelocity(j,bucket[i].GetVelocity(j)*ixvar[j]);
        else if (treetype==TPHS||treetype==TMETRIC)
            for (Int_t i=0;i<numparts;i++)
    	        for (int j=0;j<ND;j++) bucket[i].SetPhase(j,bucket[i].GetPhase(j)*ixvar[j]);
    }

    int KDTree::TreeTypeCheck(){
        if (treetype<TPHYS||treetype>TMETRIC) {
            printf("Error in type of tree specified\n");
            printf("Returning null pointer as root node and leaving particle list unchanged.\n");
            root=NULL; return 0;
        }
        else if (treetype==TMETRIC&&metric==NULL) {
            printf("For treetype=%d must pass metric\n",TMETRIC);
            printf("Returning null pointer as root node and leaving particle list unchanged.\n");
            root=NULL; return 0;
        }
        else {
        if (treetype==TPHYS)
        {
            bmfunc=&NBody::KDTree::BoundaryandMeanPos;
            dispfunc=&NBody::KDTree::DispersionPos;
            spreadfunc=&NBody::KDTree::SpreadestPos;
            entropyfunc=&NBody::KDTree::EntropyPos;
            medianfunc=&NBody::KDTree::MedianPos;
        }
        else if (treetype==TVEL) {
            bmfunc=&NBody::KDTree::BoundaryandMeanVel;
            dispfunc=&NBody::KDTree::DispersionVel;
            spreadfunc=&NBody::KDTree::SpreadestVel;
            entropyfunc=&NBody::KDTree::EntropyVel;
            medianfunc=&NBody::KDTree::MedianVel;
        }
        else if (treetype==TPHS) {
            bmfunc=&NBody::KDTree::BoundaryandMeanPhs;
            dispfunc=&NBody::KDTree::DispersionPhs;
            spreadfunc=&NBody::KDTree::SpreadestPhs;
            entropyfunc=&NBody::KDTree::EntropyPhs;
            medianfunc=&NBody::KDTree::MedianPhs;
        }
        else if (treetype==TPROJ) {
            bmfunc=&NBody::KDTree::BoundaryandMeanPos;
            dispfunc=&NBody::KDTree::DispersionPos;
            spreadfunc=&NBody::KDTree::SpreadestPos;
            entropyfunc=&NBody::KDTree::EntropyPos;
            medianfunc=&NBody::KDTree::MedianPos;
        }
        return 1;
        }
    }
    ///Calculate kernel quantities
    void KDTree::KernelConstruction(){
        if (kernres<100) {
            printf("Error in kernel resolution, <100, setting to 100\n");
            kernres=100;
        }
        Kernel=new Double_t[kernres];
        derKernel=new Double_t[kernres];

        //determine dimension for kernel normalization
        if (treetype==TPHYS) ND=3;
        else if (treetype==TVEL) ND=3;
    	else if (treetype==TPHS) ND=6;
        else if (treetype==TPROJ) ND=2;
        else if (treetype==TMETRIC) ND=6;
        kernnorm=1.0/ND*pow(M_1_PI,ND/2.)*gsl_sf_gamma(ND/2.+1.);

        //smoothing kernel type
        if (kernfunctype==KSPH) {
            kernfunc.function=WSPH;
            kernnorm*=ND*(ND+1.)*(ND+2.)*(ND+3.)/(6*(pow(2.,ND+1)-1.));
        }
        else if (kernfunctype==KGAUSS) {
            kernfunc.function=WGauss;
            kernnorm*=pow(0.5*M_1_PI,ND/2.);
        }
        else if (kernfunctype==KEPAN) {
            kernfunc.function=WEpan;
            kernnorm*=ND*(ND+2.)*pow(0.5,ND+1.);
        }
        else if (kernfunctype==KTH) kernfunc.function=WTH;
        else {
            printf("Error in kernel choice, using default SPH B-spline kernel\n");
            kernfunc.function=WSPH;;
            kernnorm*=ND*(ND+1.)*(ND+2.)*(ND+3.)/(6*(pow(2.,ND+1)-1.));
        }
        Double_t delta=2.0/(Double_t)(kernres-1);
        for (int i=0;i<kernres;i++) {
            Kernel[i]=kernnorm*kernfunc.function(i*delta,1.0);
        }
    }

    //-- End of private functions used to build the tree

    //-- Public constructors

    KDTree::KDTree(Particle *p, Int_t nparts, Int_t bucket_size, int ttype, int smfunctype, int smres, int criterion, int aniso, int scale, Double_t *Period, Double_t **m)
    {
        numparts = nparts;
        numleafnodes=numnodes=0;
        bucket = p;
        b = bucket_size;
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (Period!=NULL)
        {
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=Period[k];
        }
        else period=NULL;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            root=BuildNodes(0,numparts);
            //else if (treetype==TMETRIC) root = BuildNodesDim(0, numparts,metric);
            if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
    }

    KDTree::KDTree(System &s, Int_t bucket_size, int ttype, int smfunctype, int smres, int criterion, int aniso, int scale, Double_t **m)
    {
//        KDTree(s.Parts(),s.GetNumParts(),bucket_size,ttype,smfunctype,smres,ecalc,aniso,scale,s.GetPeriod().GetCoord(),m);

        numparts = s.GetNumParts();
        numleafnodes=numnodes=0;
        bucket = s.Parts();
        b = bucket_size;
        treetype = ttype;
        kernfunctype = smfunctype;
        kernres = smres;
        splittingcriterion = criterion;
        anisotropic=aniso;
        scalespace = scale;
        metric = m;
        if (s.GetPeriod()[0]>0&&s.GetPeriod()[1]>0&&s.GetPeriod()[2]>0){
            period=new Double_t[3];
            for (int k=0;k<3;k++) period[k]=s.GetPeriod()[k];
        }
        else period=NULL;
        if (TreeTypeCheck()) {
            KernelConstruction();
            for (Int_t i = 0; i < numparts; i++) bucket[i].SetID(i);
            vol=1.0;ivol=1.0;
            for (int j=0;j<ND;j++) {xvar[j]=1.0;ixvar[j]=1.0;}
            if (scalespace) ScaleSpace();
            for (int j=0;j<ND;j++) {vol*=xvar[j];ivol*=ixvar[j];}
            if (splittingcriterion==1) for (int j=0;j<ND;j++) nientropy[j]=new Double_t[numparts];
            root=BuildNodes(0,numparts);
            if (splittingcriterion==1) for (int j=0;j<ND;j++) delete[] nientropy[j];
        }
    }
    KDTree::~KDTree()
    {
	    if (root!=NULL) {
            delete root;
            delete[] Kernel;
            delete[] derKernel;
            if (period!=NULL) delete[] period;
            qsort(bucket, numparts, sizeof(Particle), IDCompare);
            if (scalespace) {
            for (Int_t i=0;i<numparts;i++)
                for (int j=0;j<3;j++) {
                    bucket[i].SetPosition(j,bucket[i].GetPosition(j)*xvar[j]);
                    bucket[i].SetVelocity(j,bucket[i].GetVelocity(j)*xvar[j+3]);
                }
            }
        }
    }
}
