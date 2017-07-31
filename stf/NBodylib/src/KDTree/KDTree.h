/*! \file KDTree.h
 *  \brief header file for the \ref NBody::KDTree Class

    \class NBody::KDTree
    \brief A binary kd tree class, based on the implementation described in Friedman, Bentley, and Finkel (1977) and Arya and Mount (1993).\n

    A kd tree, based on the implementation described in Friedman, Bentley, and Finkel (1977) and Arya and Mount (1993).\n
    Now tree can span position, velocity and full phase space (also can calculate distance using a metric). \n
    This tree uses the \ref NBody::Node class (see \ref KDNode.h), \ref NBody::Particle class (see \ref Particle.h),
    and can also use the \ref NBody::System class (see \ref System.h)

    The class will create a kd tree from an \ref NBody::System or \ref NBody::Particle array.
    \b Warning: this will permanently change the particle array by re-ordering the particles so that the tree
    is balanced locally. Don't worry, when the tree is deleted, the particles will be put back the way they were prior
    to the construction the tree using particle ids which are rewritten to store the original particle order upon creation of the tree

    Note that to make density/sph calculations parallel using MPI, have to split tree up into blocks
    once main tree constructed, then give thread a subregion, make new tree and search that tree minimizing
    communication. This would be quite complicated to implement but still feasible.
    Another possibility is using something like what Gadget does with the Peno-Hilbert space filling curve but that seems excessive.

    Or for openmp a simple way of parallelizing the calculation is to use a slightly biased method.
    That is instead of using true sph kernel calculation to calculate a quantity for particle \f$ i \f$
    \f[ (0.5W(R_{ij},h_i)+0.5W(R_{ij},h_j)) \f]
    use
    \f[ W(R_{ij},h_i) \f]
    This removes the determination of \f$ h_j \f$.

    \todo NOTE: STILL MUST ALTER SMOOTHQUANTITIES SO THAT PROJECTED VALUES IMPLEMENTED IN A REASONABLE FASHION
*/

#ifndef KDTREE_H
#define KDTREE_H

#include <NBody.h>
#include <PriorityQueue.h>
#include <KDNode.h>
#include <NBodyMath.h>
#include <SmoothingKernels.h>
#include <FOFFunc.h>

#include <iomanip>
#include <cstdio>
#include <iostream>

#ifdef USEOPENMP
#include <omp.h>
#define CRITPARALLELSIZE 1000000
#endif

#ifdef USEMPI
#include <mpi.h>
#define CRITPARALLELSIZE 1000000
#endif

namespace NBody
{


    class KDTree
    {

        public:

        /// \name Public variables that specify treetype
        /// 0 is physical tree, 1 is projected physical, 2 is velocity, 3 is full phase-space, 4 is a sub-space of the full space specified by a metric
        /// ie: TPHYS=0,TPROJ=1,TVEL=2,TPHS=3,TMETRIC=4
        //@{
        const static int TPHYS=0,TPROJ=1,TVEL=2,TPHS=3,TMETRIC=4;
        //@}

        /// \name Public variables that specify kernaltype
        /// 0 is standard sph B-spline, 1 is gaussian, 2 is Epanechnikov, 3 is tophat;
        /// ie: KSPH=0,KEPAN=2,KGAUSS=3,KTH=3;
        //@{
        const static int KSPH=0,KGAUSS=1,KEPAN=2,KTH=3;
        //@}

        protected:

        private:

        /// number of parts in tree
        Int_t numparts;
        ///pointer to particles
        Particle *bucket;
        ///pointer to root node pointer.
        Node *root;
        ///if system is periodic, period in each direction
        Double_t *period;
        ///number of nodes and leafnodes
        Int_t numnodes,numleafnodes;
        ///bucket or leaf node size
        Int_t b;
        ///max number of dimensions of tree
        const static int MAXND=6;

        ///for an arbitrary tree spanning some space one would have offsets in the dimensional space to use
        ///something like \code int startdim,enddim; \endcode \n
        ///but here more appropriate to specify the type of tree that dictates the dimension
        int treetype, ND;

        ///if scale the space, then calculate the mean and variance in a given dimension and the change in the volume and inverse due to scaling
        int scalespace;
        Double_t xmean[MAXND],xvar[MAXND],ixvar[MAXND],vol,ivol;

        ///0 if using most spread dimension as criterion, 1 if use entropy along with entropy array, 2 if using largest dispersion
        int splittingcriterion;
        Double_t *nientropy[MAXND];

        ///kernel construction
        ///resolution in kernel array and type
        int kernres,kernfunctype;
        ///contains arrays used to calculate smoothing
        Double_t *Kernel,*derKernel;
        ///function pointer to smoothing kernel
        smoothfunc kernfunc;
        ///normalization of kernel which depends on dimension of tree
        Double_t kernnorm;

        ///flag if in treetype 2 and wish to account for anisotropic distribution about particle
        ///-1 don't use local phase-space metric, 0 use phase-space metric but assume isotropy, 1 anisotropic
        int anisotropic;

        ///array for user specified metric or metric calculated based on distribution of particles in the tree space
        ///this array is necessary for calculating the phase-space density corectly.
        Double_t **metric;

        /// \name Private function pointers used in building tree
        //@{
        Double_t(NBody::KDTree::*bmfunc)(int , Int_t , Int_t , Double_t *);
        Double_t(NBody::KDTree::*dispfunc)(int , Int_t, Int_t, Double_t);
        Double_t(NBody::KDTree::*spreadfunc)(int , Int_t , Int_t , Double_t *);
        Double_t(NBody::KDTree::*entropyfunc)(int , Int_t , Int_t , Double_t , Double_t, Double_t, Double_t *);
        Double_t(NBody::KDTree::*medianfunc)(int , Int_t , Int_t, Int_t, bool);
        //@}

        /// \name Private arrays used to build tree
        //@{
        Double_t spreada[MAXND];
        Double_t meana[MAXND];
        Double_t vara[MAXND];
        Double_t entropya[MAXND];
        //@}

        /// \name Tree construction methods
        /// Private methods used in constructing the tree
        //@{
        ///uses prviate function pointers to recursive build the tree
        Node* BuildNodes(Int_t start, Int_t end);
        ///scales the space if necessary by the variance in each dimension
        void ScaleSpace();
        ///checks to see if tree is of proper type
        int TreeTypeCheck();
        ///build the table of kernel values
        void KernelConstruction();
        //@}

        public :

        /// \name Constructors/Destructors
        //@{
        ///Creates tree from an NBody::Particle array
        KDTree(Particle *p, Int_t numparts, Int_t bucket_size = 16, int TreeType=TPHYS, int KernType=KEPAN, int KernRes=1000, int SplittingCriterion=0, int Aniso=0, int ScaleSpace=0, Double_t *Period=NULL, Double_t **metric=NULL);
        ///Creates tree from NBody::System
        KDTree(System &s, Int_t bucket_size = 16, int TreeType=TPHYS, int KernType=KEPAN, int KernRes=1000, int SplittingCriterion=0, int Aniso=0, int ScaleSpace=0, Double_t **metric=NULL);
        ///resets particle order
        ~KDTree();
        //@}

        /// \name Simple Get functions
        //@{
        Int_t GetNumNodes(){return numnodes;}
        Int_t GetNumLeafNodes(){return numleafnodes;}
        Int_t GetBucketSize(){return b;}
        Int_t GetTreeType(){return treetype;}
        Int_t GetKernType(){return kernfunctype;}
        Double_t GetKernNorm(){return kernnorm;}
        Node * GetRoot(){return root;}
        Double_t GetPeriod(int j){return period[j];}
        //@}

        /// \name Find Nearest Neighbour functions
        /// using tree, for each particle find nearest.
        /// FindNearest is using treetype, rest are explicit search using position, velocity, phase regardless of treetype.
        /// The nn array stores the indices of the particles, not their ids. \n
        /// These have multiple interfaces: \n
        /// using tree, find nearest particles for entire system, or find nearest to a particle bucket[tt] or position x, v and (x,v) \n
        /// Implementation in \ref KDFindNearest.cxx
        //@{
        void FindNearest(Int_t **nn, Double_t **dist2, Int_t Nsearch=64);
        void FindNearestPos(Int_t **nn, Double_t **dist2, Int_t Nsearch=64);
        void FindNearestVel(Int_t **nn, Double_t **dist2, Int_t Nsearch=64);
        void FindNearestPhase(Int_t **nn, Double_t **dist2, Int_t Nsearch=64);
        //using tree nearest particles that meet a criterion
        void FindNearestCriterion(FOFcompfunc cmp, Double_t *params,Int_t **nn, Double_t **dist2, Int_t Nsearch=64);

        void FindNearest(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestPos(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestVel(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestPhase(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);

        void FindNearest(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestPos(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestVel(Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestPhase(Double_t *x, Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);

        void FindNearestPos(Coordinate x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestVel(Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        void FindNearestPhase(Coordinate x, Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64);

        //using tree nearest particles that meet a criterion relative to particle bucket[tt]
        void FindNearestCriterion(Int_t tt, FOFcompfunc cmp, Double_t *params,Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        //using tree nearest particles that meet a criterion relative to particle passed
        void FindNearestCriterion(Particle p, FOFcompfunc cmp, Double_t *params,Int_t *nn, Double_t *dist2, Int_t Nsearch=64);
        //@}

        /// \name Search for all particles within a given distance
        /// using tree find all particles within a distance fdist2 to particle bucket[tt], or position
        /// the return array nn and dist2 need to be size of bucket and all particles that are within a distance
        /// dist2 are marked with imark in nn array which is in particle id order. \n
        /// Implementation in \ref KDFindNearest.cxx
        //@{

        void SearchBall(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallPos(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallVel(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallPhase(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);

        //same but search centered on position
        void SearchBall(Double_t *x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallPos(Double_t *x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallVel(Double_t *v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallPhase(Double_t *x, Double_t *v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);

        void SearchBall(Coordinate x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallPos(Coordinate x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallVel(Coordinate v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchBallPhase(Coordinate x, Coordinate v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2);

        //using tree find all particles that meet a criterion relative to particle bucket[tt]
        void SearchCriterion(Int_t tt, FOFcompfunc p, Double_t *params, Int_t imark, Int_t *nn, Double_t *dist2);
        void SearchCriterion(Int_t tt, FOFcompfunc p, Double_t *params, Int_t imark, Int_t *nn);
        //similar to above but also returns array of tagged and numtagged, where tagged is an array of index of tagged items (unsorted)
        //containing numtagged particles (note array must be of size of particle array)
        //Int_t SearchCriterion(Int_t tt, FOFcompfunc p, Double_t *params, Int_t imark, Int_t *nn, Double_t *dist2, Int_t *tagged);
        //here SearchBallPos but specifically RETURN array of size Ntagged of index of particles tagged.
        Int_t SearchBallPosTagged(Int_t tt, Double_t fdist2, Int_t *tagged);
        Int_t SearchBallPosTagged(Coordinate x, Double_t fdist2, Int_t *tagged);
        Int_t SearchBallPosTagged(Double_t *x, Double_t fdist2, Int_t *tagged);
        //return number of tagged particles meeting a criterion
        Int_t SearchCriterionTagged(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t *tagged);
        Int_t SearchCriterionTagged(Particle &p, FOFcompfunc cmp, Double_t *params, Int_t *tagged);
        //@}


        /// \name Find leaf nodes
        /// return the leaf node containing target particle, or position, or within boundaries of search. \n
        /// Implementation in \ref KDFindNearest.cxx
        //@{
        Node* FindLeafNode(Int_t tt);
        Node* FindLeafNode(Double_t *x);
        Node* FindLeafNode(Double_t search[6][2]);
        //@}

        /*! \name Calculates SPH quantities
            Using tree to calculate SPH quantities using kernel smoothing (ie: SPH).
            To calculate smooth quantities other than density, one needs the density. \n
            There are multiple interfaces: \n
            \li One can calculate a quantity using standard kernel technique for entire system \n
            \li Or calculate it for a single target particle or at a specific position.
            note that for these calculations, do not take 0.5 *smoothing kernel and add
            the contribution of W(Rij,hjmax) where hjmax is the jths particle most distnace neighbour
            instead simply smooth over near neighbours and take 1.0 *smoothing kernel \n

            \b NOTE: The arrays that are returned are ordered according to particle ids that is value[i] is the value of particle with ID=i.\n
            Extra routines can be declared here and added to \ref KDCalcSmoothQuantities.cxx
        */
        //@{

        /// Calculate the density based on the tree constructed.
        /// That is treetype=0 then physical density, 1 then velocity and 2 is phase, 3 is projected physical
        /// 4 is the phase density specified using a metric.
        void CalcDensity(Int_t Nsmooth=64);
        /// find nearest Nsearch physical neighbours (treetype must be 0 or 3) and get velocity density
        /// by smoothing over the nearest Nsmooth velocity neighbours from the Nsearch set.
        void CalcVelDensity(Int_t Nsmooth=64, Int_t Nsearch=64);
        /// similar to above but return the physical density weighted value of the velocity density
        /// that is fvel weighted by the overlap of the particles in physical space given by W(rij)/rho(i)*M(i)
        Double_t *CalcVelDensityWithPhysDensity(Int_t Nsmooth=64, Int_t Nsearch=64,int densityset=1);

        /// Calculates the smoothed mean velocity
        Coordinate *CalcSmoothVel(Int_t Nsmooth=64, int densityset=1);
        /// Calculates the velocity dispersion tensor, if meanvelset!=1 then must pass appropriate smvel array
        Matrix *CalcSmoothVelDisp(Coordinate *smvel, Int_t Nsmooth=64, int densityset=1, int meanvelset=1);
        /// Calculates the smoothed velocity skewness, checks if various calculations done
        Coordinate *CalcSmoothVelSkew(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth=64, int densityset=1, int meanvelset=1, int veldispset=1);
        /// Calculates the smoothed velocity kurtosis
        Coordinate *CalcSmoothVelKurtosis(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth=64, int densityset=1, int meanvelset=1, int veldispset=1);

        Double_t CalcDensityParticle(Int_t target, Int_t Nsmooth=64);
        Double_t CalcVelDensityParticle(Int_t target, Int_t Nsmooth=64, Int_t Nsearch=64, int iflag=0, PriorityQueue *pq=NULL, PriorityQueue *pq2=NULL, Int_t *nnIDs=NULL, Double_t *vdist=NULL);
        Double_t CalcVelDensityWithPhysDensityParticle(Int_t target, Int_t Nsmooth=64, Int_t Nsearch=64,int densityset=1);
        Coordinate CalcSmoothVelParticle(Int_t target, Int_t Nsmooth=64, int densityset=1);
        Matrix CalcSmoothVelDispParticle(Int_t target, Coordinate smvel, Int_t Nsmooth=64, int densityset=1);

        //Same as above but for a particular coordinate with two interfaces (Double_t pointers and Coordinate)
        //note for density if phase-space density then must specify velocity position.
        Double_t CalcDensityPosition(Double_t *x, Int_t Nsmooth=64, Double_t *v=NULL);
        Double_t CalcVelDensityPosition(Double_t *x, Double_t *v, Int_t Nsmooth=64, Int_t Nsearch=64);
        //Double_t CalcVelDensityWithPhysDensityPosition(Double_t *x, Double_t *v, Int_t Nsmooth=64, Int_t Nsearch=64,int densityset=1);
        Coordinate CalcSmoothVelPosition(Double_t *x, Int_t Nsmooth=64, int densityset=1);
        Matrix CalcSmoothVelDispPosition(Double_t *x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1);
        Double_t CalcDensityPosition(Coordinate x, Int_t Nsmooth=64, Coordinate v=Coordinate(0.));
        Double_t CalcVelDensityPosition(Coordinate x, Int_t Nsmooth=64, Int_t Nsearch=64);
        //Double_t CalcVelDensityWithPhysDensityPosition(Coordinate x, Int_t Nsmooth=64, Int_t Nsearch=64,int densityset=1);
        Coordinate CalcSmoothVelPosition(Coordinate x, Int_t Nsmooth=64, int densityset=1);
        Matrix CalcSmoothVelDispPosition(Coordinate x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1);

        Double_t CalcDensityPosition(Real_t *x, Int_t Nsmooth=64, Real_t *v=NULL);
        Double_t CalcVelDensityPosition(Real_t *x, Real_t *v, Int_t Nsmooth=64, Int_t Nsearch=64);
        Coordinate CalcSmoothVelPosition(Real_t *x, Int_t Nsmooth=64, int densityset=1);
        Matrix CalcSmoothVelDispPosition(Real_t *x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1);

        // Calculates the local value of a quantity using Sph kernel. Weights and return values are ordered according to original IDs
        // prior to tree construction.
        Double_t *CalcSmoothLocalMean(Double_t *weight, Int_t Nsmooth=64, int densityset=1);
        Double_t *CalcSmoothLocalDisp(Double_t *weight, Double_t *localmean, Int_t Nsmooth=64, int densityset=1);

        // Same as above but for individual particle or position
        Double_t CalcSmoothLocalMeanParticle(Int_t target, Double_t *weight, Int_t Nsmooth=64, int densityset=1);
        Double_t CalcSmoothLocalDispParticle(Int_t target, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1);
        Double_t CalcSmoothLocalMeanPosition(Double_t *x, Double_t *weight, Int_t Nsmooth=64, int densityset=1);
        Double_t CalcSmoothLocalDispPosition(Double_t *x, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1);
        Double_t CalcSmoothLocalMeanPosition(Coordinate x, Double_t *weight, Int_t Nsmooth=64, int densityset=1);
        Double_t CalcSmoothLocalDispPosition(Coordinate x, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1);

        ///Calculate quantity having already found distances stored in a PriorityQueue along with weights
        Double_t CalcSmoothLocalValue(Int_t Nsmooth, PriorityQueue *pq, Double_t* weight);

        //@}

        /*! \name FOF search routines
            Using tree to quickly search all particles that meet the criteria of the FOF comparison fuction
            There are multiple interfaces: \n
            \li One can search entire system \n
            \li Or search system to see what particles meet the criteria relative to a given particle

            \b NOTE: The arrays that are returned are ordered according to particle ids that is fof[i] is the group value of particle with ID=i.\n
            Extra routines can be declared here and added to \ref KDFOF.cxx
        */
        //@{

        /// simple physical FOF search, velocity FOF, and 6D phase FOF.
        Int_t *FOF(Double_t fdist, Int_t &numgroup, Int_t minnum=8, int order=0, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL);
        /// this searches tree for particles that meed some criterion given by some comparison function
        /// NOTE: comparison function is used only for leaf nodes.
        /// NOTE: parameters for tree search and comparison are contained in the params variable
        /// The first 5 indices of the params variable indicating how tree is searched.
        /// The [0] address is used to store treetype. The second value and third value are the
        /// associated distance^2 measures for physical and velocity space, the next 3 spaces are
        /// left for possible changes.
        /// Parameters used in comparison function start at params[6] & the use depends on the comparison function
        /// note that numgroups is total number of groups found but in fof.grp array returned, group zero is not part of this list and
        /// corresponds to particles that were not part of any group
        /*
        Int_t *FOFCriterion(int p(Particle, Particle, Double_t *), Double_t *params, Int_t &numgroups, Int_t minnum=8, int order=0);
        Int_t FOFCriterionParticle(int p(Particle, Particle, Double_t *), Int_t *pfof, Int_t target, Int_t iGroup, Double_t *params);
        /// Similar to above but only search a specific particle list for each particle given by [numparts][numberNN] array.
        Int_t *FOFNNCriterion(int p(Particle, Particle, Double_t *), Double_t *params, Int_t numNN, Int_t **nnIDs, Int_t &numgroups, Int_t minnum=8);
        /// Similar to above but uses distances as well by passing distance array to a function that returns a double, this value is then used to alter parameters
        /// of comparison function where npc is number of parameters to be changed and npca is array containing indices of parameters to be altered.
        Int_t *FOFNNDistCriterion(int p(Particle, Particle, Double_t *), Double_t *params, Int_t numNN, Int_t **nnIDs, Double_t **dist2,
                                  Double_t disfunc(Int_t , Double_t *), Int_t npc, Int_t *npca, Int_t &numgroups, Int_t minnum=8);
                                  */
        Int_t *FOFCriterion(FOFcompfunc p, Double_t *params, Int_t &numgroups, Int_t minnum=8, int order=0, int ipcheckflag=0, FOFcheckfunc check=Pnocheck, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL);
        ///Link partiles based on criterion and also only allow specific particles to be used as a basis for links using FOFcheckfunc
        Int_t* FOFCriterionSetBasisForLinks(FOFcompfunc cmp, Double_t *params, Int_t &numgroup, Int_t minnum=8, int order=0, int ipcheckflag=0, FOFcheckfunc check=Pnocheck, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL);
        //Int_t FOFCriterionParticle(FOFcompfunc p, Int_t *pfof, Int_t target, Int_t iGroup, Double_t *params);
        Int_t FOFCriterionParticle(FOFcompfunc p, Int_t *pfof, Int_t target, Int_t iGroup, Double_t *params, Int_tree_t *pGroupHead, Int_tree_t *Fifo, Int_tree_t *pHead, Int_tree_t *pTail, Int_tree_t *pNext, Int_tree_t *pLen);
        /// Similar to above but only search a specific particle list for each particle given by [numparts][numberNN] array.
        Int_t *FOFNNCriterion(FOFcompfunc p, Double_t *params, Int_t numNN, Int_t **nnIDs, Int_t &numgroups, Int_t minnum=8);
        /// Similar to above but uses distances as well by passing distance array to a function that returns a double, this value is then used to alter parameters
        /// of comparison function where npc is number of parameters to be changed and npca is array containing indices of parameters to be altered.
        Int_t *FOFNNDistCriterion(FOFcompfunc p, Double_t *params, Int_t numNN, Int_t **nnIDs, Double_t **dist2,
                                  Double_t disfunc(Int_t , Double_t *), Int_t npc, Int_t *npca, Int_t &numgroups, Int_t minnum=8);
        //@}

        private:

        //-- private inline functions declarations

        /// \name Splitting criteria methods
        /// How to split the system
        //@{

        /// Find the dimension of which the data has the most spread
        /// for positions
        inline Double_t SpreadestPos(int j, Int_t start, Int_t end, Double_t *bnd);
        /// and velocities
        inline Double_t SpreadestVel(int j, Int_t start, Int_t end, Double_t *bnd);
        /// and phase
        inline Double_t SpreadestPhs(int j, Int_t start, Int_t end, Double_t *bnd);
        /// Find the boundary of the data and return mean
        /// for positions
        inline Double_t BoundaryandMeanPos(int j, Int_t start, Int_t end, Double_t *bnd);
        /// and velocities
        inline Double_t BoundaryandMeanVel(int j, Int_t start, Int_t end, Double_t *bnd);
        /// and phs
        inline Double_t BoundaryandMeanPhs(int j, Int_t start, Int_t end, Double_t *bnd);
        /// Find the dispersion in a dimension
        /// for positions
        inline Double_t DispersionPos(int j, Int_t start, Int_t end, Double_t mean);
        /// and velocities
        inline Double_t DispersionVel(int j, Int_t start, Int_t end, Double_t mean);
        /// and phase
        inline Double_t DispersionPhs(int j, Int_t start, Int_t end, Double_t mean);
        /// Calculate the entropy in a given dimension. This can be used as a node splitting criterion
        /// instead of most spread dimension
        /// for positions
        inline Double_t EntropyPos(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni);
        /// and for velocities
        inline Double_t EntropyVel(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni);
        /// and for phase
        inline Double_t EntropyPhs(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni);
        //@}

        /// \name Rearrange and balance the tree
        /// Rearrange the particle order such that all particles with a d'th coordinate value less than
        /// the k'th particle's are lower in index, and vice versa. This function permanently alters
        /// the NBody::System, but it keeps track of the changes.
        //@{
        inline Double_t MedianPos(int d, Int_t k, Int_t start, Int_t end, bool balanced=true);
        /// same as above but with velocities
        inline Double_t MedianVel(int d, Int_t k, Int_t start, Int_t end, bool balanced=true);
        /// same as above but with full phase-space
        inline Double_t MedianPhs(int d, Int_t k, Int_t start, Int_t end, bool balanced=true);
        /// same as above but with possibly a subset of dimensions of full phase space
        /// NOTE Dim DOES NOT DO ANYTHING SPECIAL YET
        //inline Double_t MedianDim(int d, Int_t k, Int_t start, Int_t end, bool balanced=true, Double_t **metric=NULL);
        //@}

        /// \name Density and NN methods
        /// Inline functions used in searches and density calculations
        //@{
        /// smoothing kernel
        inline Double_t Wsm(Double_t r, Int_t i, Int_t size, Double_t delta, Double_t*x);
        /// inline function to get metric spacing around a target particle or position using nearby neighbours
        /// and leaf node boundary
        inline void CalculateMetricSpacing(Int_t target, int treetype, Double_t *smetric);
        inline void CalculateMetricSpacing(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric);
        inline void CalculateMetricSpacing(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric);
        /// Calculate the distribution of mass in the full phase-space around target particle, position.
        inline void CalculateMetricTensor(Int_t target, int treetype, Double_t *smetric, Double_t *metric, GMatrix &gmetric);
        inline void CalculateMetricTensor(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &gmetric);
        inline void CalculateMetricTensor(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &gmetric);
        //inline void CalculateMetricTensor(int target, int treetype, PriorityQueue *pq, GMatrix gmetric);
        ///load data from queue to array also check if search failed.
        inline void LoadNN(const Int_t ns, PriorityQueue *pq, Int_t *nn, Double_t *dist);
        //@}
    };

}
#endif
// KDTREE_H
