/*! \file KDNode.h
 *  \brief header file for the \ref NBody::Node Class

*/

#ifndef KDNODE_H
#define KDNODE_H

#include <NBody.h>
#include <NBodyMath.h>
#include <PriorityQueue.h>
#include <DistFunc.h>
#include <FOFFunc.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

///this defines the number of spatial dimensions and is useful for phase-space searches where one might only consider projected
///coordinates
#define NSPACEDIM 3
#define NPHASEDIM 6

///if large tree is declared then assumes tree will span > MAXINT particles and
///need to use long ints to store information. Otherwise, no need to have long ints
#ifdef LARGETREE
typedef long int Int_tree_t;
typedef unsigned long int UInt_tree_t;
#else
typedef int Int_tree_t;
typedef unsigned int UInt_tree_t;
#endif

namespace NBody
{

/*!
    \class NBody::Node
    \brief Base virtual class for a node used by \ref NBody::KDTree.

    This is a virtual class from which \ref NBody::SplitNode and \ref NBody::LeafNode are derived from.
*/
    class Node
    {
        protected:
        /// id of Node
        UInt_tree_t nid;
        /// boundaries of node
        DoublePos_t xbnd[6][2];
        /// number of particles associated with node
        UInt_tree_t count;

        UInt_tree_t bucket_start;
        UInt_tree_t bucket_end;
        unsigned short numdim;
        public:
        virtual ~Node() {};

        /// \name Simple Get functions
        //@{
        ///Get Node ID
        virtual Int_t GetID(){return nid;}
        ///Get boundary of volume enclosed by node
        virtual Double_t GetBoundary(int i, int j){return xbnd[i][j];}
        ///Get total number of particles in node
        virtual Int_t GetCount(){return count;}
        ///Get start index in particle array of particles enclosed by node
        virtual Int_t GetStart(){return bucket_start;}
        ///Get end index in particle array of particles enclosed by node
        virtual Int_t GetEnd(){return bucket_end;}
        //@}

        /// \name Simple Set functions
        //@{
        /// set Id --- use with caution
        virtual void SetID(Int_tree_t id){nid=id;}
        //@}

        /// \name Find Nearest routines:
        /// Find nearest particles to target t or position x and fills up the priority queue. Note velocities and positions can be projections
        /// The basic idea of this calls is to first see if the node is a split node. If its a split node, then the cut value is compared to
        /// the coordinates given. If it is smaller, then left node is searched otherwise, right node is searched.
        /// By moving up the tree, eventually find a leaf node at which point particles in the tree are compared to the coordinate.
        /// Key to finding nearest is to also compare the maximum current distance relative to the position is to the left or right
        /// of the current split nodes cut value. This form of walking the tree is also used in \ref SearchBallPos and \ref FOFSearchBall
        /// and like routines.
        //@{
        virtual void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t) = 0;
        //Same as above but find nearest to position x (or v or x,v) with interfaces for both Double_t pointers and Coordinates
        virtual void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, int dim=3) = 0;
        virtual void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *v, int dim=3) = 0;
        virtual void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v) = 0;
        virtual void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, int dim=3) = 0;
        virtual void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate v, int dim=3) = 0;
        virtual void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v) = 0;
        //Here also apply Criterion to NN under assumption that looking at physical search
        virtual void FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Particle &p, int dim=3) = 0;

        ///phase-space search but first is metric scaling
        virtual void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, Double_t * metric) = 0;
        ///phase-space search with full metric tensor to get distance with g_(mu,nu) dx^mu dx^nu
        virtual void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm) = 0;
        //same as above but different iterface
        virtual void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t * metric) = 0;
        virtual void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t * metric) = 0;
        virtual void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm) = 0;
        virtual void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm) = 0;

        //Same as above but for periodic system. Reason for duplication is to not have to check if system is periodic every time a node is checked
        virtual void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3) = 0;
        virtual void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3) = 0;
        virtual void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t) = 0;
        virtual void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3) = 0;
        virtual void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3) = 0;
        virtual void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v) = 0;
        virtual void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3) = 0;
        virtual void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3) = 0;
        virtual void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v) = 0;
        virtual void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t * metric) = 0;
        virtual void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm) = 0;
        virtual void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t * metric) = 0;
        virtual void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t * metric) = 0;
        virtual void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm) = 0;
        virtual void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm) = 0;
        virtual void FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *period, Int_t t, int dim=3)=0;
        virtual void FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *period, Particle &p, int dim=3)=0;
        //@}

        /// \name Ball Searches
        /// Searches nodes about a point or particle using a given distance. To be used with \ref NBody::KDTree::SearchBall like routines
        //@{
        virtual void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *x, int dim=3) = 0;
        virtual void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Coordinate x, int dim=3) = 0;

        //Same as above but for periodic system. Reason for duplication is to not have to check if system is periodic every time a node is checked
        virtual void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Int_t t, int dim=3) = 0;
        virtual void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Double_t *x, int dim=3) = 0;
        virtual void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Coordinate x, int dim=3) = 0;

        //Same as above but just returns number of tagged particles and an unordered array storing the index of tagged particles
        virtual void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Int_t t, Int_t &nt, int dim=3) = 0;
        virtual void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *x, Int_t &nt, int dim=3) = 0;
        virtual void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Coordinate x, Int_t &nt, int dim=3) = 0;

        //Same as above but for periodic system. Reason for duplication is to not have to check if system is periodic every time a node is checked
        virtual void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Int_t t, Int_t &nt, int dim=3) = 0;
        virtual void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Double_t *x, Int_t &nt, int dim=3) = 0;
        virtual void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Coordinate x, Int_t &nt, int dim=3) = 0;


        //same as above but for criterion
        //as we do not know the criterion a priori, for this to work also implement one with particle interface along with index
        virtual void SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Particle &p, int dim=3) = 0;
        virtual void SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Int_t t, int dim=3) = 0;
        virtual void SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Particle &p, int dim=3) = 0;
/*
        //same but also stores unsorted array of indices of tagged particles
        virtual void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Int_t &nt, Int_t *tagged, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Int_t &nt, Int_t *tagged, Double_t* off, Particle &p, int dim=3) = 0;
        virtual void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Int_t t, int dim=3) = 0;
        virtual void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Particle &p, int dim=3) = 0;
*/
        //same but also stores unsorted array of indices of tagged particles
        virtual void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Particle &p, int dim=3) = 0;
        virtual void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Int_t t, int dim=3) = 0;
        virtual void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Particle &p, int dim=3) = 0;

        //same as above but does not store a distance, just sets nn value to iGroup
        virtual void SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Int_t t, int dim=3) = 0;
        virtual void SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Particle &p, int dim=3) = 0;
        virtual void SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Double_t *period, Int_t t, int dim=3) = 0;
        virtual void SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Double_t *period, Particle &p, int dim=3) = 0;

        //@}

        /// \name FOF Searches
        /// Like Ball searches but doesn't store distances and is meant to be used with the \ref NBody::KDTree::FOF routines. \n
        /// The node routine finds all nActive particles within a distance^2 fdist2, marks all particles using
        /// their IDS and Group array.
        /// Also pass a Head, Tail and Next array which are adjusted to point to the head and tail of group and Next array
        /// is used to traverse the group. BucketFlag reduces the number of buckets searched
        /// There is FOF using Criterion and also one where Criterion is passed and another
        /// where criterion is passed and only certain particles can be used to generate links to other particles. This means additional checks
        /// as particles can meet the criterion of belonging to two different groups and these groups will not be joined as the particle bridging
        /// the groups cannot be used to form new links.
        //@{

        virtual void FOFSearchBall(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target) = 0;
        virtual void FOFSearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target) = 0;
        virtual void FOFSearchCriterionSetBasisForLinks(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target) = 0;

        //same as above but periodic
        virtual void FOFSearchBallPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target) = 0;
        virtual void FOFSearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target) = 0;
        virtual void FOFSearchCriterionSetBasisForLinksPeriodic(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target) = 0;
        //@}
    };

/*!
    \class NBody::SplitNode
    \brief A node in the kd-tree used by \ref NBody::KDTree  which contains splitting information.

    Minimally, this is a splitting dimension and value, as well as
    pointers to this nodes children.
*/
    class SplitNode : public Node
    {
        private:
        int cut_dim;
        Double_t cut_val;
        Node *left;
        Node *right;
        public:
        SplitNode(Int_t id, int d, Double_t p, Int_t Count, Double_t bnd[6][2], Int_t new_bucket_start, Int_t new_bucket_end, unsigned short ndim, Node *initial_left = NULL,
                  Node *initial_right = NULL)
        {
            nid=id;
            cut_dim = d;
            cut_val = p;
            count=Count;
            bucket_start = new_bucket_start;
            bucket_end = new_bucket_end;
            left = initial_left;
            right = initial_right;
            numdim=ndim;
            for (int j=0;j<numdim;j++) {xbnd[j][0]=bnd[j][0];xbnd[j][1]=bnd[j][1];}
        }
        ~SplitNode() { delete left; delete right; }

        /// \name Simple Get functions
        //@{
        ///Get the cut value
        Double_t GetCutValue(){return cut_val;}
        ///and get the cut dimension
        int GetCutDim(){return cut_dim;}
        ///get the child node to left of cut value
        Node *GetLeft(){return left;}
        ///get the child node to the right of the cut value
        Node *GetRight(){return right;}
        //@}

        //implementations of Find functions
        void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3);
        void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3);
        void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t);
        void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, int dim=3);
        void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *v, int dim=3);
        void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v);
        void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, int dim=3);
        void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate v, int dim=3);
        void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v);
        void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, Double_t * metric);
        void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t * metric);
        void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t * metric);
        void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3);
        void FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Particle &p, int dim=3);

        void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3);
        void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3);
        void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t);
        void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3);
        void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3);
        void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v);
        void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3);
        void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3);
        void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v);
        void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t * metric);
        void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t * metric);
        void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t * metric);
        void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *period, Particle &p, int dim=3);

        //implementation of Ball and FOF searches
        void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t t, int dim=3);
        void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *x, int dim=3);
        void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Coordinate x, int dim=3);

        void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Double_t *x, int dim=3);
        void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Coordinate x, int dim=3);

        void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Int_t t, Int_t &nt, int dim=3);
        void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *x, Int_t &nt, int dim=3);
        void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Coordinate x, Int_t &nt, int dim=3);

        void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Int_t t, Int_t &nt, int dim=3);
        void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Double_t *x, Int_t &nt, int dim=3);
        void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Coordinate x, Int_t &nt, int dim=3);

        void SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t t, int dim=3);
        void SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Particle &p, int dim=3);
        void SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Particle &p, int dim=3);

        void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Int_t t, int dim=3);
        void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Particle &p, int dim=3);
        void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Particle &p, int dim=3);

        void SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Int_t t, int dim=3);
        void SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Particle &p, int dim=3);
        void SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Double_t *period, Particle &p, int dim=3);

        //fof searches
        void FOFSearchBall(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target);
        void FOFSearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target);
        void FOFSearchCriterionSetBasisForLinks(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target);

        void FOFSearchBallPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target);
        void FOFSearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target);
        void FOFSearchCriterionSetBasisForLinksPeriodic(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target);
    };

/*!
    \class NBody::LeafNode
    \brief A node in the kd-tree used by \ref NBody::KDTree representing the end of a branch of a tree.

    It contains pointers to the particles in contains (effectively amounts in instance of the node virtual class)
*/
    class LeafNode : public Node
    {
        private:
        public:
        LeafNode(Int_t id, Int_t new_bucket_start, Int_t new_bucket_end, Double_t bnd[6][2], unsigned short ndim)
        {
            nid=id;
            bucket_start = new_bucket_start;
            bucket_end = new_bucket_end;
            count=bucket_end-bucket_start;
            numdim=ndim;
            for (int j=0;j<numdim;j++) {xbnd[j][0]=bnd[j][0];xbnd[j][1]=bnd[j][1];}
        }
        ~LeafNode() { }

        //implementations of Find functions
        void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3);
        void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3);
        void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t);
        void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, int dim=3);
        void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *v, int dim=3);
        void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v);
        void FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, int dim=3);
        void FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate v, int dim=3);
        void FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v);
        void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, Double_t * metric);
        void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t * metric);
        void FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t * metric);
        void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t t, int dim=3);
        void FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Particle &p, int dim=3);

        void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3);
        void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3);
        void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t);
        void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3);
        void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3);
        void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v);
        void FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3);
        void FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3);
        void FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v);
        void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t * metric);
        void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t * metric);
        void FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t * metric);
        void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm);
        void FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *period, Particle &p, int dim=3);

        //implementation of Ball and FOF searches
        void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t t, int dim=3);
        void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *x, int dim=3);
        void SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Coordinate x, int dim=3);

        void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Double_t *x, int dim=3);
        void SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Coordinate x, int dim=3);

        void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Int_t t, Int_t &nt, int dim=3);
        void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *x, Int_t &nt, int dim=3);
        void SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Coordinate x, Int_t &nt, int dim=3);

        void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Int_t t, Int_t &nt, int dim=3);
        void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Double_t *x, Int_t &nt, int dim=3);
        void SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *period, Coordinate x, Int_t &nt, int dim=3);

        void SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Int_t t, int dim=3);
        void SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Particle &p, int dim=3);
        void SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t* off, Double_t *period, Particle &p, int dim=3);

        void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Int_t t, int dim=3);
        void SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Particle &p, int dim=3);
        void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Double_t *period, Particle &p, int dim=3);

        void SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Int_t t, int dim=3);
        void SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Particle &p, int dim=3);
        void SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Double_t *period, Int_t t, int dim=3);
        void SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Double_t *period, Particle &p, int dim=3);

        //fof searches
        void FOFSearchBall(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target);
        void FOFSearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target);
        void FOFSearchCriterionSetBasisForLinks(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target);

        void FOFSearchBallPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target);
        void FOFSearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target);
        void FOFSearchCriterionSetBasisForLinksPeriodic(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *period, Int_t target);
    };

}
#endif // KDNODE_H
