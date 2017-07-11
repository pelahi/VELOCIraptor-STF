/*! \file KDSplitNode.cxx
 *  \brief This file contains implementation of Split Node functions

*/

#include <iostream>
#include <KDNode.h>

namespace NBody
{
    ///Split Node Functions
    ///\todo adjust periodic searches so that reflection is calculated faster.
    ///Non periodic functions
    //@{
    void SplitNode::FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPosition(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FindNearestPos(rd,bucket,pq,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestPos(rd,bucket,pq,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestPos(rd,bucket,pq,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestPos(rd,bucket,pq,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetVelocity(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FindNearestVel(rd,bucket,pq,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestVel(rd,bucket,pq,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestVel(rd,bucket,pq,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestVel(rd,bucket,pq,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FindNearestPhase(rd,bucket,pq,off,target);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestPhase(rd,bucket,pq,off,target);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestPhase(rd,bucket,pq,off,target);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestPhase(rd,bucket,pq,off,target);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, Double_t *metric)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FindNearestMetric(rd,bucket,pq,off,target,metric);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestMetric(rd,bucket,pq,off,target,metric);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestMetric(rd,bucket,pq,off,target,metric);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestMetric(rd,bucket,pq,off,target,metric);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FindNearestMetricwithTensor(rd,bucket,pq,off,target,m0,m1,gm);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestMetricwithTensor(rd,bucket,pq,off,target,m0,m1,gm);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestMetricwithTensor(rd,bucket,pq,off,target,m0,m1,gm);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestMetricwithTensor(rd,bucket,pq,off,target,m0,m1,gm);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Int_t target, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPosition(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FindNearestCriterion(rd,cmp,params,bucket,pq,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestCriterion(rd,cmp,params,bucket,pq,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestCriterion(rd,cmp,params,bucket,pq,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestCriterion(rd,cmp,params,bucket,pq,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }

    void SplitNode::FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = x[cut_dim] - cut_val;
        if (new_off < 0)
        {
            left->FindNearestPos(rd,bucket,pq,off,x, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestPos(rd,bucket,pq,off,x, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestPos(rd,bucket,pq,off,x, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestPos(rd,bucket,pq,off,x, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *v, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = v[cut_dim] - cut_val;
        if (new_off < 0)
        {
            left->FindNearestVel(rd,bucket,pq,off,v, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestVel(rd,bucket,pq,off,v, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestVel(rd,bucket,pq,off,v, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestVel(rd,bucket,pq,off,v, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v)
    {
        Double_t old_off = off[cut_dim];
        Double_t value;
        if (cut_dim<3)value=x[cut_dim];
        else value=v[cut_dim-3];
        Double_t new_off = value - cut_val;
        if (new_off < 0)
        {
            left->FindNearestPhase(rd,bucket,pq,off,x,v);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestPhase(rd,bucket,pq,off,x,v);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestPhase(rd,bucket,pq,off,x,v);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestPhase(rd,bucket,pq,off,x,v);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestPos(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, int dim)
    {
        FindNearestPos(rd,bucket,pq,off,x.GetCoord(),dim);
    }
    void SplitNode::FindNearestVel(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate v, int dim)
    {
        FindNearestVel(rd,bucket,pq,off,v.GetCoord(),dim);
    }
    void SplitNode::FindNearestPhase(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v)
    {
        FindNearestPhase(rd,bucket,pq,off,x.GetCoord(),v.GetCoord());
    }

    void SplitNode::FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t *metric)
    {
        Double_t old_off = off[cut_dim];
        Double_t value;
        if (cut_dim<3)value=x[cut_dim];
        else value=v[cut_dim-3];
        Double_t new_off = value - cut_val;
        if (new_off < 0)
        {
            left->FindNearestMetric(rd,bucket,pq,off,x,v,metric);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestMetric(rd,bucket,pq,off,x,v,metric);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestMetric(rd,bucket,pq,off,x,v,metric);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestMetric(rd,bucket,pq,off,x,v,metric);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        Double_t old_off = off[cut_dim];
        Double_t value;
        if (cut_dim<3)value=x[cut_dim];
        else value=v[cut_dim-3];
        Double_t new_off = value - cut_val;
        if (new_off < 0)
        {
            left->FindNearestMetricwithTensor(rd,bucket,pq,off,x,v,m0,m1,gm);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestMetricwithTensor(rd,bucket,pq,off,x,v,m0,m1,gm);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestMetricwithTensor(rd,bucket,pq,off,x,v,m0,m1,gm);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestMetricwithTensor(rd,bucket,pq,off,x,v,m0,m1,gm);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::FindNearestMetric(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t *metric)
    {
        FindNearestMetric(rd, bucket, pq, off, x.GetCoord(), v.GetCoord(), metric);
    }
    void SplitNode::FindNearestMetricwithTensor(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t* off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        FindNearestMetricwithTensor(rd, bucket, pq, off, x.GetCoord(), v.GetCoord(), m0, m1, gm);
    }
    void SplitNode::FindNearestCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t* off, Particle &p, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = p.GetPosition(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FindNearestCriterion(rd,cmp,params,bucket,pq,off,p, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                right->FindNearestCriterion(rd,cmp,params,bucket,pq,off,p, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FindNearestCriterion(rd,cmp,params,bucket,pq,off,p, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < pq->TopPriority())
            {
                off[cut_dim] = new_off;
                left->FindNearestCriterion(rd,cmp,params,bucket,pq,off,p, dim);
                off[cut_dim] = old_off;
            }
        }
    }

    void SplitNode::SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t* off, Int_t target, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                right->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                left->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }

    void SplitNode::SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t* off, Double_t *x, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = x[cut_dim] - cut_val;
        if (new_off < 0)
        {
            left->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,x,dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                right->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,x,dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,x,dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                left->SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,x,dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::SearchBallPos(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t* off, Coordinate x, int dim)
    {
		SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,x.GetCoord(),dim);
	}


    void SplitNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Int_t target, Int_t &nt, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,target,nt,dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                right->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,target,nt,dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,target,nt,dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                left->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,target,nt,dim);
                off[cut_dim] = old_off;
            }
        }
    }

    void SplitNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Double_t *x, Int_t &nt, int dim)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = x[cut_dim] - cut_val;
        if (new_off < 0)
        {
            left->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x,nt,dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                right->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x,nt,dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x,nt,dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                left->SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x,nt,dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::SearchBallPosTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t* off, Coordinate x, Int_t &nt,int dim)
    {
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x.GetCoord(),nt,dim);
    }


	void SplitNode::SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t* off, Int_t target, int dim)
    {
        //assume some distance measure is stored in params[1]
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                right->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                left->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::SearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t* off, Particle &target, int dim)
    {
        //assume some distance measure is stored in params[1]
        Double_t old_off = off[cut_dim];
        Double_t new_off = target.GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                right->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                left->SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }

    void SplitNode::SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Int_t target, int dim)
    {
        //assume some distance measure is stored in params[1]
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                right->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                left->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::SearchCriterionTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t* off, Particle &target, int dim)
    {
        //assume some distance measure is stored in params[1]
        Double_t old_off = off[cut_dim];
        Double_t new_off = target.GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                right->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                left->SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Int_t target, int dim)
    {
        //assume some distance measure is stored in params[1]
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                right->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                left->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }
    void SplitNode::SearchCriterionNoDist(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t* off, Particle &target, int dim)
    {
        //assume some distance measure is stored in params[1]
        Double_t old_off = off[cut_dim];
        Double_t new_off = target.GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                right->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < params[1])
            {
                off[cut_dim] = new_off;
                left->SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,target, dim);
                off[cut_dim] = old_off;
            }
        }
    }

    void SplitNode::FOFSearchBall(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FOFSearchBall(rd,fdist2,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                right->FOFSearchBall(rd,fdist2,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FOFSearchBall(rd,fdist2,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            rd += -old_off*old_off + new_off*new_off;
            if (rd < fdist2)
            {
                off[cut_dim] = new_off;
                left->FOFSearchBall(rd,fdist2,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
                off[cut_dim] = old_off;
            }
        }
    }

    //key here is params which tell one how to search the tree
    void SplitNode::FOFSearchCriterion(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            if ((int)params[0]==0) rd += (-old_off*old_off + new_off*new_off)/params[1];
            else if ((int)params[0]==1) rd += (-old_off*old_off + new_off*new_off)/params[2];
            else if ((int)params[0]==2) rd += (-old_off*old_off + new_off*new_off)/params[(cut_dim<3)*1+(cut_dim>=3)*2];
            if (rd < 1)
            {
                off[cut_dim] = new_off;
                right->FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            if ((int)params[0]==0) rd += (-old_off*old_off + new_off*new_off)/params[1];
            else if ((int)params[0]==1) rd += (-old_off*old_off + new_off*new_off)/params[2];
            else if ((int)params[0]==2) rd += (-old_off*old_off + new_off*new_off)/params[(cut_dim<3)*1+(cut_dim>=3)*2];
            if (rd < 1)
            {
                off[cut_dim] = new_off;
                left->FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
                off[cut_dim] = old_off;
            }
        }
    }

    //key here is params which tell one how to search the tree
    void SplitNode::FOFSearchCriterionSetBasisForLinks(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Int_t target)
    {
        Double_t old_off = off[cut_dim];
        Double_t new_off = bucket[target].GetPhase(cut_dim) - cut_val;
        if (new_off < 0)
        {
            left->FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            if ((int)params[0]==0) rd += (-old_off*old_off + new_off*new_off)/params[1];
            else if ((int)params[0]==1) rd += (-old_off*old_off + new_off*new_off)/params[2];
            else if ((int)params[0]==2) rd += (-old_off*old_off + new_off*new_off)/params[(cut_dim<3)*1+(cut_dim>=3)*2];
            if (rd < 1)
            {
                off[cut_dim] = new_off;
                right->FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
                off[cut_dim] = old_off;
            }
        }
        else
        {
            right->FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            if ((int)params[0]==0) rd += (-old_off*old_off + new_off*new_off)/params[1];
            else if ((int)params[0]==1) rd += (-old_off*old_off + new_off*new_off)/params[2];
            else if ((int)params[0]==2) rd += (-old_off*old_off + new_off*new_off)/params[(cut_dim<3)*1+(cut_dim>=3)*2];
            if (rd < 1)
            {
                off[cut_dim] = new_off;
                left->FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
                off[cut_dim] = old_off;
            }
        }
    }

    //@}

    ///Periodic version that generate reflections and use non periodic searches
    //@{
    //there are several reflections to search (sum_i^ND choose(ND,i), so for 3d have 7 possible reflections)
    void SplitNode::FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Double_t sval;
        Coordinate x0,xp;
        x0=Coordinate(bucket[target].GetPosition());
        FindNearestPosPeriodic(rd,bucket,pq,off,p,x0,dim);
    }
    void SplitNode::FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        FindNearestVel(rd,bucket,pq,off,target,dim);
    }
    void SplitNode::FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target)
    {
        Coordinate x0,xp,v;
        x0=Coordinate(bucket[target].GetPosition());
        v=Coordinate(bucket[target].GetVelocity());
        FindNearestPhasePeriodic(rd,bucket,pq,off,p,x0,v);
    }
    void SplitNode::FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, Double_t *metric)
    {
        Coordinate x0,v;
        x0=Coordinate(bucket[target].GetPosition());
        v=Coordinate(bucket[target].GetVelocity());
        FindNearestMetricPeriodic(rd,bucket,pq,off,p,x0,v,metric);
    }
    void SplitNode::FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        Coordinate x0,xp,v;
        x0=Coordinate(bucket[target].GetPosition());
        v=Coordinate(bucket[target].GetVelocity());
        FindNearestMetricwithTensorPeriodic(rd,bucket,pq,off,p,x0,v,m0,m1,gm);
    }
    void SplitNode::FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle p0;
        p0=bucket[target];
        FindNearestCriterionPeriodic(rd,cmp,params,bucket,pq,off,p,p0,dim);
    }

    void SplitNode::FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, int dim)
    {
        //there are several reflections to search (sum_i^ND choose(ND,i), so for 3d have 7 possible reflections)
        //reset offset search for each new search
        Double_t sval;
        Coordinate x0(x),xp;
        FindNearestPos(rd,bucket,pq,off,x0,dim);
        for (int k=0;k<dim;k++) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection1D(x0,xp,p,k);
            if (sqrt(pq->TopPriority())>sval) FindNearestPos(rd,bucket,pq,off,xp,dim);
        }
        if (dim==3) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,1);
            if (pq->TopPriority()>sval) FindNearestPos(rd,bucket,pq,off,xp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,2);
            if (pq->TopPriority()>sval) FindNearestPos(rd,bucket,pq,off,xp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,1,2);
            if (pq->TopPriority()>sval) FindNearestPos(rd,bucket,pq,off,xp,dim);
        }
        // search all axis if current max dist less than search radius
        if (dim>1) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflectionND(x0,xp,p,dim);
            if (pq->TopPriority()>sval) FindNearestPos(rd,bucket,pq,off,xp,dim);
        }
    }
    void SplitNode::FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *v, int dim)
    {
        FindNearestVel(rd,bucket,pq,off,v,dim);
    }
    void SplitNode::FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, Double_t *v0)
    {
        Coordinate x0(x),xp,v(v0);
        FindNearestPhase(rd,bucket,pq,off,x0,v);
        for (int k=0;k<NSPACEDIM;k++) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection1D(x0,xp,p,k);
            FindNearestPhase(rd,bucket,pq,off,xp,v);
        }
        if (NSPACEDIM==3) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,0,1);
            FindNearestPhase(rd,bucket,pq,off,xp,v);
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,0,2);
            FindNearestPhase(rd,bucket,pq,off,xp,v);
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,1,2);
            FindNearestPhase(rd,bucket,pq,off,xp,v);
        }
        // search all axis if current max dist less than search radius
        if (NSPACEDIM>1) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflectionND(x0,xp,p,NSPACEDIM);
            FindNearestPhase(rd,bucket,pq,off,xp,v);
        }
    }
    void SplitNode::FindNearestPosPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, int dim)
    {
        FindNearestPosPeriodic(rd,bucket,pq,off,p,x.GetCoord(),dim);
    }
    void SplitNode::FindNearestVelPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate v, int dim)
    {
        FindNearestVelPeriodic(rd,bucket,pq,off,p,v.GetCoord(),dim);
    }
    void SplitNode::FindNearestPhasePeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, Coordinate v)
    {
        FindNearestPhasePeriodic(rd,bucket,pq,off,p,x.GetCoord(),v.GetCoord());
    }

    void SplitNode::FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, Double_t *v0, Double_t *metric)
    {
        Coordinate x0(x),xp,v(v0);
        FindNearestMetric(rd,bucket,pq,off,x0,v,metric);
        for (int k=0;k<NSPACEDIM;k++) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection1D(x0,xp,p,k);
            FindNearestMetric(rd,bucket,pq,off,xp,v,metric);
        }
        if (NSPACEDIM==3) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,0,1);
            FindNearestMetric(rd,bucket,pq,off,xp,v,metric);
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,0,2);
            FindNearestMetric(rd,bucket,pq,off,xp,v,metric);
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,1,2);
            FindNearestMetric(rd,bucket,pq,off,xp,v,metric);
        }
        // search all axis if current max dist less than search radius
        if (NSPACEDIM>1) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflectionND(x0,xp,p,NSPACEDIM);
            FindNearestMetric(rd,bucket,pq,off,xp,v,metric);
        }
    }
    void SplitNode::FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Double_t *x, Double_t *v0, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        Coordinate x0(x),xp,v(v0);
        FindNearestMetricwithTensor(rd,bucket,pq,off,x0,v,m0,m1,gm);
        for (int k=0;k<NSPACEDIM;k++) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection1D(x0,xp,p,k);
            FindNearestMetricwithTensor(rd,bucket,pq,off,xp,v,m0,m1,gm);
        }
        if (NSPACEDIM==3) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,0,1);
            FindNearestMetricwithTensor(rd,bucket,pq,off,xp,v,m0,m1,gm);
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,0,2);
            FindNearestMetricwithTensor(rd,bucket,pq,off,xp,v,m0,m1,gm);
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflection2D(x0,xp,p,1,2);
            FindNearestMetricwithTensor(rd,bucket,pq,off,xp,v,m0,m1,gm);
        }
        // search all axis if current max dist less than search radius
        if (NSPACEDIM>1) {
            for (int j = 0; j < NPHASEDIM; j++) off[j] = 0.0;
            PeriodicReflectionND(x0,xp,p,NSPACEDIM);
            FindNearestMetricwithTensor(rd,bucket,pq,off,xp,v,m0,m1,gm);
        }
    }
    void SplitNode::FindNearestMetricPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, Coordinate v, Double_t *metric)
    {
        FindNearestMetricPeriodic(rd, bucket, pq, off, p, x.GetCoord(), v.GetCoord(), metric);
    }
    void SplitNode::FindNearestMetricwithTensorPeriodic(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)
    {
        FindNearestMetricwithTensorPeriodic(rd, bucket, pq, off, p, x.GetCoord(), v.GetCoord(), m0, m1, gm);
    }
    void SplitNode::FindNearestCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        Particle pp;
        pp=p0;
        FindNearestCriterion(rd,cmp,params,bucket,pq,off,p0,dim);
        for (int k=0;k<dim;k++) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection1D(p0,pp,p,k);
            FindNearestCriterion(rd,cmp,params,bucket,pq,off,pp,dim);
        }
        if (dim==3) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,1);
            FindNearestCriterion(rd,cmp,params,bucket,pq,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,2);
            FindNearestCriterion(rd,cmp,params,bucket,pq,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,1,2);
            FindNearestCriterion(rd,cmp,params,bucket,pq,off,pp,dim);
        }
        // search all axis if current max dist less than search radius
        if (dim>1) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflectionND(p0,pp,p,dim);
            FindNearestCriterion(rd,cmp,params,bucket,pq,off,pp,dim);
        }
    }

    void SplitNode::SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Coordinate x0;
        x0=Coordinate(bucket[target].GetPosition());
        SearchBallPosPeriodic(rd,fdist2,iGroup,bucket,Group,dist2,off,p,x0,dim);
    }

    void SplitNode::SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t *off, Double_t *p, Double_t *x, int dim)
    {
        Double_t sval;
        Coordinate x0(x),xp;
        SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,x0,dim);
        for (int k=0;k<dim;k++) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection1D(x0,xp,p,k);
            if (fdist2>sval*sval) SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,xp,dim);
        }
        if (dim==3) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,1);
            if (fdist2>sval*sval) SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,xp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,2);
            if (fdist2>sval*sval) SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,xp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,1,2);
            if (fdist2>sval*sval) SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,xp,dim);
        }
        // search all axis if current max dist less than search radius
        if (dim>1) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflectionND(x0,xp,p,dim);
            if (fdist2>sval*sval) SearchBallPos(rd,fdist2,iGroup,bucket,Group,dist2,off,xp,dim);
        }
    }
    void SplitNode::SearchBallPosPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t *off, Double_t *p, Coordinate x, int dim)
    {
        SearchBallPosPeriodic(rd,fdist2,iGroup,bucket,Group,dist2,off,p,x.GetCoord(),dim);
    }

    void SplitNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *p, Int_t target, Int_t &nt, int dim)
    {
        Coordinate x0;
        x0=Coordinate(bucket[target].GetPosition());
        SearchBallPosPeriodicTagged(rd,fdist2,bucket,tagged,off,p,x0,nt,dim);
    }

    void SplitNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *p, Double_t *x, Int_t &nt, int dim)
    {
        Double_t sval;
        Coordinate x0(x),xp;
        SearchBallPosTagged(rd,fdist2,bucket,tagged,off,x0,nt,dim);
        for (int k=0;k<dim;k++) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection1D(x0,xp,p,k);
            if (fdist2>sval*sval) SearchBallPosTagged(rd,fdist2,bucket,tagged,off,xp,nt,dim);
        }
        if (dim==3) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,1);
            if (fdist2>sval*sval) SearchBallPosTagged(rd,fdist2,bucket,tagged,off,xp,nt,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,2);
            if (fdist2>sval*sval) SearchBallPosTagged(rd,fdist2,bucket,tagged,off,xp,nt,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,1,2);
            if (fdist2>sval*sval) SearchBallPosTagged(rd,fdist2,bucket,tagged,off,xp,nt,dim);
        }
        // search all axis if current max dist less than search radius
        if (dim>1) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            sval=PeriodicReflectionND(x0,xp,p,dim);
            if (fdist2>sval*sval) SearchBallPosTagged(rd,fdist2,bucket,tagged,off,xp,nt,dim);
        }
    }
    void SplitNode::SearchBallPosPeriodicTagged(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *p, Coordinate x, Int_t &nt, int dim)
    {
        SearchBallPosPeriodicTagged(rd,fdist2,bucket,tagged,off,p,x.GetCoord(),nt,dim);
    }

    void SplitNode::SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle x0;
        x0=bucket[target];
        SearchCriterionPeriodic(rd,cmp,params,iGroup,bucket,Group,dist2,off,p,x0,dim);
    }
    void SplitNode::SearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *dist2, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        Particle pp;
        pp=p0;
        SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,p0,dim);
        for (int k=0;k<dim;k++) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection1D(p0,pp,p,k);
            SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,pp,dim);
        }
        if (dim==3) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,1);
            SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,2);
            SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,1,2);
            SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,pp,dim);
        }
        // search all axis if current max dist less than search radius
        if (dim>1) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflectionND(p0,pp,p,dim);
            SearchCriterion(rd,cmp,params,iGroup,bucket,Group,dist2,off,pp,dim);
        }
    }
    void SplitNode::SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle x0;
        x0=bucket[target];
        SearchCriterionPeriodicTagged(rd,cmp,params,bucket,nt,tagged,off,p,x0,dim);
    }
    void SplitNode::SearchCriterionPeriodicTagged(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &nt, Int_t *tagged, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        Particle pp;
        pp=p0;
        SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,p0,dim);
        for (int k=0;k<dim;k++) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection1D(p0,pp,p,k);
            SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,pp,dim);
        }
        if (dim==3) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,1);
            SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,2);
            SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,1,2);
            SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,pp,dim);
        }
        // search all axis if current max dist less than search radius
        if (dim>1) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflectionND(p0,pp,p,dim);
            SearchCriterionTagged(rd,cmp,params,bucket,nt,tagged,off,pp,dim);
        }
    }

    void SplitNode::SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *p, Int_t target, int dim)
    {
        Particle x0;
        x0=bucket[target];
        SearchCriterionNoDistPeriodic(rd,cmp,params,iGroup,bucket,Group,off,p,x0,dim);
    }

    void SplitNode::SearchCriterionNoDistPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *p, Particle &p0, int dim)
    {
        Particle pp;
        pp=p0;
        SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,p0,dim);
        for (int k=0;k<dim;k++) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection1D(p0,pp,p,k);
            SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,pp,dim);
        }
        if (dim==3) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,1);
            SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,0,2);
            SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,pp,dim);
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflection2D(p0,pp,p,1,2);
            SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,pp,dim);
        }
        // search all axis if current max dist less than search radius
        if (dim>1) {
            for (int j = 0; j < dim; j++) off[j] = 0.0;
            PeriodicReflectionND(p0,pp,p,dim);
            SearchCriterionNoDist(rd,cmp,params,iGroup,bucket,Group,off,pp,dim);
        }
    }
    //here code is effectively like that of FOFsearchBallPeriodic but adjust the particle's position in the left/right search
    void SplitNode::FOFSearchBallPeriodic(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *p, Int_t target)
    {
        //first search normal particle
        FOFSearchBall(rd, fdist2, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
        Coordinate x0(bucket[target].GetPosition()),xp;
        Double_t sval;
        for (int k=0;k<NSPACEDIM;k++) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection1D(x0,xp,p,k);
            if (fdist2>sval*sval) {
                for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
                FOFSearchBall(rd, fdist2, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
            }
        }
        if (NSPACEDIM==3) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,1);
            if (fdist2>sval*sval) {
                for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
                FOFSearchBall(rd, fdist2, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
            }
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,2);
            if (fdist2>sval*sval) {
                for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
                FOFSearchBall(rd, fdist2, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
            }
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,1,2);
            if (fdist2>sval*sval) {
                for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
                FOFSearchBall(rd, fdist2, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
            }
        }
        // search all axis if current max dist less than search radius
        if (NSPACEDIM>1) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflectionND(x0,xp,p,NSPACEDIM);
            if (fdist2>sval*sval) {
                for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
                FOFSearchBall(rd, fdist2, iGroup, nActive, bucket, Group, Len, Head, Tail, Next, BucketFlag, Fifo, iTail, off, target);
            }
        }
        for (int j=0;j<3;j++) bucket[target].SetPosition(j,x0[j]);
    }

    //key here is params which tell one how to search the tree
    void SplitNode::FOFSearchCriterionPeriodic(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *p, Int_t target)
    {
        FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        Coordinate x0(bucket[target].GetPosition()),xp;
        Double_t sval;
        for (int k=0;k<NSPACEDIM;k++) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection1D(x0,xp,p,k);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        }
        if (NSPACEDIM==3) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,1);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,2);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,1,2);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        }
        // search all axis if current max dist less than search radius
        if (NSPACEDIM>1) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflectionND(x0,xp,p,NSPACEDIM);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterion(rd,cmp,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        }
        for (int j=0;j<3;j++) bucket[target].SetPosition(j,x0[j]);
    }

    void SplitNode::FOFSearchCriterionSetBasisForLinksPeriodic(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &iTail, Double_t* off, Double_t *p, Int_t target)
    {
        FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        Coordinate x0(bucket[target].GetPosition()),xp;
        Double_t sval;
        for (int k=0;k<NSPACEDIM;k++) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection1D(x0,xp,p,k);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        }
        if (NSPACEDIM==3) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,1);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,0,2);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflection2D(x0,xp,p,1,2);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        }
        // search all axis if current max dist less than search radius
        if (NSPACEDIM>1) {
            for (int j = 0; j < NSPACEDIM; j++) off[j] = 0.0;
            sval=PeriodicReflectionND(x0,xp,p,NSPACEDIM);
            for (int j=0;j<3;j++) bucket[target].SetPosition(j,xp[j]);
            FOFSearchCriterionSetBasisForLinks(rd,cmp,check,params,iGroup,nActive,bucket,Group,Len,Head,Tail,Next,BucketFlag,Fifo,iTail,off,target);
        }
        for (int j=0;j<3;j++) bucket[target].SetPosition(j,x0[j]);
    }

    //@}
}
