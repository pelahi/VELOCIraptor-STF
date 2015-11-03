/*! \file PriorityQueue.h
 *  \brief header file for the PriorityQueue
*/

#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <cstdlib>
#include <iostream>
#include <NBodyMath.h>

namespace NBody
{
    class PriorityQueue
    {
/*!
    \class NBody::PriorityQueue
    \brief A simple array priority queue that orders the values in the queue
*/
        protected:

        /// structure used to order items via a priority queue
        struct queue_member
        {
            Int_t queue;
            Double_t priority;
        };
        Int_t n;
        Int_t max_size;
        queue_member *pq;
        public:
        PriorityQueue(Int_t max)
        {
            n = 0;
            max_size = max;
            pq = new queue_member[max+1];
        }
        ~PriorityQueue() { delete [] pq; }
        bool Empty()
        {
            if (n == 0) return true;
            else return false;
        }
        Int_t Size() { return n; }
        Int_t MaxSize() { return max_size; }
        void Reset() { n = 0; }
        void Push(Int_t p, Double_t dist)
        {
            if (++n > max_size)
            {
                std::cerr << "Priority queue full!" << std::endl;
                exit(EXIT_FAILURE);
            }
            Int_t r = n;
            while (r > 1)
            {
                Int_t p = r>>1;
                if (pq[p].priority >= dist)
                    break;
                pq[r] = pq[p];
                r = p;
            }
            pq[r].priority = dist;
            pq[r].queue = p;
        }
        void Pop()
        {
            Double_t kn = pq[n--].priority;
            Int_t p = 1;
            Int_t r = p<<1;
            while (r <= n)
            {
                if (r < n && pq[r].priority < pq[r+1].priority)
                    r++;
                if (kn >= pq[r].priority)
                    break;
                pq[p] = pq[r];
                p = r;
                r = p<<1;
            }
            pq[p] = pq[n+1];
        }
        Double_t TopPriority() { return pq[1].priority; }
        Int_t TopQueue() { return pq[1].queue; }
    };

    class NPriorityQueue : public PriorityQueue
    {
/*!
    \class NBody::NPriorityQueue
    \brief A priority queue which has additional N associated values besides priority (doesn't quite work yet don't know why)
*/
        private :
        /// array structure used to order items via a priority queue
        struct nqueue {
            Double_t *values;
        };
        Int_t nvalues;
        nqueue *npq;
        public :
        NPriorityQueue(Int_t max, Int_t Nval) : PriorityQueue(max)
        {
            npq = new nqueue[max+1];
            nvalues=Nval;
            for (Int_t i=0;i<max+1;i++) npq[i].values=new Double_t[Nval];
        }
        ~NPriorityQueue()
        { 
            for (Int_t i=0;i<max_size+1;i++) delete[] npq[i].values;
            delete[] npq;
        }
        void Push(Int_t p, Double_t dist, Double_t *pvalue) 
        {
            if (++n > max_size)
            {
                std::cerr << "Priority queue full!" << std::endl;
                exit(EXIT_FAILURE);
            }
            Int_t r = n;
            while (r > 1)
            {
                Int_t p = r>>1;
                if (pq[p].priority >= dist)
                    break;
                pq[r] = pq[p];
                r = p;
                for (int i=0;i<nvalues;i++) npq[r].values[i]=npq[p].values[i];
            }
            pq[r].priority = dist;
            pq[r].queue = p;
            for (int i=0;i<nvalues;i++) npq[r].values[i]=pvalue[i];
        }
        void Pop()
        {
            Double_t kn = pq[n--].priority;
            Int_t p = 1;
            Int_t r = p<<1;
            while (r <= n)
            {
                if (r < n && pq[r].priority < pq[r+1].priority)
                    r++;
                if (kn >= pq[r].priority)
                    break;
                pq[p] = pq[r];
                p = r;
                for (int i=0;i<nvalues;i++) npq[r].values[i]=npq[p].values[i];
                r = p<<1;
            }
            pq[p] = pq[n+1];
            for (int i=0;i<nvalues;i++) npq[p].values[i]=npq[n+1].values[i];
        }
        Double_t TopValue(Int_t N) {
            if (N<nvalues) return npq[1].values[N]; 
            else return 0;
        }
    };
}

#endif // PRIORITYQUEUE_H
