/*! \file KDFindNearest.cxx
 *  \brief This file contains subroutines involving finding the nearest neighbours 
 
*/

#include <KDTree.h>

namespace NBody
{

    inline void KDTree::LoadNN(const Int_t ns, PriorityQueue *pq, Int_t *nn, Double_t *d)
    {
        for (Int_t i = ns-1; i >=0; i--)
        {
            nn[i] = pq->TopQueue();
            if (nn[i] == -1)
            {
                printf("FindNearest failed for unknown reasons\n");
//                exit(1);
            }
            d[i] = pq->TopPriority();
            pq->Pop();
        }
    }

    /// \name Inline functions used when examining phase-space volume
    //@{
    /// Calculate the local metric describing the mass distribution of the local phase-space volume using the leaf node containing particle position
    /// For this to work effectively, leaf nodes should contain few particles (but more than one).
    /// In Enbid, where leaf nodes contain one particle, the idea is to find the enclosing leaf node and then expand volume around point till the ratio 
    /// of the search volume to the node boundary is greater than one, that is move out from the particle till the search radius is greater than the node enclosing the mid point
    /// once the product of this ratio for all dimensions is > 1, then Enbid kicks out of the loop and the local phase space volume is given by the boundaries of the node enclosing the point.
    /// Another way of doing this is to find the leaf node containing the particle, find its nearest neighbour and use that information to determine the metric spacing.
    inline void KDTree::CalculateMetricSpacing(Int_t target, int treetype, Double_t *smetric){
        CalculateMetricSpacing(bucket[target].GetPosition(),bucket[target].GetVelocity(),treetype,smetric);
    }
    inline void KDTree::CalculateMetricSpacing(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric){
        Double_t xyz[6];
        Int_t nsmooth=8,ivol;//search for nsmooth most distant nearest neighbour
        Int_t ID[nsmooth];
        Double_t dist2[nsmooth],off[ND];
        Double_t sr0[ND],vsmooth=0.,kvol;//? what is vsmooth ? see below where sr0 is normalized by vsmooth
        for (int i=0;i<3;i++) {xyz[i]=x[i];xyz[i+3]=v[i];}
        //FindNearestPhase(x, v, ID, dist2, nsmooth);
        PriorityQueue *pq=new PriorityQueue(nsmooth); 
        for (Int_t i = 0; i < nsmooth; i++) pq->Push(-1, MAXVALUE);
        for (int i = 0; i < ND; i++) off[i] = 0.0;
        if (period==NULL) root->FindNearestPhase(0.0,bucket,pq,off,x,v);
        else root->FindNearestPhasePeriodic(0.0,bucket,pq,off,period,x,v);
        LoadNN(nsmooth,pq,ID,dist2);
        delete pq;
        //move outward from 2nd most nearest neighbour to nsmooth till all metric distances are not zero (ie product is nolonger zero
        kvol=0.0;
        ivol=1;
        //note that this has an issue since code is setup in such a fashsion that it will fail if ivol==nsmooth and kvol still ==0.
        while (kvol==0.0&&ivol<nsmooth) {
            kvol=1.0;
            for (int j=0;j<ND;j++) {
                sr0[j]=fabs(bucket[ID[ivol]].GetPhase(j)-xyz[j]);
                kvol*=sr0[j];
            }
            ivol++;
        }
        //if cubic cells then have to get for physical and velocity volume separately then take vol^(1/3) and take dist^2
        //otherwise take 1./sr0[j]/sr0[j] as smetric
        for (int j=0;j<ND;j++) smetric[j]=1.0/(sr0[j]*sr0[j]);
    }
    inline void KDTree::CalculateMetricSpacing(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric){
        Double_t xyz[6];
        Int_t nsmooth=8,ivol;//search for nsmooth most distant nearest neighbour
        Int_t ID[nsmooth];
        Double_t dist2[nsmooth],off[ND];
        Double_t sr0[ND],vsmooth=0.,kvol;//? what is vsmooth ? see below where sr0 is normalized by vsmooth
        for (int i=0;i<3;i++) {xyz[i]=x[i];xyz[i+3]=v[i];}
        //FindNearestPhase(x, v, ID, dist2, nsmooth);
        PriorityQueue *pq=new PriorityQueue(nsmooth); 
        for (Int_t i = 0; i < nsmooth; i++) pq->Push(-1, MAXVALUE);
        for (int i = 0; i < ND; i++) off[i] = 0.0;
        if (period==NULL) root->FindNearestPhase(0.0,bucket,pq,off,x,v);
        else root->FindNearestPhasePeriodic(0.0,bucket,pq,off,period,x,v);
        LoadNN(nsmooth,pq,ID,dist2);
        
        /*for (Int_t i = nsmooth-1; i >=0; i--)
        {
            ID[i] = pq->TopQueue();
            if (ID[i] == -1)
            {
                printf("CalculateMetricSpacing failed for unknown reasons\n");
                exit(1);
            }
            dist2[i] = pq->TopPriority();
            pq->Pop();
        }*/
        
        delete pq;
        //move outward from 2nd most nearest neighbour to nsmooth till all metric distances are not zero (ie product is nolonger zero
        kvol=0.0;
        ivol=1;
        //note that this has an issue since code is setup in such a fashsion that it will fail if ivol==nsmooth and kvol still ==0.
        while (kvol==0.0&&ivol<nsmooth) {
            kvol=1.0;
            for (int j=0;j<ND;j++) {
                sr0[j]=fabs(bucket[ID[ivol]].GetPhase(j)-xyz[j]);
                kvol*=sr0[j];
            }
            ivol++;
        }
        //if cubic cells then have to get for physical and velocity volume separately then take vol^(1/3) and take dist^2
        //otherwise take 1./sr0[j]/sr0[j] as smetric
        for (int j=0;j<ND;j++) smetric[j]=1.0/(sr0[j]*sr0[j]);
    }
    /*
    //original Enbid like code
    inline void KDTree::CalculateMetricSpacing(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric){
        Double_t xyz[6];
        for (int i=0;i<3;i++) {xyz[i]=x[i];xyz[i+3]=v[i];}
        //search radius 0, 1, search factor, and search boundary
        Double_t sr0[ND],sr1[ND],sf=1.2,sx[ND][2];
        //leaf node pointer
        LeafNode *lp;
        Double_t temp1,temp2,temp3,temp4,kvol,volcheck=0.,vsmooth=0.;

        lp=(LeafNode*)FindLeafNode(xyz);
        for (int j=0;j<ND;j++) {
            sr1[j]=sf*0.25*(lp->GetBoundary(j,1)-lp->GetBoundary(j,0));
        }
        while (volcheck<1.0){
            for (int j=0;j<ND;j++) {
                sr1[j]*=2.;
                sx[j][0]=xyz[j]-sr1[j];
                sx[j][1]=xyz[j]+sr1[j];
            }
            lp=(LeafNode*)FindLeafNode(sx);
            kvol=1.0;
            for(int j=0; j<ND; j++)
            {
                temp1=sx[j][0];             temp2=sx[j][1];
                temp3=lp->GetBoundary(j,0); temp4=lp->GetBoundary(j,1);
                if(temp3>temp1) temp1=temp3;
                if(temp4<temp2) temp2=temp4;
                if (temp2<=temp1) break;//go back to start of while loop, this really should not happen
                else kvol*=(temp2-temp1)/(temp4-temp3);
            }
            if (temp2>temp1){
                for(int j=0; j<ND; j++) sr0[j]=(lp->GetBoundary(j,1)-lp->GetBoundary(j,0));
                volcheck+=kvol;//see if break out of loop
                vsmooth+=1.0;
                if(volcheck>1.) for(int j=0; j<ND; j++) sr0[j]/=vsmooth;
            }
        }
        //if cubic cells then have to get for physical and velocity volume then take vol^(1/3) and take dist^2
        //otherwise take 1./sr0[j]/sr0[j] as smetric
        for (int j=0;j<ND;j++) smetric[j]=1.0/(sr0[j]*sr0[j]);
    }
    */
    /// Calculate the distribution of mass in the full phase-space. 
    /// to ensure good statistics must use ~10 times the number of nearest neighbours used to get initial metric estimate, nsmooth here is 64 (but 128 is also okay)
    inline void KDTree::CalculateMetricTensor(Int_t target, int treetype, Double_t *smetric, Double_t *metric, GMatrix &gmetric){
        CalculateMetricTensor(bucket[target].GetPosition(), bucket[target].GetVelocity(),treetype,smetric,metric,gmetric);
    }
    inline void KDTree::CalculateMetricTensor(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &gmetric){
        Double_t xyz[6];
        int nsmooth=64;
        PriorityQueue *pq=new PriorityQueue(nsmooth); 
        Int_t qsize=pq->MaxSize();
        Int_t idlist[qsize];
        Double_t r2list[qsize];
        Double_t off[6];
        for (int i=0;i<3;i++) {xyz[i]=x[i];xyz[i+3]=v[i];}
        //use initila metric to find nearest neighbours to then calculate covariance metric
        for (Int_t i = 0; i < qsize; i++) pq->Push(-1, MAXVALUE);
        for (int i = 0; i < ND; i++) off[i] = 0.0;
        if (period==NULL) root->FindNearestMetric(0.0,bucket,pq,off,x,v,smetric);
        else root->FindNearestMetricPeriodic(0.0,bucket,pq,off,period,x,v,smetric);

        Double_t cm[ND],dx[ND],rcm=0.,mtot=0.,eigmax=0.;
        GMatrix mdisp(ND,ND),eigval(ND,1),eigvec(ND,ND);
        for (int i=0;i<ND;i++) {
            cm[i]=0.;off[i]=0.;
            for (int j=0;j<ND;j++) mdisp(i,j)=0.;
        }
        //calculate center of mass based on all nearby particles in the priority queue
        //find the nearest particles and store them in the queue
        for (Int_t i=0;i<qsize;i++){idlist[i]=pq->TopQueue();r2list[i]=pq->TopPriority();pq->Pop();}
        for (Int_t i=0;i<qsize;i++){
            for (int j=0;j<ND;j++)
                cm[j]+=(xyz[j]-bucket[idlist[i]].GetPhase(j))*bucket[idlist[i]].GetMass();
            mtot+=bucket[idlist[i]].GetMass();
        }
        for (int i=0;i<ND;i++) {cm[i]/=mtot;rcm+=cm[i]*cm[i];}

        //then calculate the mass dispersion about xyz-cm
        for (Int_t i=0;i<qsize;i++){
            for (int j=0;j<ND;j++)
                dx[j]=(xyz[j]-cm[j])-bucket[idlist[i]].GetPhase(j);
            for(int j=0;j<ND;j++)
                for(int k=0;k<ND;k++)
                    mdisp(j,k)+=dx[j]*dx[k]*bucket[idlist[i]].GetMass();
        }
        for (int i=0;i<ND;i++) 
            for(int j=0;j<ND;j++)
                mdisp(i,j)*=sqrt(smetric[i])*sqrt(smetric[j])/(mtot*r2list[0]);

        //then get the eigenvalues of the above matrix and the transformation matrix ie the eigenvectors
        mdisp.Eigenvalvec(eigval, eigvec);
        //might not really need to take sqrt of eigenvalues
        //for (int i=0;i<ND;i++)eigval(i,0)=sqrt(eigval(i,0));
        for (int i=0;i<ND;i++) if (eigmax<eigval(i,0))eigmax=eigval(i,0);
        for (int i=0;i<ND;i++) metric[i]=eigmax/eigval(i,0);
        for (int i=0;i<ND;i++)
            for (int j=0;j<ND;j++)
                gmetric(i,j)=eigvec(i,j);
    }
    inline void KDTree::CalculateMetricTensor(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &gmetric){
        Double_t xyz[6];
        int nsmooth=64;
        PriorityQueue *pq=new PriorityQueue(nsmooth); 
        Int_t qsize=pq->MaxSize();
        Int_t idlist[qsize];
        Double_t r2list[qsize];
        Double_t off[6];
        for (int i=0;i<3;i++) {xyz[i]=x[i];xyz[i+3]=v[i];}
        //use initila metric to find nearest neighbours to then calculate covariance metric
        for (Int_t i = 0; i < qsize; i++) pq->Push(-1, MAXVALUE);
        for (int i = 0; i < ND; i++) off[i] = 0.0;
        if (period==NULL) root->FindNearestMetric(0.0,bucket,pq,off,x,v,smetric);
        else root->FindNearestMetricPeriodic(0.0,bucket,pq,off,period,x,v,smetric);

        Double_t cm[ND],dx[ND],rcm=0.,mtot=0.,eigmax=0.;
        GMatrix mdisp(ND,ND),eigval(ND,1),eigvec(ND,ND);
        for (int i=0;i<ND;i++) {
            cm[i]=0.;off[i]=0.;
            for (int j=0;j<ND;j++) mdisp(i,j)=0.;
        }
        //calculate center of mass based on all nearby particles in the priority queue
        //find the nearest particles and store them in the queue
        for (Int_t i=0;i<qsize;i++){idlist[i]=pq->TopQueue();r2list[i]=pq->TopPriority();pq->Pop();}
        for (Int_t i=0;i<qsize;i++){
            for (int j=0;j<ND;j++)
                cm[j]+=(xyz[j]-bucket[idlist[i]].GetPhase(j))*bucket[idlist[i]].GetMass();
            mtot+=bucket[idlist[i]].GetMass();
        }
        for (int i=0;i<ND;i++) {cm[i]/=mtot;rcm+=cm[i]*cm[i];}

        //then calculate the mass dispersion about xyz-cm
        for (Int_t i=0;i<qsize;i++){
            for (int j=0;j<ND;j++)
                dx[j]=(xyz[j]-cm[j])-bucket[idlist[i]].GetPhase(j);
            for(int j=0;j<ND;j++)
                for(int k=0;k<ND;k++)
                    mdisp(j,k)+=dx[j]*dx[k]*bucket[idlist[i]].GetMass();
        }
        for (int i=0;i<ND;i++) 
            for(int j=0;j<ND;j++)
                mdisp(i,j)*=sqrt(smetric[i])*sqrt(smetric[j])/(mtot*r2list[0]);

        //then get the eigenvalues of the above matrix and the transformation matrix ie the eigenvectors
        mdisp.Eigenvalvec(eigval, eigvec);
        //might not really need to take sqrt of eigenvalues
        //for (int i=0;i<ND;i++)eigval(i,0)=sqrt(eigval(i,0));
        for (int i=0;i<ND;i++) if (eigmax<eigval(i,0))eigmax=eigval(i,0);
        for (int i=0;i<ND;i++) metric[i]=eigmax/eigval(i,0);
        for (int i=0;i<ND;i++)
            for (int j=0;j<ND;j++)
                gmetric(i,j)=eigvec(i,j);
//distance measure in Enbid is this where in my case, metric1=sqrt(smetric)
//so must be a way of altering this but I am not certain why this form which is effectively
//eigenvalue(k)* (sum(delta x of coord l)*eigenvector(l,k)*sqrt(smetric(l)))
/*
    for(k=0,r=0; k<ND; k++)
    {
        for(l=0,temp=0; l<ND; ++l)
        {
        temp+=(Part[i].Pos[l]-searchcenter[l])*gmatrix[l][k]*metric1[l];
        }
        temp=temp*metric[k];
        r+=temp*temp;

    }
*/
    }
    //@}

    //------- Begining Find NN routines

    /// Find the m nearest neighbours of the tt'th particle in the tree.
    /// The nearest neighbours and squared distances are stored in
    /// the arrays nn and dist2, respectively.
    void KDTree::FindNearest(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        //here if periodic and searching based on target, increase Nsearch by 1 because coord based search will find target, then remove target from list at the end
        if (period!=NULL) Nsearch+=1;
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[ND];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);
        for (int i = 0; i < ND; i++) off[i] = 0.0;

        if (period==NULL) {
        if (treetype==TPHYS) root->FindNearestPos(0.0,bucket,pq,off,tt,ND);
        else if (treetype==TVEL) root->FindNearestVel(0.0,bucket,pq,off,tt,ND);
        else if (treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,tt,ND);
        else if (treetype==TPHS) {
            //simple phase space search
            if (anisotropic==-1) root->FindNearestPhase(0.0,bucket,pq,off,tt);
            else {
            //get metric describing phase-space volume around target particle
            //first get the metric spacing of the leaf nodes
            Double_t m0[ND],m1[ND];
            GMatrix gm(ND,ND);
            CalculateMetricSpacing(tt, treetype, m0);
            //now if not accounting for anisotropic dispersion, then just use metric
            if (anisotropic==0) root->FindNearestMetric(0.0,bucket,pq,off,tt,m0);
            //otherwise use a metric tensor to calculate anisotropic kernel
            //see Sharma & Steinmetz 2006 or Shapiro et al 1996
            //kernel instead of being W(|x|) is W(|D^(-1/2) E H^(-1) x|)
            //where H is diagonal matrix with Hii=hi (hi is metric spacing)
            //C is convariance matrix about x`, where x`=H^(-1)x, that is C is mdisp
            //E is eigenvalues matrix that diagonalizes C,
            //D is diagonal matrix of eigenvalues normalized by det|D|^(1/dimensions)
            //that is take eigenvalues and normalize by product of all to the power of 1/dimensions
            //one must also adjust
            else if (anisotropic==1) {
                CalculateMetricTensor(tt, treetype, m0, m1, gm);
                //take sqrts of metrics
                for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}

                //finally find the nearest particles based on the metrics m0,m1, and metric tensor gm
                //for (Int_t j = 0; j <Nsearch; j++) pq->Push(-1, MAXVALUE);
                //root->FindNearestPhase(0.0,bucket,pq,off,target);
                root->FindNearestMetricwithTensor(0.0,bucket,pq,off,tt,m0,m1,gm);
            }
            }
        }
        //else if (treetype==TMETRIC) root->FindNearestMetric(0.0,bucket,pq,off,tt,metric);
        }
        else 
        {
        if (treetype==TPHYS) root->FindNearestPosPeriodic(0.0,bucket,pq,off,period,tt,ND);
        else if (treetype==TVEL) root->FindNearestVel(0.0,bucket,pq,off,tt,ND);
        else if (treetype==TPROJ) root->FindNearestPosPeriodic(0.0,bucket,pq,off,period,tt,ND);
        else if (treetype==TPHS) {
            //simple phase space search
            if (anisotropic==-1) root->FindNearestPhasePeriodic(0.0,bucket,pq,off,period,tt);
            else {
            Double_t m0[ND],m1[ND];
            GMatrix gm(ND,ND);
            CalculateMetricSpacing(tt, treetype, m0);
            if (anisotropic==0) root->FindNearestMetricPeriodic(0.0,bucket,pq,off,period,tt,m0);
            else if (anisotropic==1) {
                CalculateMetricTensor(tt, treetype, m0, m1, gm);
                for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                root->FindNearestMetricwithTensorPeriodic(0.0,bucket,pq,off,period,tt,m0,m1,gm);
            }
            }
        }
        }
        if (period!=NULL) {Nsearch-=1;}
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }

    void KDTree::FindNearestPos(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        if (period!=NULL) Nsearch+=1;
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[3];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 3; i++) off[i] = 0.0;

        if (period==NULL) root->FindNearestPos(0.0,bucket,pq,off,tt);
        else root->FindNearestPosPeriodic(0.0,bucket,pq,off,period,tt);
        if (period!=NULL) {pq->Pop();Nsearch-1;}
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    void KDTree::FindNearestVel(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[3];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 3; i++) off[i] = 0.0;

        root->FindNearestVel(0.0,bucket,pq,off,tt);
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    void KDTree::FindNearestPhase(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        if (period!=NULL) Nsearch+=1;
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[6];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 6; i++) off[i] = 0.0;

        if (period==NULL) root->FindNearestPhase(0.0,bucket,pq,off,tt);
        else root->FindNearestPhasePeriodic(0.0,bucket,pq,off,period,tt);
        if (period!=NULL) {Nsearch-=1;}
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    //criterion NN
    void KDTree::FindNearestCriterion(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        if (period!=NULL) Nsearch+=1;
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[3];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 3; i++) off[i] = 0.0;

        if (period==NULL) root->FindNearestCriterion(0.0,cmp,params,bucket,pq,off,tt);
        else root->FindNearestCriterionPeriodic(0.0,cmp,params,bucket,pq,off,period,tt);
        if (period!=NULL) {Nsearch-=1;}
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    //criterion NN
    void KDTree::FindNearestCriterion(Particle p, FOFcompfunc cmp, Double_t *params, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        if (period!=NULL) Nsearch+=1;
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[3];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 3; i++) off[i] = 0.0;

        if (period==NULL) root->FindNearestCriterion(0.0,cmp,params,bucket,pq,off,p);
        else root->FindNearestCriterionPeriodic(0.0,cmp,params,bucket,pq,off,period,p);
        if (period!=NULL) {Nsearch-=1;}
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }

    // Same as above but done for every particle
    void KDTree::FindNearest(Int_t **nn, Double_t **dist2, Int_t Nsearch){
        for (Int_t i=0;i<numparts;i++) FindNearest(i,nn[i],dist2[i],Nsearch);

    }
    void KDTree::FindNearestPos(Int_t **nn, Double_t **dist2, Int_t Nsearch){
        for (Int_t i=0;i<numparts;i++) FindNearestPos(i,nn[i],dist2[i],Nsearch);
    }
    void KDTree::FindNearestVel(Int_t **nn, Double_t **dist2, Int_t Nsearch){
        for (Int_t i=0;i<numparts;i++) FindNearestVel(i,nn[i],dist2[i],Nsearch);
    }
    void KDTree::FindNearestPhase(Int_t **nn, Double_t **dist2, Int_t Nsearch){
        for (Int_t i=0;i<numparts;i++) FindNearestPhase(i,nn[i],dist2[i],Nsearch);
    }
    void KDTree::FindNearestCriterion(FOFcompfunc cmp, Double_t *params, Int_t **nn, Double_t **dist2, Int_t Nsearch){
        for (Int_t i=0;i<numparts;i++) FindNearestCriterion(i,cmp, params, nn[i],dist2[i],Nsearch);
    }

    // Similar to above but find nearest particles to a coordinate
    void KDTree::FindNearest(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[ND];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < ND; i++) off[i] = 0.0;

        if (period==NULL) {
        if (treetype==TPHYS) root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        else if (treetype==TVEL) root->FindNearestVel(0.0,bucket,pq,off,x,ND);
        else if (treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        else if (treetype==TPHS) {
            Double_t *v=new Double_t[3];
            for (int i=0;i<3;i++) v[i]=x[i+3];
            if (anisotropic==-1) root->FindNearestPhase(0.0,bucket,pq,off,x,v);
            else {
                Double_t m0[ND],m1[ND];
                GMatrix gm(ND,ND);
                CalculateMetricSpacing(x,v,treetype, m0);
                if (anisotropic==0) root->FindNearestMetric(0.0,bucket,pq,off,x,v,m0);
                else if (anisotropic==1) {
                    CalculateMetricTensor(x,v, treetype, m0, m1, gm);
                    for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                    root->FindNearestMetricwithTensor(0.0,bucket,pq,off,x,v,m0,m1,gm);
                }
            }
        }
        //else if (treetype==TMETRIC) root->FindNearestMetric(0.0,bucket,pq,off,x,metric);
        }
        else {
        if (treetype==TPHYS) root->FindNearestPosPeriodic(0.0,bucket,pq,off,period,x,ND);
        else if (treetype==TVEL) root->FindNearestVelPeriodic(0.0,bucket,pq,off,period,x,ND);
        else if (treetype==TPROJ) root->FindNearestPosPeriodic(0.0,bucket,pq,off,period,x,ND);
        else if (treetype==TPHS) {
            Double_t *v=new Double_t[3];
            for (int i=0;i<3;i++) v[i]=x[i+3];
            if (anisotropic==-1) root->FindNearestPhasePeriodic(0.0,bucket,pq,off,period,x,v);
            else {
                Double_t m0[ND],m1[ND];
                GMatrix gm(ND,ND);
                CalculateMetricSpacing(x,v,treetype, m0);
                if (anisotropic==0) root->FindNearestMetricPeriodic(0.0,bucket,pq,off,period,x,v,m0);
                else if (anisotropic==1) {
                    CalculateMetricTensor(x,v, treetype, m0, m1, gm);
                    for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                    root->FindNearestMetricwithTensorPeriodic(0.0,bucket,pq,off,period,x,v,m0,m1,gm);
                }
            }
        }
        }            
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    
    void KDTree::FindNearestPos(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[3];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 3; i++) off[i] = 0.0;

        if (period==NULL) root->FindNearestPos(0.0,bucket,pq,off,x);
        else root->FindNearestPosPeriodic(0.0,bucket,pq,off,period,x);
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    void KDTree::FindNearestVel(Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[3];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 3; i++) off[i] = 0.0;

        root->FindNearestVel(0.0,bucket,pq,off,v);
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    void KDTree::FindNearestPhase(Double_t *x, Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch)
    {
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        Double_t off[6];
        for (Int_t i = 0; i < Nsearch; i++) pq->Push(-1, MAXVALUE);

        for (int i = 0; i < 6; i++) off[i] = 0.0;

        if (period==NULL) root->FindNearestPhase(0.0,bucket,pq,off,x,v);
        else root->FindNearestPhasePeriodic(0.0,bucket,pq,off,period,x,v);
        LoadNN(Nsearch,pq,nn,dist2);
        delete pq;
    }
    // Similar to above but find nearest particles to a coordinate
    void KDTree::FindNearestPos(Coordinate x, Int_t *nn, Double_t *dist2, Int_t Nsearch){
        FindNearestPos(x.GetCoord(),nn,dist2,Nsearch);
    }
    void KDTree::FindNearestVel(Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch){
        FindNearestVel(v.GetCoord(),nn,dist2,Nsearch);
    }
    void KDTree::FindNearestPhase(Coordinate x, Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch){
        FindNearestPhase(x.GetCoord(),v.GetCoord(),nn,dist2,Nsearch);
    }

    // Find particles that lie within a distance fdist2 to target particle
    void KDTree::SearchBallPos(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)
    {
        Double_t off[6];
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchBallPos(0.0,fdist2,imark,bucket,nn,dist2,off,tt);
        else root->SearchBallPosPeriodic(0.0,fdist2,imark,bucket,nn,dist2,off,period,tt);
    }
    // Find particles that lie within a distance fdist2 to target position
    void KDTree::SearchBallPos(Double_t *x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)
    {
        Double_t off[6];
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchBallPos(0.0,fdist2,imark,bucket,nn,dist2,off,x);
        else root->SearchBallPosPeriodic(0.0,fdist2,imark,bucket,nn,dist2,off,period,x);
    }

    // Find particles that lie within a distance fdist2 to target position
    void KDTree::SearchBallPos(Coordinate x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)
    {
        SearchBallPos(x.GetCoord(),fdist2,imark,nn,dist2);
    }

    // Find particles that lie within the comparison function
    void KDTree::SearchCriterion(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t imark, Int_t *nn, Double_t *dist2)
    {
        Double_t off[6];
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchCriterion(0.0,cmp,params,imark,bucket,nn,dist2,off,tt);
        else root->SearchCriterionPeriodic(0.0,cmp,params,imark,bucket,nn,dist2,off,period,tt);
    }
    void KDTree::SearchCriterion(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t imark, Int_t *nn)
    {
        Double_t off[6];
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchCriterionNoDist(0.0,cmp,params,imark,bucket,nn,off,tt);
        else root->SearchCriterionNoDistPeriodic(0.0,cmp,params,imark,bucket,nn,off,period,tt);
    }
/*
    Int_t KDTree::SearchCriterion(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t imark, Int_t *nn, Double_t *dist2, Int_t *tagged)
    {
        Double_t off[6];
        Int_t numtagged=0;
        if(nn==NULL) nn=new Int_t[numparts];
        if(dist2==NULL) dist2=new Double_t[numparts];
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchCriterionTagged(0.0,cmp,params,imark,bucket,nn,dist2,numtagged,tagged,off,tt);
        else root->SearchCriterionPeriodicTagged(0.0,cmp,params,imark,bucket,nn,dist2,numtagged,tagged,off,period,tt);
        return numtagged;
    }
*/
    // Find particles that lie within a distance fdist2 to target particle
    Int_t KDTree::SearchBallPosTagged(Int_t tt, Double_t fdist2, Int_t *tagged)
    {
        Int_t nt=0;
        Double_t off[6];
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchBallPosTagged(0.0,fdist2,bucket,tagged,off,tt,nt);
        else root->SearchBallPosPeriodicTagged(0.0,fdist2,bucket,tagged,off,period,tt,nt);
        return nt;
    }
    // Find particles that lie within a distance fdist2 to target position
    Int_t KDTree::SearchBallPosTagged(Double_t *x, Double_t fdist2, Int_t *tagged)
    {
        Int_t nt=0;
        Double_t off[6];
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchBallPosTagged(0.0,fdist2,bucket,tagged,off,x,nt);
        else root->SearchBallPosPeriodicTagged(0.0,fdist2,bucket,tagged,off,period,x,nt);
        return nt;
    }
    // Find particles that lie within a distance fdist2 to target position
    Int_t KDTree::SearchBallPosTagged(Coordinate x, Double_t fdist2, Int_t *tagged)
    {
        Int_t nt=SearchBallPosTagged(x.GetCoord(),fdist2,tagged);
        return nt;
    }
    Int_t KDTree::SearchCriterionTagged(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t *tagged)
    {
        Double_t off[6];
        Int_t numtagged=0;
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchCriterionTagged(0.0,cmp,params,bucket,numtagged,tagged,off,tt);
        else root->SearchCriterionPeriodicTagged(0.0,cmp,params,bucket,numtagged,tagged,off,period,tt);
        return numtagged;
    }
    Int_t KDTree::SearchCriterionTagged(Particle &p, FOFcompfunc cmp, Double_t *params, Int_t *tagged)
    {
        Double_t off[6];
        Int_t numtagged=0;
        for (int i = 0; i < 3; i++) off[i] = 0.0;
        if (period==NULL) root->SearchCriterionTagged(0.0,cmp,params,bucket,numtagged,tagged,off,p);
        else root->SearchCriterionPeriodicTagged(0.0,cmp,params,bucket,numtagged,tagged,off,period,p);
        return numtagged;
    }


    // Find leaf node containing Particle
    Node *KDTree::FindLeafNode(Int_t tt){
        Node* np=root;
        SplitNode *sp;
        int k;
        while(np->GetCount()>b){
            sp=(SplitNode *)np;
            k=sp->GetCutDim();
            if (bucket[tt].GetPosition(k)<(sp->GetLeft())->GetBoundary(k,1)) np=sp->GetLeft();
            else np=sp->GetRight();
        }
        return np;
    }
    // Find leaf note containing point
    Node *KDTree::FindLeafNode(Double_t *x){
        Node* np=root;
        SplitNode *sp;
        int k;
        while(np->GetCount()>b){
            sp=(SplitNode *)np;
            k=sp->GetCutDim();
            if (x[k]<(sp->GetLeft())->GetBoundary(k,1)) np=sp->GetLeft();
            else np=sp->GetRight();
        }
        return np;
    }
    // Find Node enclosing volume
    Node *KDTree::FindLeafNode(Double_t search[6][2]){
        Node* np=root;
        SplitNode *sp;
        int k;
        while(np->GetCount()>b){
            sp=(SplitNode *)np;
            k=sp->GetCutDim();
            if (search[k][0]<((SplitNode *)(sp->GetLeft()))->GetCutValue()) np=sp->GetLeft();
            if (search[k][1]>((SplitNode *)(sp->GetRight()))->GetCutValue()) np=sp->GetRight();
        }
        return np;
    }

}
