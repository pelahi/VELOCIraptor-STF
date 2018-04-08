/*! \file KDCalcSmoothQuantities.cxx
 *  \brief This file contains subroutines involving calculating kernel smooth quantities.

*/

#include <KDTree.h>

namespace NBody
{
    /// inline functions
    /// get smoothing kernel
    inline Double_t KDTree::Wsm(Double_t r, Int_t i, Int_t size, Double_t delta, Double_t*x){
        if (i<size-1) return (x[i]+(x[i+1]-x[i])*(r-delta*i)/delta);
        else return x[i];//(x[i]-(x[i]-x[i-1])*(delta*(i+1.)-r)/delta);
    }
    /// get metric describing local phase-space volume (see KDFindNearest for description) (since inline, copied again, think of better way of including the functions)
    inline void KDTree::CalculateMetricSpacing(Int_t target, int treetype, Double_t *smetric){
        CalculateMetricSpacing(bucket[target].GetPosition(),bucket[target].GetVelocity(),treetype,smetric);
    }
    inline void KDTree::CalculateMetricSpacing(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric){
        Double_t xyz[6];
        Int_t nsmooth=8,ivol;//search for nsmooth most distant nearest neighbour
        Int_t ID[nsmooth];
        Double_t off[ND];
        Double_t sr0[ND],kvol;
        for (int i=0;i<3;i++) {xyz[i]=x[i];xyz[i+3]=v[i];}
        PriorityQueue *pq=new PriorityQueue(nsmooth);
        for (Int_t i = 0; i < nsmooth; i++) pq->Push(-1, MAXVALUE);
        for (int i = 0; i < ND; i++) off[i] = 0.0;
        root->FindNearestPhase(0.0,bucket,pq,off,x,v);
        for (Int_t i = nsmooth-1; i >=0; i--)
        {
            ID[i] = pq->TopQueue();
            if (ID[i] == -1)
            {
                printf("CalculateMetricSpacing failed for unknown reasons\n");
                exit(1);
            }
            pq->Pop();
        }
        delete pq;
        //move outward from 2nd most nearest neighbour to nsmooth till all metric distances are not zero (ie product is nolonger zero
        kvol=0.0;
        ivol=1;
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
        Double_t off[ND];
        Double_t sr0[ND],kvol;
        for (int i=0;i<3;i++) {xyz[i]=x[i];xyz[i+3]=v[i];}
        PriorityQueue *pq=new PriorityQueue(nsmooth);
        for (Int_t i = 0; i < nsmooth; i++) pq->Push(-1, MAXVALUE);
        for (int i = 0; i < ND; i++) off[i] = 0.0;
        root->FindNearestPhase(0.0,bucket,pq,off,x,v);
        for (Int_t i = nsmooth-1; i >=0; i--)
        {
            ID[i] = pq->TopQueue();
            if (ID[i] == -1)
            {
                printf("CalculateMetricSpacing failed for unknown reasons\n");
                exit(1);
            }
            pq->Pop();
        }
        delete pq;
        //move outward from 2nd most nearest neighbour to nsmooth till all metric distances are not zero (ie product is nolonger zero
        kvol=0.0;
        ivol=1;
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
        root->FindNearestMetric(0.0,bucket,pq,off,x,v,smetric);

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
        root->FindNearestMetric(0.0,bucket,pq,off,x,v,smetric);

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

    /// Calculate the density of each particle in the tree by smoothing over the nearest Nsmooth particles.
    void KDTree::CalcDensity(Int_t Nsmooth)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        for (Int_t i = 0; i < numparts; i++) bucket[i].SetDensity(0);
        //create a priority queue
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE;
        Double_t m0[ND],m1[ND],off[ND];
        GMatrix gm(ND,ND);

        //for each particle, first fill queue with infinite distances
        //then search the tree starting at the root node for nearest neighbours
        //this is done by looking for the nodes that are closer than a distance furthest
        //which could be altered to say new_furthest=1000*most distant neighbour of previous
        //particle to reduce search
        for (Int_t i = 0; i < numparts; i++)
        {
            Int_t target = i;

            for (Int_t j = 0; j <Nsmooth; j++) pq->Push(-1, furthest);
            for (int j = 0; j < ND; j++) off[j] = 0.0;
            if (treetype==TPHYS||treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,target,ND);
            else if (treetype==TVEL)root->FindNearestVel(0.0,bucket,pq,off,target,ND);
            else if (treetype==TPHS){
                //simple phase space search
                if (anisotropic==-1) root->FindNearestPhase(0.0,bucket,pq,off,target);
                else {
                    //first get the metric spacing of the local phase-space volume
                    CalculateMetricSpacing(target, treetype, m0);
                    //now if not accounting for anisotropic dispersion, then just use metric
                    if (anisotropic==0) {
                        root->FindNearestMetric(0.0,bucket,pq,off,target,m0);
                        for (int j=0;j<ND;j++)m0[j]=sqrt(m0[j]);
                    }
                    //otherwise use a metric tensor to calculate anisotropic kernel
                    //see Sharma & Steinmetz 2006 or Shapiro et al 1996
                    //kernel instead of being W(|x|) is W(|D^(-1/2) E H^(-1) x|)
                    //where H is diagonal matrix with Hii=hi (hi is metric spacing)
                    //C is convariance matrix about x`, where x`=H^(-1)x, that is C is mdisp
                    //E is eigenvalues matrix that diagonalizes C,
                    //D is diagonal matrix of eigenvalues normalized by det|D|^(1/dimensions)
                    //that is take eigenvalues and normalize by product of all to the power of 1/dimensions
                    //one must also adjust
                    else if (anisotropic==1){
                        //then use that metric spacing to examine the dispersion about the particle
                        CalculateMetricTensor(target, treetype, m0, m1, gm);
                        //take sqrts of metrics
                        for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                        //finally find the nearest particles based on the metrics m0,m1, and metric tensor gm
                        //root->FindNearestPhase(0.0,bucket,pq,off,target);
                        root->FindNearestMetricwithTensor(0.0,bucket,pq,off,target,m0,m1,gm);
                    }
                }
            }
            Double_t hi = 0.5 * sqrt(pq->TopPriority());
			//Normalizing by most distant neighbour
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            //for phase with anisotropic kernel must account for the fact that metric used in finding
            //near neighbours is not unity. see Sharma & Steinmetz for details.
            if (treetype==TPHS) {
                Double_t temp=1.0;
                if (anisotropic==0) for (int k=0;k<ND;k++)temp*=m0[k];
                else for (int k=0;k<ND;k++)temp*=m0[k]*m1[k];
                norm*=temp;
            }
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcDensity failed for some reason\n");
                    exit(1);
                }
                Double_t rij = sqrt(pq->TopPriority());
                //Double_t Wij = 0.5 * WSPH(rij, hi);
                //smoothing kernel used to get weight of particle in SPH calculation
                Double_t Wij = 0.5 * Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

                Int_t id = pq->TopQueue();

                bucket[i].SetDensity(bucket[i].GetDensity() + Wij * bucket[id].GetMass());
                bucket[id].SetDensity(bucket[id].GetDensity() + Wij * bucket[i].GetMass());
                //if want to reduce search by setting a maximum distance
                /*if (i < numparts - 1)
                {
                    Double_t d;
                    if (treetype==0) d = DistanceSqd(bucket[i+1].GetPosition(),bucket[pq->TopQueue()].GetPosition());
                    else if (treetype==1) d = VelDistSqd(bucket[i+1].GetVelocity(), bucket[pq->TopQueue()].GetVelocity());
                    else if (treetype==2) d = PhaseDistSqd(bucket[i+1].GetPosition(),bucket[pq->TopQueue()].GetPosition(),
                                        bucket[i+1].GetVelocity(),bucket[pq->TopQueue()].GetVelocity());
                    else if (treetype==3) d = DimDistSqd(bucket[i+1].GetPosition(),bucket[pq->TopQueue()].GetPosition(),
                                        bucket[i+1].GetVelocity(),bucket[pq->TopQueue()].GetVelocity(),metric);
                    if (d > new_furthest) new_furthest = d;
                }*/
                pq->Pop();
            }
            //furthest = new_furthest*1000;
        }
        if (scalespace) for (Int_t i = 0; i < numparts; i++) bucket[i].SetDensity(bucket[i].GetDensity()*ivol);
        delete pq;
    }

    /// Calculate the velocity density of each particle in the tree by smoothing over
    /// the nearest Nsmooth velocity neighbours from set of Nsearch physical neighbours
    void KDTree::CalcVelDensity(Int_t Nsmooth, Int_t Nsearch)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (!(treetype==TPHYS||treetype==TPROJ||treetype==TPHS)) {
            printf("CalcVelDensity is only relvant if tree is physical or phase-space tree\n");
            exit(1);
        }
        if (Nsmooth>Nsearch) {
            printf("CalcVelDensity Nsmooth must be < Nsearch, setting Nsmooth=Nsearch\n");
            Nsmooth=Nsearch;
        }

        for (Int_t i = 0; i < numparts; i++) bucket[i].SetDensity(0);

        //create priority queues
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        PriorityQueue *pq2=new PriorityQueue(Nsmooth);
        Int_t *nnIDs=new Int_t[Nsearch];
        Double_t *vdist=new Double_t[Nsearch];
        Double_t furthest = MAXVALUE,off[ND];

        for (Int_t i = 0; i < numparts; i++)
        {
            for (Int_t j = 0; j <Nsearch; j++) pq->Push(-1, furthest);
            for (int j = 0; j < ND; j++) off[j] = 0.0;
            //find nearest physical neighbours or phase-space neighbours
            if (treetype==TPHYS||treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,i,ND);
            else if (treetype==TPHS) {
                if (anisotropic==-1) root->FindNearestPhase(0.0,bucket,pq,off,i);
                else {
                    Double_t m0[ND],m1[ND];
                    GMatrix gm(ND,ND);
                    CalculateMetricSpacing(i, treetype, m0);
                    if (anisotropic==0) root->FindNearestMetric(0.0,bucket,pq,off,i,m0);
                    else if (anisotropic==1) {
                        CalculateMetricTensor(i, treetype, m0, m1, gm);
                        //take sqrts of metrics
                        for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                        root->FindNearestMetricwithTensor(0.0,bucket,pq,off,i,m0,m1,gm);
                    }
                }
            }
            for (Int_t j = 0; j <Nsearch; j++) {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcVelDensity failed for some reason\n");
                    exit(1);
                }
                nnIDs[j]=pq->TopQueue();
                vdist[j]=sqrt(VelDistSqd(bucket[i].GetVelocity(),bucket[nnIDs[j]].GetVelocity(),ND));
                pq->Pop();
            }
            for (Int_t j = 0; j <Nsmooth; j++) pq2->Push(-1, furthest);
            //from this set find nearest velocity neighbours
            for (Int_t j=0;j<Nsearch;j++)
                if (vdist[j] < pq2->TopPriority()){
                    pq2->Pop();
                    pq2->Push(nnIDs[j], vdist[j]);
                }
            //and calculate velocity density from subset of physical near neighbours using nearest velocity neighbours
            Double_t hi=0.5*pq2->TopPriority();
            //Normalizing by most distant velocity neighbour
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                Double_t rij = pq2->TopPriority();
                //Int_t id=pq2->TopQueue();
                //smoothing kernel used to get weight of particle in SPH calculation
                Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
                bucket[i].SetDensity(bucket[i].GetDensity() + Wij);
                pq2->Pop();
            }
        }
        delete pq;
        delete pq2;
        delete[] vdist;
        delete[] nnIDs;
    }

    /// Calculate the velocity density of each particle in the tree by smoothing over the nearest Nsmooth velocity neighbours from set of Nsearch physical neighbours
    /// where the velocity density is weighted by the overlap of the particles in physical space given by W(rij)/rho(i)*M(i)
    Double_t *KDTree::CalcVelDensityWithPhysDensity(Int_t Nsmooth, Int_t Nsearch, int densityset)
    {
        Double_t *vden;
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (!(treetype==TPHYS||treetype==TPROJ)) {
            printf("CalcVelDensity is only relvant if tree is physical tree\n");
            exit(1);
        }
        if (Nsmooth>Nsearch) {
            printf("CalcVelDensity Nsmooth must be < Nsearch, setting Nsmooth=Nsearch\n");
            Nsmooth=Nsearch;
        }
        if (densityset!=1) CalcDensity(Nsmooth);

        //create priority queues
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        //would like to use NPriorityQueue but not working yet
        //NPriorityQueue *pq2=new NPriorityQueue(Nsmooth,1);
        PriorityQueue *pq2=new PriorityQueue(Nsmooth);
        Int_t *nnIDs=new Int_t[Nsearch];
        //Double_t *xdist=new Double_t[Nsearch];
        Double_t *vdist=new Double_t[Nsearch];
        Double_t furthest = MAXVALUE,off[ND];
		vden=new Double_t[numparts];
		for (Int_t i=0;i<numparts;i++) vden[i]=0;

        for (Int_t i = 0; i < numparts; i++)
        {
            for (Int_t j = 0; j <Nsearch; j++) pq->Push(-1, furthest);
            for (int j = 0; j < ND; j++) off[j] = 0.0;
            //find nearest physical neighbours
            root->FindNearestPos(0.0,bucket,pq,off,i,ND);
            Double_t hi=0.5*pq->TopPriority();
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            for (Int_t j = 0; j <Nsearch; j++) {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcVelDensitywithPhysDensity failed for some reason\n");
                    exit(1);
                }
                nnIDs[j]=pq->TopQueue();
                //xdist[j]=pq->TopPriority();
                vdist[j]=sqrt(VelDistSqd(bucket[i].GetVelocity(),bucket[nnIDs[j]].GetVelocity(),ND));
                pq->Pop();
            }
            //for (Int_t j = 0; j <Nsmooth; j++) pq2->Push(-1, furthest,0);
            for (Int_t j = 0; j <Nsmooth; j++) pq2->Push(-1, furthest);
            //from this set find nearest velocity neighbours
            for (Int_t j=0;j<Nsearch;j++)
                if (vdist[j] < pq2->TopPriority()){
                    pq2->Pop();
                    //pq2->Push(nnIDs[j], vdist[j],&xdist[j]);
                    pq2->Push(nnIDs[j], vdist[j]);
                }
            //and calculate velocity density from subset of physical near neighbours using nearest velocity neighbours
            Double_t vhi=0.5*pq2->TopPriority();
			//Normalizing by most distant velocity neighbour
            Double_t vnorm=1.0/pow(vhi,(Double_t)(ND*1.));
            Int_t id = bucket[i].GetID();
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                Double_t vij = pq2->TopPriority();
                Int_t id2=pq2->TopQueue();
                Double_t rij=sqrt(DistanceSqd(bucket[i].GetPosition(),bucket[id2].GetPosition(),ND));
                //Double_t rij=pq2->TopValue(0);
                //smoothing kernel used to get weight of particle in SPH calculation
                Double_t vWij = Wsm(vij/vhi, (int)(vij/vhi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*vnorm;
                Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
                //here weight velocity distance by the physical overlap between the particles
                vden[id]+=vWij*(Wij/bucket[i].GetDensity() * bucket[i].GetMass());
                pq2->Pop();
            }
        }
        delete pq;
        delete pq2;
        //delete[] xdist;
        delete[] vdist;
        delete[] nnIDs;

        return vden;
    }

    /// Calculates mean smoothed velocity
    /// Requires that physical density be calculated, ie: treetype==0
    Coordinate * KDTree::CalcSmoothVel(Int_t Nsmooth, int densityset)
    {
        Coordinate *smvel;

		if (root==NULL) {
			printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
		}
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVel is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
		smvel=new Coordinate[numparts];

		for (Int_t i=0;i<numparts;i++)
            for (int j=0;j<ND;j++) smvel[i][j]=0;

        for (Int_t i = 0; i < numparts; i++)
        {
            Int_t target = i;

            for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
            for (int j = 0; j < ND; j++) off[j] = 0.0;
            root->FindNearestPos(0.0,bucket,pq,off,target,ND);

            Double_t hi = 0.5 * sqrt(pq->TopPriority());
			//Normalizing by most distant neighbour
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcSmoothVel failed for some reason\n");
                    exit(1);
                }
                Double_t rij = sqrt(pq->TopPriority());
                Double_t Wij = 0.5 * Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

                Int_t index = pq->TopQueue();
                Int_t id1 = bucket[i].GetID();
                Int_t id2 = bucket[index].GetID();
                Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
                for (int k=0;k<ND;k++)
                    smvel[id1][k]+=temp* bucket[index].GetVelocity(k);
                temp=Wij/bucket[i].GetDensity() * bucket[i].GetMass();
                for (int k=0;k<ND;k++)
                    smvel[id2][k]+=temp * bucket[i].GetVelocity(k);

                pq->Pop();
            }
        }
        if (scalespace) for (Int_t i = 0; i < numparts; i++) for (int j=0;j<ND;j++) smvel[i][j]*=ivol;
        delete pq;
        return smvel;
    }

    /// Calculates velocity dispersion
    /// like CalcSmoothVel but also requires mean velocities be calculated.
    Matrix * KDTree::CalcSmoothVelDisp(Coordinate *smvel, Int_t Nsmooth, int densityset, int meanvelset)
    {
        Matrix *smveldisp;

		if (root==NULL) {
			printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
		}
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVelDisp is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);
        if (meanvelset!=1){
			if (smvel!=NULL) delete[] smvel;
			smvel=CalcSmoothVel(Nsmooth);
		}
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        smveldisp=new Matrix[numparts];

        for (Int_t i=0;i<numparts;i++)
            for (int j=0;j<ND;j++)for (int k=0;k<ND;k++) smveldisp[i](j,k)=0;

        for (Int_t i=0;i<numparts;i++)
        {
            Int_t target = i;

            for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
            for (int j = 0; j < 6; j++) off[j] = 0.0;
            root->FindNearestPos(0.0,bucket,pq,off,target,ND);
            Double_t dv1,dv2;

            Double_t hi = 0.5 * sqrt(pq->TopPriority());
			//Normalizing by most distant neighbour
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcSmoothVelDisp failed for some reason\n");
                    exit(1);
                }
                Double_t rij = sqrt(pq->TopPriority());
                Double_t Wij = 0.5 * Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

                Int_t index = pq->TopQueue();
                Int_t id1 = bucket[i].GetID();
                Int_t id2 = bucket[index].GetID();

                Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
                for (int k=0;k<ND;k++)
                for (int l=0;l<ND;l++)
                {
                    dv1=bucket[index].GetVelocity(k)-smvel[id1][k];
                    dv2=bucket[index].GetVelocity(l)-smvel[id1][l];
                    smveldisp[id1](k,l)+= temp*dv1*dv2;
                }
                temp=Wij/bucket[i].GetDensity() * bucket[i].GetMass();
                for (int k=0;k<ND;k++)
                for (int l=0;l<ND;l++)
                {
                    dv1=bucket[i].GetVelocity(k)-smvel[id2][k];
                    dv2=bucket[i].GetVelocity(l)-smvel[id2][l];
                    smveldisp[id2](k,l)+=temp*dv1*dv2;
                }
                pq->Pop();
            }
        }
        if (scalespace) for (Int_t i = 0; i < numparts; i++) for (int j=0;j<ND;j++) for (int k=0;k<ND;k++) smveldisp[i](j,k)*=ivol;
        delete pq;
        return smveldisp;
    }

    /// Calculates velocity skewness (here just getting it as function of velocity
    Coordinate * KDTree::CalcSmoothVelSkew(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth, int densityset, int meanvelset, int veldispset)
    {
        Coordinate *smvelskew;

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVelSkew is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);
        if (meanvelset!=1){
            if (smvel!=NULL) delete[] smvel;
            smvel=CalcSmoothVel(Nsmooth);
        }
        if (veldispset!=1){
            if (smveldisp!=NULL) delete[] smveldisp;
            smveldisp=CalcSmoothVelDisp(smvel,Nsmooth);
        }
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        smvelskew=new Coordinate[numparts];

        for (Int_t i=0;i<numparts;i++)
            for (int j=0;j<ND;j++) smvelskew[i][j]=0.;

        for (Int_t i=0;i<numparts;i++)
        {
            Int_t target = i;

            for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
            for (int j = 0; j < ND; j++) off[j] = 0.0;
            root->FindNearestPos(0.0,bucket,pq,off,target,ND);
            Double_t dv1;

            Double_t hi = 0.5 * sqrt(pq->TopPriority());
			//Normalizing by most distant neighbour
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcSmoothVelDisp failed for some reason\n");
                    exit(1);
                }
                Double_t rij = sqrt(pq->TopPriority());
                Double_t Wij = 0.5 * Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

                Int_t index = pq->TopQueue();
                Int_t id1 = bucket[i].GetID();
                Int_t id2 = bucket[index].GetID();

                Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
                for (int k=0;k<ND;k++)
                {
                    dv1=bucket[index].GetVelocity(k)-smvel[id1][k];
                    smvelskew[id1][k]+= temp*dv1*dv1*dv1/pow(smveldisp[id1](k,k),(Double_t)1.5);
                }
                temp=Wij/bucket[i].GetDensity() * bucket[i].GetMass();
                for (int k=0;k<ND;k++)
                {
                    dv1=bucket[i].GetVelocity(k)-smvel[id2][k];
                    smvelskew[id2][k]+=temp*dv1*dv1*dv1/pow(smveldisp[id2](k,k),(Double_t)1.5);
                }
                pq->Pop();
            }
        }
        if (scalespace) for (Int_t i = 0; i < numparts; i++) for (int j=0;j<ND;j++) smvelskew[i][j]*=ivol;
        delete pq;
        return smvelskew;
    }

    /// Calculates velocity kurtosis (here just getting it as function of velocity)
    Coordinate * KDTree::CalcSmoothVelKurtosis(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth, int densityset, int meanvelset, int veldispset)
    {
        Coordinate *smvelkurt;

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVelKurtosis is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);
        if (meanvelset!=1){
            if (smvel!=NULL) delete[] smvel;
            smvel=CalcSmoothVel(Nsmooth);
        }
        if (veldispset!=1){
            if (smveldisp!=NULL) delete[] smveldisp;
            smveldisp=CalcSmoothVelDisp(smvel,Nsmooth);
        }
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        smvelkurt=new Coordinate[numparts];

        for (Int_t i=0;i<numparts;i++)
            for (int j=0;j<ND;j++) smvelkurt[i][j]=0.;

        for (Int_t i=0;i<numparts;i++)
        {
            Int_t target = i;

            for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
            for (int j = 0; j < ND; j++) off[j] = 0.0;
            root->FindNearestPos(0.0,bucket,pq,off,target,ND);
            Double_t dv1;

            Double_t hi = 0.5 * sqrt(pq->TopPriority());
			//Normalizing by most distant neighbour
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcSmoothVelDisp failed for some reason\n");
                    exit(1);
                }
                Double_t rij = sqrt(pq->TopPriority());
                Double_t Wij = 0.5 * Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

                Int_t index = pq->TopQueue();
                Int_t id1 = bucket[i].GetID();
                Int_t id2 = bucket[index].GetID();

                Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
                for (int k=0;k<ND;k++)
                {
                    dv1=bucket[index].GetVelocity(k)-smvel[id1][k];
                    smvelkurt[id1][k]+= temp*dv1*dv1*dv1*dv1/(smveldisp[id1](k,k)*smveldisp[id1](k,k))-3.;
                }
                temp=Wij/bucket[i].GetDensity() * bucket[i].GetMass();
                for (int k=0;k<ND;k++)
                {
                    dv1=bucket[i].GetVelocity(k)-smvel[id2][k];
                    smvelkurt[id2][k]+=temp*dv1*dv1*dv1*dv1/(smveldisp[id2](k,k)*smveldisp[id2](k,k))-3.;
                }
                pq->Pop();
            }
        }
        if (scalespace) for (Int_t i = 0; i < numparts; i++) for (int j=0;j<ND;j++) smvelkurt[i][j]*=ivol;
        delete pq;
        return smvelkurt;
    }

    //-- Like the above routines but for single particle

    Double_t KDTree::CalcDensityParticle(Int_t target, Int_t Nsmooth)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE, den=0.;
        Double_t m0[ND],m1[ND],off[ND];
        GMatrix gm(ND,ND);

        for (Int_t j = 0; j <Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        if (treetype==TPHYS||treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,target,ND);
        else if (treetype==TVEL)root->FindNearestVel(0.0,bucket,pq,off,target,ND);
        else if (treetype==TPHS){
            //simple phase space search
            if (anisotropic==-1) root->FindNearestPhase(0.0,bucket,pq,off,target);
            else {
                //first get the metric spacing of the local phase-space volume
                CalculateMetricSpacing(target, treetype, m0);
                if (anisotropic==0) {
                    //now if not accounting for anisotropic dispersion, then just use metric
                    root->FindNearestMetric(0.0,bucket,pq,off,target,m0);
                    for (int j=0;j<ND;j++)m0[j]=sqrt(m0[j]);
                }
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
                    //then use that metric spacing to examine the dispersion about the particle
                    CalculateMetricTensor(target, treetype, m0, m1, gm);
                    //take sqrts of metrics
                    for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                    //finally find the nearest particles based on the metrics m0,m1, and metric tensor gm
                    root->FindNearestMetricwithTensor(0.0,bucket,pq,off,target,m0,m1,gm);
                }
            }
        }
        Double_t hi = 0.5 * sqrt(pq->TopPriority());
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        //for phase with anisotropic kernel must account for the fact that metric used in finding
        //near neighbours is not unity. see Sharma & Steinmetz for details.

        if (treetype==TPHS) {
            Double_t temp=1.0;
            if (anisotropic==0) for (int k=0;k<ND;k++)temp*=m0[k];
            else for (int k=0;k<ND;k++)temp*=m0[k]*m1[k];
            norm*=temp;
        }

        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcDensity failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij =  Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

            Int_t id = pq->TopQueue();

            den+= Wij * bucket[id].GetMass();
            pq->Pop();
        }
        if (scalespace) den*=ivol;
        delete pq;
        return den;
    }
    Double_t KDTree::CalcVelDensityParticle(Int_t target, Int_t Nsmooth, Int_t Nsearch, int iflag, PriorityQueue *pq, PriorityQueue *pq2, Int_t *nnIDs, Double_t *vdist)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (!(treetype==TPHYS||treetype==TPROJ||treetype==TPHS)) {
            printf("CalcVelDensity is only relvant if tree is physical or phase-space tree\n");
            exit(1);
        }
        if (Nsmooth>Nsearch) {
            printf("CalcVelDensity Nsmooth must be < Nsearch, setting Nsmooth=Nsearch\n");
            Nsmooth=Nsearch;
        }

        if (iflag==0) {
        pq=new PriorityQueue(Nsearch);
        pq2=new PriorityQueue(Nsmooth);
        nnIDs=new Int_t[Nsearch];
        vdist=new Double_t[Nsearch];
        }
        Double_t furthest = MAXVALUE,off[ND];
        Double_t vden=0.;

        for (Int_t j = 0; j <Nsearch; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        if (treetype==TPHYS||treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,target,ND);
        else if (treetype==TPHS) {
            if (anisotropic==-1) root->FindNearestPhase(0.0,bucket,pq,off,target);
            else {
                Double_t m0[ND],m1[ND];
                GMatrix gm(ND,ND);
                CalculateMetricSpacing(target, treetype, m0);
                if (anisotropic==0) root->FindNearestMetric(0.0,bucket,pq,off,target,m0);
                else if (anisotropic==1) {
                    CalculateMetricTensor(target, treetype, m0, m1, gm);
                    //take sqrts of metrics
                    for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                    root->FindNearestMetricwithTensor(0.0,bucket,pq,off,target,m0,m1,gm);
                }
            }
        }
        for (Int_t j = 0; j <Nsearch; j++) {
            if (pq->TopQueue() == -1)
            {
                printf("CalcDensity failed for some reason\n");
                exit(1);
            }
            nnIDs[j]=pq->TopQueue();
            vdist[j]=sqrt(VelDistSqd(bucket[target].GetVelocity(),bucket[nnIDs[j]].GetVelocity(),ND));
            pq->Pop();
        }
        for (Int_t j = 0; j <Nsmooth; j++) pq2->Push(-1, furthest);
        for (Int_t j=0;j<Nsearch;j++)
            if (vdist[j] < pq2->TopPriority()){
                pq2->Pop();
                pq2->Push(nnIDs[j], vdist[j]);
            }
        Double_t hi=0.5*pq2->TopPriority();
        //Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            Double_t rij = pq2->TopPriority();
            //smoothing kernel used to get weight of particle in SPH calculation
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            vden+=Wij;
            pq2->Pop();
        }
        if (iflag==0) {
        delete pq;
        delete pq2;
        delete[] vdist;
        delete[] nnIDs;
        }
        return vden;
    }
    Double_t KDTree::CalcVelDensityWithPhysDensityParticle(Int_t target, Int_t Nsmooth, Int_t Nsearch, int densityset)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (!(treetype==TPHYS||treetype==TPROJ)) {
            printf("CalcVelDensity is only relvant if tree is physical tree\n");
            exit(1);
        }
        if (Nsmooth>Nsearch) {
            printf("CalcVelDensity Nsmooth must be < Nsearch, setting Nsmooth=Nsearch\n");
            Nsmooth=Nsearch;
        }
        if (densityset!=1) CalcDensity(Nsmooth);

        //create priority queues
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        PriorityQueue *pq2=new PriorityQueue(Nsmooth);
        //NPriorityQueue *pq2=new NPriorityQueue(Nsmooth,1);
        Int_t *nnIDs=new Int_t[Nsearch];
        Double_t *vdist=new Double_t[Nsearch];
        Double_t furthest = MAXVALUE,off[ND];
        Double_t vden=0;

            for (Int_t j = 0; j <Nsearch; j++) pq->Push(-1, furthest);
            for (int j = 0; j < ND; j++) off[j] = 0.0;
            //find nearest physical neighbours
            root->FindNearestPos(0.0,bucket,pq,off,target,ND);
            Double_t hi=0.5*pq->TopPriority();
            Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
            for (Int_t j = 0; j <Nsearch; j++) {
                if (pq->TopQueue() == -1)
                {
                    printf("CalcVelDensitywithPhysDensity failed for some reason\n");
                    exit(1);
                }
                nnIDs[j]=pq->TopQueue();
                vdist[j]=sqrt(VelDistSqd(bucket[target].GetVelocity(),bucket[nnIDs[j]].GetVelocity(),ND));
                pq->Pop();
            }
            for (Int_t j = 0; j <Nsmooth; j++) pq2->Push(-1, furthest);
            //from this set find nearest velocity neighbours
            for (Int_t j=0;j<Nsearch;j++)
                if (vdist[j] < pq2->TopPriority()){
                    pq2->Pop();
                    pq2->Push(nnIDs[j], vdist[j]);
                }
            //and calculate velocity density from subset of physical near neighbours using nearest velocity neighbours
            Double_t vhi=0.5*pq2->TopPriority();
			//Normalizing by most distant velocity neighbour
            Double_t vnorm=1.0/pow(vhi,(Double_t)(ND*1.));
            for (Int_t j = 0; j < Nsmooth; j++)
            {
                Double_t vij = pq2->TopPriority();
                Int_t id=pq2->TopQueue();
                Double_t rij=sqrt(DistanceSqd(bucket[target].GetPosition(),bucket[id].GetPosition(),ND));
                //smoothing kernel used to get weight of particle in SPH calculation
                Double_t vWij = Wsm(vij/vhi, (int)(vij/vhi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*vnorm;
                Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
                //here weight velocity distance by the physical overlap between the particles
                vden+=vWij*(Wij/bucket[target].GetDensity() * bucket[target].GetMass());
                pq2->Pop();
            }
        delete pq;
        delete pq2;
        delete[] vdist;
        delete[] nnIDs;
        return vden;
    }
    Coordinate KDTree::CalcSmoothVelParticle(Int_t target, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVel is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Coordinate smvel;

        for (int j=0;j<ND;j++) smvel[j]=0;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,target,ND);

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothVel failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();
            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            for (int k=0;k<3;k++)
                smvel[k]+= temp* bucket[index].GetVelocity(k);

            pq->Pop();
        }
        if (scalespace) for (int j=0;j<ND;j++) smvel[j]*=ivol;
        delete pq;
        return smvel;
    }
    Matrix KDTree::CalcSmoothVelDispParticle(Int_t target, Coordinate smvel, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVelDisp is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Matrix smveldisp;

        for (int j=0;j<ND;j++) for (int k=0;k<ND;k++) smveldisp(j,k)=0;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,target,ND);
        Double_t dv1,dv2;

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothVelDisp failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();

            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            for (int k=0;k<ND;k++)
            for (int l=0;l<ND;l++)
            {
                dv1=bucket[index].GetVelocity(k)-smvel[k];
                dv2=bucket[index].GetVelocity(l)-smvel[l];
                smveldisp(k,l)+= temp*dv1*dv2;
            }
            pq->Pop();
        }
        if (scalespace) for (int j=0;j<ND;j++) for (int k=0;k<ND;k++) smveldisp(j,k)*=ivol;
        delete pq;
        return smveldisp;
    }

    //-- Like above but for a specific coordinate

    Double_t KDTree::CalcDensityPosition(Double_t *x, Int_t Nsmooth, Double_t *v)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE, den=0.;
        Double_t m0[ND],m1[ND],off[ND];
        GMatrix gm(ND,ND);

        for (Int_t j = 0; j <Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        if (treetype==TPHYS||treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        else if (treetype==TVEL)root->FindNearestVel(0.0,bucket,pq,off,x,ND);
        else if (treetype==TPHS){
            CalculateMetricSpacing(x,v, treetype, m0);
            if (anisotropic==0) {
                root->FindNearestMetric(0.0,bucket,pq,off,x,v,m0);
                for (int j=0;j<ND;j++)m0[j]=sqrt(m0[j]);
            }
            //otherwise use a metric tensor
            else if (anisotropic==1){
                CalculateMetricTensor(x,v, treetype, m0, m1, gm);
                for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                root->FindNearestMetricwithTensor(0.0,bucket,pq,off,x,v,m0,m1,gm);
            }
        }
        Double_t hi = 0.5 * sqrt(pq->TopPriority());
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        //for phase with anisotropic kernel must account for the fact that metric used in finding
        //near neighbours is not unity. see Sharma & Steinmetz for details.
        if (treetype>=2) {
            Double_t temp=1.0;
            if (anisotropic==0) for (int k=0;k<ND;k++)temp*=m0[k];
            else for (int k=0;k<ND;k++)temp*=m0[k]*m1[k];
            norm*=temp;
        }
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcDensity failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij =  Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

            Int_t id = pq->TopQueue();

            den+= Wij * bucket[id].GetMass();
            pq->Pop();
        }
        if (scalespace) den*=ivol;
        delete pq;
        return den;
    }
    Double_t KDTree::CalcVelDensityPosition(Double_t *x, Double_t *v, Int_t Nsmooth, Int_t Nsearch)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcVelDensity is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (Nsmooth>Nsearch) {
            printf("CalcVelDensity Nsmooth must be < Nsearch, setting Nsmooth=Nsearch\n");
            Nsmooth=Nsearch;
        }

        //create a priority queue
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        PriorityQueue *pq2=new PriorityQueue(Nsmooth);
        Int_t *nnIDs=new Int_t[Nsearch];
        Double_t *vdist=new Double_t[Nsearch];
        Double_t furthest = MAXVALUE,off[ND];
        Double_t vden=0.;

        for (Int_t j = 0; j <Nsearch; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        for (Int_t j = 0; j <Nsearch; j++) {
            if (pq->TopQueue() == -1)
            {
                printf("CalcDensity failed for some reason\n");
                exit(1);
            }
            nnIDs[j]=pq->TopQueue();
            vdist[j]=sqrt(VelDistSqd(v,bucket[nnIDs[j]].GetVelocity(),ND));
            pq->Pop();
        }
        for (Int_t j = 0; j <Nsmooth; j++) pq2->Push(-1, furthest);
        for (Int_t j=0;j<Nsearch;j++)
            if (vdist[j] < pq2->TopPriority()){
                pq2->Pop();
                pq2->Push(nnIDs[j], vdist[j]);
            }
        Double_t hi=0.5*pq2->TopPriority();
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            Double_t rij = pq2->TopPriority();
            //smoothing kernel used to get weight of particle in SPH calculation
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            vden+=Wij;
            pq2->Pop();
        }
        delete pq;
        delete pq2;
        delete[] vdist;
        delete[] nnIDs;
        return vden;
    }
    Coordinate KDTree::CalcSmoothVelPosition(Double_t *x, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVel is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Coordinate smvel;

        for (int j=0;j<ND;j++) smvel[j]=0;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothVel failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();
            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            for (int k=0;k<3;k++)
                smvel[k]+= temp* bucket[index].GetVelocity(k);

            pq->Pop();
        }
        if (scalespace) for (int j=0;j<ND;j++) smvel[j]*=ivol;
        delete pq;
        return smvel;
    }
    Matrix KDTree::CalcSmoothVelDispPosition(Double_t *x, Coordinate smvel, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
		if (!(treetype==TPHYS||treetype==TPROJ)) {
			printf("CalcSmoothVelDisp is only relvant if physical density was calculated and thus requires physical tree\n");
			exit(1);
		}
        if (densityset!=1) CalcDensity(Nsmooth);
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Matrix smveldisp;

        for (int j=0;j<ND;j++) for (int k=0;k<ND;k++) smveldisp(j,k)=0;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        Double_t dv1,dv2;

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothVelDisp failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();

            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            for (int k=0;k<ND;k++)
            for (int l=0;l<ND;l++)
            {
                dv1=bucket[index].GetVelocity(k)-smvel[k];
                dv2=bucket[index].GetVelocity(l)-smvel[l];
                smveldisp(k,l)+= temp*dv1*dv2;
            }
            pq->Pop();
        }
        if (scalespace) for (int j=0;j<ND;j++) for (int k=0;k<ND;k++) smveldisp(j,k)*=ivol;
        delete pq;
        return smveldisp;
    }
///
    Double_t KDTree::CalcDensityPosition(Real_t *x, Int_t Nsmooth, Real_t *v)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE, den=0.;
        Double_t m0[ND],m1[ND],off[ND];
        GMatrix gm(ND,ND);

        for (Int_t j = 0; j <Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        if (treetype==TPHYS||treetype==TPROJ) root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        else if (treetype==TVEL)root->FindNearestVel(0.0,bucket,pq,off,x,ND);
        else if (treetype==TPHS){
            CalculateMetricSpacing(x,v, treetype, m0);
            if (anisotropic==0) {
                root->FindNearestMetric(0.0,bucket,pq,off,x,v,m0);
                for (int j=0;j<ND;j++)m0[j]=sqrt(m0[j]);
            }
            //otherwise use a metric tensor
            else if (anisotropic==1){
                CalculateMetricTensor(x,v, treetype, m0, m1, gm);
                for (int j=0;j<ND;j++){m0[j]=sqrt(m0[j]);m1[j]=sqrt(m1[j]);}
                root->FindNearestMetricwithTensor(0.0,bucket,pq,off,x,v,m0,m1,gm);
            }
        }
        Double_t hi = 0.5 * sqrt(pq->TopPriority());
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        //for phase with anisotropic kernel must account for the fact that metric used in finding
        //near neighbours is not unity. see Sharma & Steinmetz for details.
        if (treetype>=2) {
            Double_t temp=1.0;
            if (anisotropic==0) for (int k=0;k<ND;k++)temp*=m0[k];
            else for (int k=0;k<ND;k++)temp*=m0[k]*m1[k];
            norm*=temp;
        }
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcDensity failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij =  Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;

            Int_t id = pq->TopQueue();

            den+= Wij * bucket[id].GetMass();
            pq->Pop();
        }
        if (scalespace) den*=ivol;
        delete pq;
        return den;
    }
    Double_t KDTree::CalcVelDensityPosition(Real_t *x, Real_t *v, Int_t Nsmooth, Int_t Nsearch)
    {
        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (!(treetype==TPHYS||treetype==TPROJ)) {
            printf("CalcVelDensity is only relvant if physical density was calculated and thus requires physical tree\n");
            exit(1);
        }
        if (Nsmooth>Nsearch) {
            printf("CalcVelDensity Nsmooth must be < Nsearch, setting Nsmooth=Nsearch\n");
            Nsmooth=Nsearch;
        }

        //create a priority queue
        PriorityQueue *pq=new PriorityQueue(Nsearch);
        PriorityQueue *pq2=new PriorityQueue(Nsmooth);
        Int_t *nnIDs=new Int_t[Nsearch];
        Double_t *vdist=new Double_t[Nsearch];
        Double_t furthest = MAXVALUE,off[ND];
        Double_t vden=0.;

        for (Int_t j = 0; j <Nsearch; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        for (Int_t j = 0; j <Nsearch; j++) {
            if (pq->TopQueue() == -1)
            {
                printf("CalcDensity failed for some reason\n");
                exit(1);
            }
            nnIDs[j]=pq->TopQueue();
            vdist[j]=sqrt(VelDistSqd(v,bucket[nnIDs[j]].GetVelocity(),ND));
            pq->Pop();
        }
        for (Int_t j = 0; j <Nsmooth; j++) pq2->Push(-1, furthest);
        for (Int_t j=0;j<Nsearch;j++)
            if (vdist[j] < pq2->TopPriority()){
                pq2->Pop();
                pq2->Push(nnIDs[j], vdist[j]);
            }
        Double_t hi=0.5*pq2->TopPriority();
        //Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            Double_t rij = pq2->TopPriority();
            //smoothing kernel used to get weight of particle in SPH calculation
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            vden+=Wij;
            pq2->Pop();
        }
        delete pq;
        delete pq2;
        delete[] vdist;
        delete[] nnIDs;
        return vden;
    }
    Coordinate KDTree::CalcSmoothVelPosition(Real_t *x, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (!(treetype==TPHYS||treetype==TPROJ)) {
            printf("CalcSmoothVel is only relvant if physical density was calculated and thus requires physical tree\n");
            exit(1);
        }
        if (densityset!=1) CalcDensity(Nsmooth);

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Coordinate smvel;

        for (int j=0;j<ND;j++) smvel[j]=0;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
        //Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothVel failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();
            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            for (int k=0;k<3;k++)
                smvel[k]+= temp* bucket[index].GetVelocity(k);

            pq->Pop();
        }
        if (scalespace) for (int j=0;j<ND;j++) smvel[j]*=ivol;
        delete pq;
        return smvel;
    }
    Matrix KDTree::CalcSmoothVelDispPosition(Real_t *x, Coordinate smvel, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (!(treetype==TPHYS||treetype==TPROJ)) {
            printf("CalcSmoothVelDisp is only relvant if physical density was calculated and thus requires physical tree\n");
            exit(1);
        }
        if (densityset!=1) CalcDensity(Nsmooth);
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Matrix smveldisp;

        for (int j=0;j<ND;j++) for (int k=0;k<ND;k++) smveldisp(j,k)=0;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);
        Double_t dv1,dv2;

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
        //Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothVelDisp failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();

            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            for (int k=0;k<ND;k++)
            for (int l=0;l<ND;l++)
            {
                dv1=bucket[index].GetVelocity(k)-smvel[k];
                dv2=bucket[index].GetVelocity(l)-smvel[l];
                smveldisp(k,l)+= temp*dv1*dv2;
            }
            pq->Pop();
        }
        if (scalespace) for (int j=0;j<ND;j++) for (int k=0;k<ND;k++) smveldisp(j,k)*=ivol;
        delete pq;
        return smveldisp;
    }
    Double_t KDTree::CalcDensityPosition(Coordinate x, Int_t Nsmooth, Coordinate v)
    {
        return CalcDensityPosition(x.GetCoord(),Nsmooth,v.GetCoord());
    }
    Double_t KDTree::CalcVelDensityPosition(Coordinate x, Int_t Nsmooth, Int_t Nsearch)
    {
        return CalcVelDensityPosition(x.GetCoord(),Nsmooth,Nsearch);
    }
    Coordinate KDTree::CalcSmoothVelPosition(Coordinate x, Int_t Nsmooth, int densityset)
    {
        return CalcSmoothVelPosition(x.GetCoord(),Nsmooth,densityset);
    }
    Matrix KDTree::CalcSmoothVelDispPosition(Coordinate x, Coordinate smvel, Int_t Nsmooth, int densityset)
    {
        return CalcSmoothVelDispPosition(x.GetCoord(), smvel, Nsmooth,densityset);
    }

    ///Calculate some arbitrary local quantity using sph kernel
    Double_t KDTree::CalcSmoothLocalMeanParticle(Int_t target, Double_t *weight, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (densityset!=1) CalcDensity(Nsmooth);

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Double_t mean=0.;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,target,ND);

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothLocalMean failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();
            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            mean+= temp* weight[bucket[index].GetID()];
            pq->Pop();
        }
        if (scalespace) mean*=ivol;
        delete pq;
        return mean;
    }

    ///Calculate the dispersion of some arbitrary local quantity using sph kernel
    Double_t KDTree::CalcSmoothLocalDispParticle(Int_t target, Double_t *weight, Double_t mean, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (densityset!=1) CalcDensity(Nsmooth);
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Double_t disp=0.,delta;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,target,ND);

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothLocalDisp failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();
            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            delta=(weight[bucket[index].GetID()]-mean);
            disp+=temp*delta*delta;
            pq->Pop();
        }
        if (scalespace) disp*=ivol;
        delete pq;
        return disp;
    }

    //-- Like above but calculate it at a position
    Double_t KDTree::CalcSmoothLocalMeanPosition(Double_t *x, Double_t *weight, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (densityset!=1) CalcDensity(Nsmooth);

        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Double_t mean=0.;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothLocalMean failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();
            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            mean+= temp* weight[bucket[index].GetID()];
            pq->Pop();
        }
        if (scalespace) mean*=ivol;
        delete pq;
        return mean;
    }
    Double_t KDTree::CalcSmoothLocalDispPosition(Double_t *x, Double_t *weight, Double_t mean, Int_t Nsmooth, int densityset)
    {

        if (root==NULL) {
            printf("Error in tree construction, rootNode==NULL. Nothing Done.\n");
            exit(1);
        }
        if (densityset!=1) CalcDensity(Nsmooth);
        PriorityQueue *pq=new PriorityQueue(Nsmooth);
        Double_t furthest = MAXVALUE,off[ND];
        Double_t disp=0.,delta;

        for (Int_t j = 0; j < Nsmooth; j++) pq->Push(-1, furthest);
        for (int j = 0; j < ND; j++) off[j] = 0.0;
        root->FindNearestPos(0.0,bucket,pq,off,x,ND);

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
    	//Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            if (pq->TopQueue() == -1)
            {
                printf("CalcSmoothLocalDisp failed for some reason\n");
                exit(1);
            }
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            Int_t index = pq->TopQueue();
            Double_t temp=Wij/bucket[index].GetDensity() * bucket[index].GetMass();
            delta=(weight[bucket[index].GetID()]-mean);
            disp+=temp*delta*delta;
            pq->Pop();
        }
        if (scalespace) disp*=ivol;
        delete pq;
        return disp;
    }
    Double_t KDTree::CalcSmoothLocalMeanPosition(Coordinate x, Double_t *weight, Int_t Nsmooth, int densityset)
    {
        return CalcSmoothLocalMeanPosition(x.GetCoord(),weight,Nsmooth,densityset);
    }
    Double_t KDTree::CalcSmoothLocalDispPosition(Coordinate x, Double_t *weight, Double_t localmean, Int_t Nsmooth, int densityset)
    {
        return CalcSmoothLocalDispPosition(x.GetCoord(),weight,localmean,Nsmooth,densityset);
    }

    Double_t KDTree::CalcSmoothLocalValue(Int_t Nsmooth, PriorityQueue *pq, Double_t* weight)
    {

        Double_t hi = 0.5 * sqrt(pq->TopPriority());
        //Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        Double_t value=0;
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            Double_t rij = sqrt(pq->TopPriority());
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            value+= Wij* weight[j];
            pq->Pop();
        }
        return value;
    }

    Double_t KDTree::CalcSmoothLocalValue(Int_t Nsmooth, Double_t *dist, Double_t* weight)
    {

        Double_t hi = 0.5 * dist[0];
        //Normalizing by most distant neighbour
        Double_t norm=1.0/pow(hi,(Double_t)(ND*1.));
        Double_t value=0;
        for (Int_t j = 0; j < Nsmooth; j++)
        {
            Double_t rij = dist[j];
            Double_t Wij = Wsm(rij/hi, (int)(rij/hi*0.5*(kernres-1)), kernres, 2.0/(Double_t)(kernres-1), Kernel)*norm;
            value+= Wij* weight[j];
        }
        return value;
    }
}
