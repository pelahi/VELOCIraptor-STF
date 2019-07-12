/*! \file bgfield.cxx
 *  \brief this file contains routines to characterize the background velocity field
 */

//--  Background Velocity Routines

#include "stf.h"

///\name Cell construction using binary kd-tree
//@{

///Decompose Halo/system using KDTree
/*!
    Serial Decomposition of halo/system local (to processor if mpi is invoked).
    The binary kd-tree decomposition using anything other than a tree in configuration space constructed
    using the shannon entropy criteron is allowed but is not \e \b advised.
    Building a physical tree with shannon entropy ensures that regions have uniform M(R)/r or more specifically uniform inter particle spacing
    \todo need to alter pglist array. Much smoother if use particles themselves to find nearest cells
    Instead of storing pglist which is memory intensive (nbodies*ncells) just calculate it when needed.
*/
KDTree* InitializeTreeGrid(Options &opt, const Int_t nbodies, Particle *Part){
    //First rotate into eigenvector frame.
#ifdef SCALING
    Double_t q=1,s=1;
    Coordinate vel;
    Matrix eigvec(0.);
    Int_t i;
    GetGlobalMorphologyWithMass(nbodies,Part, q, s, 1e-3, eigvec);
    //GetGlobalMorphologyWithMass(nbodies, Part, q, s, 1e-3);
    //and rotate velocities to new frame
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,vel)
{
#pragma omp for
#endif
    for (i=0;i<nbodies;i++) {
        for (int j=0;j<3;j++) {
            vel[j]=0;
            for (int k=0;k<3;k++) vel[j]+=eigvec(j,k)*Part[i].GetVelocity(k);
        }
        Part[i].SetVelocity(vel[0],vel[1],vel[2]);
    }
#ifdef USEOPENMP
}
#endif
#endif

    //then build tree
    KDTree *tree;
    if (opt.iverbose>=2) cout<<"Grid system using leaf nodes with maximum size of "<<opt.Ncell<<endl;
    if (opt.gridtype==PHYSGRID) {
        if (opt.iverbose>=2) cout<<"Building Physical Tree using simple spatial extend as splitting criterion"<<endl;
        tree=new KDTree(Part,nbodies,opt.Ncell,tree->TPHYS);
    }
    else if (opt.gridtype==PHYSENGRID) {
        if (opt.iverbose>=2) cout<<"Building physical Tree using minimum shannon entropy as splitting criterion"<<endl;
        //tree=new KDTree(*S,opt.Ncell,tree->TPHYS,tree->KEPAN,100,1);
        tree=new KDTree(Part,nbodies,opt.Ncell,tree->TPHYS,tree->KEPAN,100,1);
    }
    else if (opt.gridtype==PHASEENGRID) {
        if (opt.iverbose>=2) cout<<"Building Phase-space Tree using minimum shannon entropy as splitting criterion"<<endl;
        //tree=new KDTree(*S,opt.Ncell,tree->TPHS,tree->KEPAN,100,1,1);//if phase tree, use entropy criterion with anisotropic kernel
        //if phase tree, use entropy criterion with anisotropic kernel
        tree=new KDTree(Part,nbodies,opt.Ncell,tree->TPHS,tree->KEPAN,100,1,1);
    }
    return tree;
}

///Fills the GridCell struct using KD-Tree initialized by \ref InitializeTreeGrid
void FillTreeGrid(Options &opt, const Int_t nbodies, const Int_t ngrid, KDTree *&tree, Particle *Part, GridCell* &grid)
//void FillTreeGrid(Options &opt, const Int_t nbodies, const Int_t ngrid, KDTree *tree, Particle *Part, GridCell* grid, PartCellNum *pglist)
{
    Int_t gridcount=0,ncount=0;
    //this is used to create tree and search for near neighbours
    Particle *ptemp=new Particle[ngrid];
    int treetype=tree->GetTreeType();
    int ND;
    if (treetype==tree->TPHYS) ND=3;
    else if (treetype==tree->TPHS) ND=6;

    if (opt.iverbose>=2) cout<<"Filling KD-Tree Grid"<<endl;

    //this loop works well for large cells but for small cells may have to replace this with a for loop which goes through every particle
    //and stores nid and size of node, then builds grid based on these values and then places the particles in the cells.
    //for leaf nodes, use starts and ends indices (ie what particles in the system are in the node)
    while (ncount<nbodies) {
        //check if particle is within start and end, if not just means that splitting is not perfect
        //and this particle is left alone. Go to the next particle
        Node *np=(tree->FindLeafNode((Int_t)ncount));
        Int_t start=((LeafNode*)np)->GetStart();
        Int_t end=((LeafNode*)np)->GetEnd();
        //this while loop ensures that if there is a particle out of place, one still proceeds to move through
        //the tree appropriately.
        while (!(ncount>=start&&ncount<end))
        {
            ncount++;
            np=(tree->FindLeafNode((Int_t)ncount));
            start=((LeafNode*)np)->GetStart();
            end=((LeafNode*)np)->GetEnd();
        }

        grid[gridcount].ndim=ND;
        //get center of mass and boundaries of grid cell
        for (int j=0;j<ND;j++) {
            grid[gridcount].xm[j]=0.;
            grid[gridcount].xbl[j]=np->GetBoundary(j,0);
            grid[gridcount].xbu[j]=np->GetBoundary(j,1);
        }
        grid[gridcount].nparts=np->GetCount();
        grid[gridcount].gid=np->GetID();
        grid[gridcount].nindex=new Int_t[grid[gridcount].nparts];

        Double_t mtot=0.;
        for (Int_t k=start,l=0;k<end;k++,l++){
            Int_t id=Part[k].GetID();
            grid[gridcount].nindex[l]=id;
            for (int j=0;j<ND;j++)
                grid[gridcount].xm[j]+=Part[k].GetPosition(j)*Part[k].GetMass();
            mtot+=Part[k].GetMass();
        }
        grid[gridcount].mass=mtot;
        ptemp[gridcount].SetMass(mtot);
        mtot=1.0/mtot;
        for (int j=0;j<ND;j++) grid[gridcount].xm[j]=grid[gridcount].xm[j]*mtot;
        for (int j=0;j<ND;j++) ptemp[gridcount].SetPhase(j,grid[gridcount].xm[j]);
        gridcount++; ncount=end;
    }
    //resets particle order
    delete tree;
    if (opt.iverbose>=2) cout<<"Done."<<endl;
}

//@}

///\name Calculate mean velocity distribution quantities
//@{

///Calculate CM vel of cell
Coordinate* GetCellVel(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngrid, GridCell *grid)
{
    Int_t i;
    Double_t mtot;
    Coordinate *gvel;
    gvel=new Coordinate[ngrid];
    if (opt.iverbose>=2) cout<<"Calculating Grid Mean Velocity"<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,mtot)
{
#pragma omp for
#endif
    for (i=0;i<ngrid;i++) {
        for (int k=0;k<3;k++) gvel[i][k]=0.;
        for (int j=0;j<grid[i].nparts;j++) {
            for (int k=0;k<3;k++) gvel[i][k]+=Part[grid[i].nindex[j]].GetVelocity(k)*Part[grid[i].nindex[j]].GetMass();
        }
        mtot=1.0/grid[i].mass;
        for (int k=0;k<3;k++)gvel[i][k]*=mtot;
    }
#ifdef USEOPENMP
}
#endif
    if (opt.iverbose>=2) cout<<"Done"<<endl;
    return gvel;
}

///Calculate velocity dispersion tensor of cell
Matrix* GetCellVelDisp(Options &opt, const Int_t nbodies, Particle *Part, Int_t ngrid, GridCell *grid, Coordinate *gvel)
{
    Int_t i;
    Double_t mtot;
    Matrix *gveldisp;
    gveldisp=new Matrix[ngrid];
    if (opt.iverbose>=2) cout<<"Calculating Grid Velocity Dispersion"<<endl;
#ifdef USEOPENMP
#pragma omp parallel default(shared) \
private(i,mtot)
{
#pragma omp for
#endif
    for (i=0;i<ngrid;i++) {
        for (int k=0;k<3;k++) for (int l=0;l<3;l++) gveldisp[i](k,l)=0.;
        for (int j=0;j<grid[i].nparts;j++) {
            for (int k=0;k<3;k++) for (int l=0;l<3;l++) gveldisp[i](k,l)+=(Part[grid[i].nindex[j]].GetVelocity(k)-gvel[i][k])*(Part[grid[i].nindex[j]].GetVelocity(l)-gvel[i][l])*Part[grid[i].nindex[j]].GetMass();
        }
        mtot=1.0/grid[i].mass;
        for (int k=0;k<3;k++)for (int l=0;l<3;l++)gveldisp[i](k,l)*=mtot;
    }
#ifdef USEOPENMP
}
#endif
    if (opt.iverbose>=2) cout<<"Done"<<endl;
    return gveldisp;
}

//@}
