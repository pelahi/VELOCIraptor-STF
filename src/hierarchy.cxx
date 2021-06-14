/*! \file hierarchy.cxx
 *  \brief this file contains routines determine the hierarchy of objects 
 */

#include "stf.h"

/// \name Routines used to determine substructure hierarchy
//@{
Int_t GetHierarchy(Options &opt, Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t* stype)
{
    if (opt.iverbose) cout<<"Getting Hierarchy "<<ngroups<<endl;
    Int_t nhierarchy=1;
    StrucLevelData *ppsldata,**papsldata;
    ppsldata=psldata;
    while (ppsldata->nextlevel!=NULL){nhierarchy++;ppsldata=ppsldata->nextlevel;}
    for (Int_t i=1;i<=ngroups;i++) nsub[i]=0;
    for (Int_t i=1;i<=ngroups;i++) parentgid[i]=GROUPNOPARENT;
    for (Int_t i=1;i<=ngroups;i++) uparentgid[i]=GROUPNOPARENT;
    ppsldata=psldata;
    papsldata=new StrucLevelData*[nhierarchy];
    nhierarchy=0;
    while (ppsldata!=NULL) {papsldata[nhierarchy++]=ppsldata;ppsldata=ppsldata->nextlevel;}
    for (int i=nhierarchy-1;i>=1;i--){
        //store number of substructures
        for (int j=1;j<=papsldata[i]->nsinlevel;j++) {
               if (papsldata[i]->gidparenthead[j]!=NULL&&papsldata[i]->gidparenthead[j]!=papsldata[i]->gidhead[j]) nsub[*(papsldata[i]->gidparenthead[j])]++;
        }
        //then add these to parent substructure
        for (int j=1;j<=papsldata[i]->nsinlevel;j++) {
            if (papsldata[i]->gidparenthead[j]!=NULL&&papsldata[i]->gidparenthead[j]!=papsldata[i]->gidhead[j]){
                nsub[*(papsldata[i]->gidparenthead[j])]+=nsub[*(papsldata[i]->gidhead[j])];
                parentgid[*(papsldata[i]->gidhead[j])]=*(papsldata[i]->gidparenthead[j]);
            }
            if (papsldata[i]->giduberparenthead[j]!=NULL &&papsldata[i]->gidparenthead[j]!=papsldata[i]->gidhead[j]) uparentgid[*(papsldata[i]->gidhead[j])]=*(papsldata[i]->giduberparenthead[j]);
            stype[*(papsldata[i]->gidhead[j])]=(papsldata[i]->stypeinlevel[j]);
        }
    }
    //store field structures (top treel level) types
    for (int j=1;j<=papsldata[0]->nsinlevel;j++) {
        stype[*(papsldata[0]->gidhead[j])]=(papsldata[0]->stypeinlevel[j]);
    }
    for (int i=0;i<nhierarchy;i++)papsldata[i]=NULL;
    delete[] papsldata;
    if(opt.iverbose) cout<<"Done"<<endl;
    return nhierarchy;
}

void CopyHierarchyToPropData(Options &opt, PropData *pdata, Int_t ngroups, Int_t *nsub, Int_t *parentgid, Int_t *uparentgid, Int_t* stype)
{
    Int_t i,haloidoffset=0;
#ifdef USEMPI
    for (int j=0;j<ThisTask;j++)haloidoffset+=mpi_ngroups[j];
#endif
#ifdef USEOPENMP
#pragma omp parallel default(shared)  \
private(i)
{
    #pragma omp for nowait
#endif
    for (i=1;i<=ngroups;i++) {
        pdata[i].haloid=opt.snapshotvalue+i;
        //if using mpi than ids must be offset
#ifdef USEMPI
        pdata[i].haloid+=haloidoffset;
#endif
        if (parentgid[i]!=GROUPNOPARENT) {
            pdata[i].hostid=opt.snapshotvalue+uparentgid[i];
            pdata[i].directhostid=opt.snapshotvalue+parentgid[i];
            if(opt.iKeepFOF && uparentgid[i]<=opt.num3dfof) pdata[i].hostfofid=opt.snapshotvalue+uparentgid[i];
            if(opt.iKeepFOF && uparentgid[i]>opt.num3dfof) pdata[i].hostfofid=GROUPNOPARENT;
#ifdef USEMPI
            pdata[i].hostid+=haloidoffset;
            pdata[i].directhostid+=haloidoffset;
            if(opt.iKeepFOF && uparentgid[i]<=opt.num3dfof) pdata[i].hostfofid+=haloidoffset;
#endif
        }
        else {
            pdata[i].hostid=GROUPNOPARENT;
            pdata[i].directhostid=GROUPNOPARENT;
            pdata[i].hostfofid=GROUPNOPARENT;
        }
        pdata[i].stype=stype[i];
        pdata[i].numsubs=nsub[i];
    }
#ifdef USEOPENMP
}
#endif
}
//@}
