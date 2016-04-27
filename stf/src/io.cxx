/*! \file io.cxx
 *  \brief this file contains routines for io
 */

//-- IO

#include "stf.h"

#include "gadgetitems.h"
#include "tipsy_structs.h"
#include "endianutils.h"
#ifdef USEHDF
#include "hdfitems.h"
#endif
#include "ramsesitems.h"

///Checks if file exits by attempting to get the file attributes
///If success file obviously exists.
///If failure may mean that we don't have permission to access the folder which contains this file or doesn't exist. 
///If we need to do that level of checking, lookup return values of stat which will give you more details on why stat failed.
bool FileExists(const char *fname) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;
  intStat = stat(fname,&stFileInfo);
  if(intStat == 0) return  true;
  else return false;
}

///\name Read particle data files
//@{

///Reads the header structure if its a tipsy file or get_nbodies if gadget
Int_t ReadHeader(Options &opt){
    InitEndian();
    if (opt.inputtype==IOTIPSY) {
        struct dump tipsyheader;
        fstream Ftip(opt.fname, ios::in | ios::binary);
        if (!Ftip){cerr<<"ERROR: Unable to open " <<opt.fname<<endl;exit(8);}
        else cout<<"Reading tipsy format from "<<opt.fname<<endl;
        Ftip.read((char*)&tipsyheader,sizeof(dump));
        if (opt.partsearchtype==PSTALL) return tipsyheader.nbodies;
        else if (opt.partsearchtype==PSTDARK) return tipsyheader.ndark;
        else if (opt.partsearchtype==PSTGAS) return tipsyheader.nsph;
        else if (opt.partsearchtype==PSTSTAR) return tipsyheader.nstar;
    }
    else if (opt.inputtype==IOGADGET) {
        if (opt.partsearchtype==PSTALL) return get_nbodies(opt.fname);
        else if (opt.partsearchtype==PSTDARK) return get_nbodies(opt.fname,-2);
        else if (opt.partsearchtype==PSTGAS) return get_nbodies(opt.fname,GGASTYPE);
        else if (opt.partsearchtype==PSTSTAR) return get_nbodies(opt.fname,GSTARTYPE);
        else if (opt.partsearchtype==PSTBH) return get_nbodies(opt.fname,GBHTYPE);
    }
    else if (opt.inputtype==IORAMSES) return RAMSES_get_nbodies(opt.fname,opt.partsearchtype,opt);
#ifdef USEHDF
    else if (opt.inputtype==IOHDF) return HDF_get_nbodies(opt.fname,opt.partsearchtype);
#endif
    return 0;
}

///Reads particle data
///To add a new interface simply alter this to include the appropriate user written call
void ReadData(Options &opt, Particle *&Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons)
{
    InitEndian();
    if(opt.inputtype==IOTIPSY) ReadTipsy(opt,Part,nbodies, Pbaryons, nbaryons);
    else if (opt.inputtype==IOGADGET) ReadGadget(opt,Part,nbodies, Pbaryons, nbaryons);
    else if (opt.inputtype==IORAMSES) ReadRamses(opt,Part,nbodies, Pbaryons, nbaryons);
#ifdef USEHDF
    else if (opt.inputtype==IOHDF) ReadHDF(opt,Part,nbodies, Pbaryons, nbaryons);
#endif

#ifdef USEMPI
    MPIAdjustDomain(opt);
#endif
}

//@}

///\name Read STF data files
//@{

///Read local velocity density
void ReadLocalVelocityDensity(Options &opt, const Int_t nbodies, Particle * Part){
    Int_t tempi;
    Double_t tempd;
    fstream Fin;
    char fname[1000];
    //set filename appropriate to mpi thread if necessary
#ifdef USEMPI
    sprintf(fname,"%s.%d",opt.smname,ThisTask);
#else
    sprintf(fname,"%s",opt.smname);
#endif

    cout<<"Reading smooth density data from "<<fname<<endl;
    if (opt.ibinaryout) {
        Fin.open(fname,ios::in| ios::binary);
        Fin.read((char*)&tempi,sizeof(Int_t));
        if (tempi!=nbodies) {
            cerr<<"File "<<fname<<" contains incorrect number of particles. Exiting\n";
            exit(9);
        }
        for(Int_t i=0;i<nbodies;i++) {Fin.read((char*)&tempd,sizeof(Double_t));Part[i].SetDensity(tempd);}
    }
    else {
        Fin.open(fname,ios::in);
        Fin>>tempi;
        if (tempi!=nbodies) {
            cerr<<"File "<<fname<<" contains incorrect number of particles. Exiting\n";
            exit(9);
        }
        for(Int_t i=0;i<nbodies;i++) {Fin>>tempd;Part[i].SetDensity(tempd);}
    }
    cout<<"Done"<<endl;
    Fin.close();
}

//@}

/// \name Write STF data files for intermediate steps 
//@{

///Writes local velocity density of each particle to a file
void WriteLocalVelocityDensity(Options &opt, const Int_t nbodies, Particle * Part){
    fstream Fout;
    char fname[1000];
#ifdef USEMPI
    if(opt.smname==NULL) sprintf(fname,"%s.smdata.%d",opt.outname,ThisTask);
    else sprintf(fname,"%s.%d",opt.smname,ThisTask);
#else
    if(opt.smname==NULL) sprintf(fname,"%s.smdata",opt.outname);
    else sprintf(fname,"%s",opt.smname);
#endif
    if (opt.ibinaryout) {
        Fout.open(fname,ios::out|ios::binary);
        Fout.write((char*)&nbodies,sizeof(Int_t));
        Double_t tempd;
        for(Int_t i=0;i<nbodies;i++) {tempd=Part[i].GetDensity();Fout.write((char*)&tempd,sizeof(Double_t));}
    }
    if (opt.ibinaryout) {
        Fout.open(fname,ios::out);
        Fout<<nbodies<<endl;
        Fout<<scientific<<setprecision(10);
        for(Int_t i=0;i<nbodies;i++)Fout<<Part[i].GetDensity()<<endl;
    }
    Fout.close();
}

//@}

///\name FOF outputs
//@{

/*! Writes a tipsy formatted fof.grp array file that contains the number of particles first then for each particle the group id of that particle
    group zero is untagged particles. \n
*/
void WriteFOF(Options &opt, const Int_t nbodies, Int_t *pfof){
    fstream Fout;
    char fname[1000];
    sprintf(fname,"%s.fof.grp",opt.outname);
    cout<<"saving fof data to "<<fname<<endl;
    Fout.open(fname,ios::out);
    if (opt.partsearchtype==PSTALL) {
        Fout<<nbodies<<endl;
        for (Int_t i=0;i<nbodies;i++) Fout<<pfof[i]<<endl;
    }
    else if (opt.partsearchtype==PSTDARK) {
        Int_t nt=0;
        for (int i=0;i<NPARTTYPES;i++) nt+=opt.numpart[i];
        Fout<<nt<<endl;
        for (Int_t i=0;i<opt.numpart[GASTYPE];i++) Fout<<0<<endl;
        for (Int_t i=0;i<nbodies;i++) Fout<<pfof[i]<<endl;
        for (Int_t i=0;i<opt.numpart[STARTYPE];i++) Fout<<0<<endl;
    }
    else if (opt.partsearchtype==PSTSTAR) {
        Int_t nt=0;
        for (int i=0;i<NPARTTYPES;i++) nt+=opt.numpart[i];
        Fout<<nt<<endl;
        for (Int_t i=0;i<opt.numpart[GASTYPE];i++) Fout<<0<<endl;
        for (Int_t i=0;i<opt.numpart[DARKTYPE];i++) Fout<<0<<endl;
        for (Int_t i=0;i<nbodies;i++) Fout<<pfof[i]<<endl;
    }
    else if (opt.partsearchtype==PSTGAS) {
        Int_t nt=0;
        for (int i=0;i<NPARTTYPES;i++) nt+=opt.numpart[i];
        Fout<<nt<<endl;
        for (Int_t i=0;i<nbodies;i++) Fout<<pfof[i]<<endl;
        for (Int_t i=0;i<opt.numpart[DARKTYPE];i++) Fout<<0<<endl;
        for (Int_t i=0;i<opt.numpart[STARTYPE];i++) Fout<<0<<endl;
    }
    Fout.close();
    cout<<"Done"<<endl;
}

/*! Writes a particle group list array file that contains the total number of groups,
    local number of groups (if using MPI) and group id followed by number of particles
    in that group and particle ids in the group
    For MPI, each thread writes it its own file (ie: parallel write)
    \todo Must check if parallel output is not an issue since on shared memory machine, write would probably write
    to the same hard drive, whereas on a cluster, each system could in principle write to different drive.
*/
void WritePGListIndex(Options &opt, const Int_t ngroups, const Int_t ng, Int_t *numingroup, Int_t **pglist){
    fstream Fout;
    char fname[1000];
    Int_t noffset=0,ngtot=0;
#ifdef USEMPI
    sprintf(fname,"%s.pglist.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.pglist",opt.outname);
#endif
    cout<<"saving fof data to "<<fname<<endl;
    Fout.open(fname,ios::out);
#ifdef USEMPI
    for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
    Fout<<ngtot<<" "<<ngroups<<endl;
    for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
#else
    Fout<<ngroups<<" "<<ngroups<<endl;
#endif
    for (Int_t i=1;i<=ngroups;i++) {
#ifdef USEMPI
        Fout<<i+noffset<<" "<<numingroup[i]<<" ";
#else
        Fout<<i<<" "<<numingroup[i]<<" ";
#endif
        for (Int_t j=0;j<numingroup[i];j++) 
#ifdef USEMPI
            //Fout<<mpi_indexlist[pglist[i][j]]<<" ";
            Fout<<pglist[i][j]<<" ";
#else
            Fout<<pglist[i][j]<<" ";
#endif
        Fout<<endl;
    }
    for (Int_t i=1;i<ng;i++) delete[] pglist;
    delete[] pglist;
    delete[] numingroup;
    cout<<"Done"<<endl;
    Fout.close();
}
void WritePGList(Options &opt, const Int_t ngroups, const Int_t ng, Int_t *numingroup, Int_t **pglist, Int_t *ids){
    fstream Fout;
    char fname[1000];
    Int_t noffset=0,ngtot=0;
#ifdef USEMPI
    sprintf(fname,"%s.pglist.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.pglist",opt.outname);
#endif
    cout<<"saving fof data to "<<fname<<endl;
    Fout.open(fname,ios::out);

#ifdef USEMPI
    for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
    Fout<<ngtot<<" "<<ngroups<<endl;
    for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
#else
    Fout<<ngroups<<" "<<ngroups<<endl;
#endif
    for (Int_t i=1;i<=ngroups;i++) {
#ifdef USEMPI
        Fout<<i+noffset<<" "<<numingroup[i]<<" ";
#else
        Fout<<i<<" "<<numingroup[i]<<" ";
#endif
        for (Int_t j=0;j<numingroup[i];j++) 
#ifdef USEMPI
            //Fout<<mpi_idlist[pglist[i][j]]<<" ";
            Fout<<ids[pglist[i][j]]<<" ";
#else
            Fout<<ids[pglist[i][j]]<<" ";
#endif
        Fout<<endl;
    }
    for (Int_t i=1;i<ng;i++) delete[] pglist;
    delete[] pglist;
    delete[] numingroup;
    cout<<"Done"<<endl;
    Fout.close();
}

void WriteGroupCatalog(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, Particle *Part, Int_t nadditional){
    fstream Fout,Fout2,Fout3;
    char fname[500];
    char fname2[500];
    char fname3[500];
    Int_t noffset=0,ngtot=0,nids=0,nidstot,nuids=0,nuidstot;
    Int_t *offset;

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

#ifdef USEMPI
    sprintf(fname,"%s.catalog_groups.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.catalog_groups",opt.outname);
#endif

    cout<<"saving group catalog to "<<fname<<endl;
    if (opt.ibinaryout) Fout.open(fname,ios::out|ios::binary);
    else Fout.open(fname,ios::out);


#ifdef USEMPI
    for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
#else
    ngtot=ngroups+nadditional;//useful if outputing field halos
#endif
    //write header
    if (opt.ibinaryout) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&ngroups,sizeof(Int_t));
        Fout.write((char*)&ngtot,sizeof(Int_t));
        //Fout.write((char*)&mpi_ngroups[ThisTask],sizeof(Int_t));
    }
    else{
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ngroups<<" "<<ngtot<<endl;
    }

    //write group size
    if (opt.ibinaryout) Fout.write((char*)&numingroup[1],sizeof(Int_t)*ngroups);
    else for (Int_t i=1;i<=ngroups;i++) Fout<<numingroup[i]<<endl;


    //Write offsets for bound and unbound particles
    offset=new Int_t[ngroups+1];
    offset[1]=0;
    //note before had offsets at numingroup but to account for unbound particles use value of pglist at numingroup 
    for (Int_t i=2;i<=ngroups;i++) offset[i]=offset[i-1]+pglist[i-1][numingroup[i-1]];

    if (opt.ibinaryout) Fout.write((char*)&offset[1],sizeof(Int_t)*ngroups);
    else for (Int_t i=1;i<=ngroups;i++) Fout<<offset[i]<<endl;

    //position of unbound particle
    for (Int_t i=2;i<=ngroups;i++) offset[i]=offset[i-1]+numingroup[i-1]-pglist[i-1][numingroup[i-1]];
    if (opt.ibinaryout) Fout.write((char*)&offset[1],sizeof(Int_t)*ngroups);
    else for (Int_t i=1;i<=ngroups;i++) Fout<<offset[i]<<endl;

    delete[] offset;
    Fout.close();

    //now write pid files
#ifdef USEMPI
    sprintf(fname,"%s.catalog_particles.%d",opt.outname,ThisTask);
    sprintf(fname3,"%s.catalog_particles.unbound.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.catalog_particles",opt.outname);
    sprintf(fname3,"%s.catalog_particles.unbound",opt.outname);
#endif
    cout<<"saving particle catalog to "<<fname<<endl;

    if (opt.ibinaryout) {
        Fout.open(fname,ios::out|ios::binary);
        Fout3.open(fname3,ios::out|ios::binary);
    }
    else {
        Fout.open(fname,ios::out);
        Fout3.open(fname3,ios::out);
    }

    //see above regarding unbound particle
    //for (Int_t i=1;i<=ngroups;i++) nids+=numingroup[i];
    for (Int_t i=1;i<=ngroups;i++) {nids+=pglist[i][numingroup[i]];nuids+=numingroup[i]-pglist[i][numingroup[i]];}
#ifdef USEMPI
#ifdef LONGINT
    MPI_Allreduce(&nids, &nidstot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nuids, &nuidstot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Allreduce(&nids, &nidstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nuids, &nuidstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
#else
    nidstot=nids;
    nuidstot=nuids;
#endif

    //write header
    if (opt.ibinaryout) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&nids,sizeof(Int_t));
        Fout.write((char*)&nidstot,sizeof(Int_t));

        Fout3.write((char*)&ThisTask,sizeof(int));
        Fout3.write((char*)&NProcs,sizeof(int));
        Fout3.write((char*)&nuids,sizeof(Int_t));
        Fout3.write((char*)&nuidstot,sizeof(Int_t));
    }
    else {
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<nids<<" "<<nidstot<<endl;

        Fout3<<ThisTask<<" "<<NProcs<<endl;
        Fout3<<nuids<<" "<<nuidstot<<endl;
    }

    Int_t *idval;
    if (nids>0) {
    idval=new Int_t[nids];
    nids=0;
    for (Int_t i=1;i<=ngroups;i++) 
        for (Int_t j=0;j<pglist[i][numingroup[i]];j++) 
            idval[nids++]=Part[pglist[i][j]].GetPID();
    if (opt.ibinaryout) Fout.write((char*)idval,sizeof(Int_t)*nids);
    else for (Int_t i=0;i<nids;i++) Fout<<idval[i]<<endl;
    delete[] idval;
    }
    Fout.close();

    if (nuids>0) {
    idval=new Int_t[nuids];
    nuids=0;
    for (Int_t i=1;i<=ngroups;i++) 
        for (Int_t j=pglist[i][numingroup[i]];j<numingroup[i];j++) 
            idval[nuids++]=Part[pglist[i][j]].GetPID();
    if (opt.ibinaryout) Fout3.write((char*)idval,sizeof(Int_t)*nuids);
    else for (Int_t i=0;i<nuids;i++) Fout3<<idval[i]<<endl;
    delete[] idval;
    }
    Fout3.close();
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

///if particles are separately searched (i.e. \ref Options.iBaryonSearch is set) then produce list of particle types
void WriteGroupPartType(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, Particle *Part){
    fstream Fout,Fout2;
    char fname[2000];
    char fname2[2000];
    Int_t noffset=0,ngtot=0,nids=0,nidstot,nuids=0,nuidstot=0;
    Int_t *offset;
    int *typeval;

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

#ifdef USEMPI
    sprintf(fname,"%s.catalog_parttypes.%d",opt.outname,ThisTask);
    sprintf(fname2,"%s.catalog_parttypes.unbound.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.catalog_parttypes",opt.outname);
    sprintf(fname2,"%s.catalog_parttypes.unbound",opt.outname);
#endif
    cout<<"saving particle type info to "<<fname<<endl;


    if (opt.ibinaryout) {
        Fout.open(fname,ios::out|ios::binary);
        Fout2.open(fname2,ios::out|ios::binary);
    }
    else {
        Fout.open(fname,ios::out);
        Fout2.open(fname2,ios::out);
    }

    for (Int_t i=1;i<=ngroups;i++) {nids+=pglist[i][numingroup[i]];nuids+=numingroup[i]-pglist[i][numingroup[i]];}
#ifdef USEMPI
#ifdef LONGINT
    MPI_Allreduce(&nids, &nidstot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nuids, &nuidstot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Allreduce(&nids, &nidstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&nuids, &nuidstot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
#else
    nidstot=nids;
    nuidstot=nuids;
#endif

    //write header
    if (opt.ibinaryout) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&nids,sizeof(Int_t));
        Fout.write((char*)&nidstot,sizeof(Int_t));

        Fout2.write((char*)&ThisTask,sizeof(int));
        Fout2.write((char*)&NProcs,sizeof(int));
        Fout2.write((char*)&nuids,sizeof(Int_t));
        Fout2.write((char*)&nuidstot,sizeof(Int_t));
    }
    else {
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<nids<<" "<<nidstot<<endl;

        Fout2<<ThisTask<<" "<<NProcs<<endl;
        Fout2<<nuids<<" "<<nuidstot<<endl;
    }

    if (nids>0) {
    typeval=new int[nids];
    nids=0;
    for (Int_t i=1;i<=ngroups;i++) 
        for (Int_t j=0;j<pglist[i][numingroup[i]];j++) 
            typeval[nids++]=Part[pglist[i][j]].GetType();
    if (opt.ibinaryout) Fout.write((char*)typeval,sizeof(int)*nids);
    else for (Int_t i=0;i<nids;i++) Fout<<typeval[i]<<endl;
    delete[] typeval;
    }
    Fout.close();

    if (nuids>0) {
    typeval=new int[nuids];
    nuids=0;
    for (Int_t i=1;i<=ngroups;i++) 
        for (Int_t j=pglist[i][numingroup[i]];j<numingroup[i];j++) 
            typeval[nuids++]=Part[pglist[i][j]].GetType();
    if (opt.ibinaryout) Fout2.write((char*)typeval,sizeof(int)*nuids);
    else for (Int_t i=0;i<nuids;i++) Fout2<<typeval[i]<<endl;
    delete[] typeval;
    }
    Fout2.close();
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//@}


///\name Final outputs such as properties and output that can be used to construct merger trees and substructure hierarchy
//@{
///Writes the bulk properties of the substructures
void WriteProperties(Options &opt, const Int_t ngroups, PropData *pdata){
    fstream Fout;
    char fname[1000];
    Int_t ngtot=0, noffset=0;
    //list the header info
    vector<string> headerdatainfo;
    headerdatainfo.push_back("ID");
    headerdatainfo.push_back("ID_mbp");
    headerdatainfo.push_back("hostHaloID");
    headerdatainfo.push_back("numSubStruct");
    headerdatainfo.push_back("npart");
    headerdatainfo.push_back("Mvir");
    headerdatainfo.push_back("Xc");
    headerdatainfo.push_back("Yc");
    headerdatainfo.push_back("Zc");
    headerdatainfo.push_back("Xcmbp");
    headerdatainfo.push_back("Ycmbp");
    headerdatainfo.push_back("Zcmbp");
    headerdatainfo.push_back("VXc");
    headerdatainfo.push_back("VYc");
    headerdatainfo.push_back("VZc");
    headerdatainfo.push_back("VXcmbp");
    headerdatainfo.push_back("VYcmbp");
    headerdatainfo.push_back("VZcmbp");
    headerdatainfo.push_back("Mass_tot");
    headerdatainfo.push_back("Mass_FOF");
    headerdatainfo.push_back("Mass_200mean");
    headerdatainfo.push_back("Mass_200crit");
    headerdatainfo.push_back("Mass_BN97");
    headerdatainfo.push_back("Efrac");
    headerdatainfo.push_back("Rvir");
    headerdatainfo.push_back("R_size");
    headerdatainfo.push_back("R_200mean");
    headerdatainfo.push_back("R_200crit");
    headerdatainfo.push_back("R_BN97");
    headerdatainfo.push_back("R_HalfMass");
    headerdatainfo.push_back("Rmax");
    headerdatainfo.push_back("Vmax");
    headerdatainfo.push_back("sigV");
    headerdatainfo.push_back("veldisp_xx");
    headerdatainfo.push_back("veldisp_xy");
    headerdatainfo.push_back("veldisp_xz");
    headerdatainfo.push_back("veldisp_yx");
    headerdatainfo.push_back("veldisp_yy");
    headerdatainfo.push_back("veldisp_yz");
    headerdatainfo.push_back("veldisp_zx");
    headerdatainfo.push_back("veldisp_zy");
    headerdatainfo.push_back("veldisp_zz");
    headerdatainfo.push_back("lambda_B");
    headerdatainfo.push_back("Lx");
    headerdatainfo.push_back("Ly");
    headerdatainfo.push_back("Lz");
    headerdatainfo.push_back("q");
    headerdatainfo.push_back("s");
    headerdatainfo.push_back("eig_xx");
    headerdatainfo.push_back("eig_xy");
    headerdatainfo.push_back("eig_xz");
    headerdatainfo.push_back("eig_yx");
    headerdatainfo.push_back("eig_yy");
    headerdatainfo.push_back("eig_yz");
    headerdatainfo.push_back("eig_zx");
    headerdatainfo.push_back("eig_zy");
    headerdatainfo.push_back("eig_zz");
    headerdatainfo.push_back("cNFW");
    headerdatainfo.push_back("Krot");
    headerdatainfo.push_back("Ekin");
    headerdatainfo.push_back("Epot");
#ifdef GASON
    headerdatainfo.push_back("n_gas");
    headerdatainfo.push_back("M_gas");
    headerdatainfo.push_back("Xc_gas");
    headerdatainfo.push_back("Yc_gas");
    headerdatainfo.push_back("Zc_gas");
    headerdatainfo.push_back("VXc_gas");
    headerdatainfo.push_back("VYc_gas");
    headerdatainfo.push_back("VZc_gas");
    headerdatainfo.push_back("Efrac_gas");
    headerdatainfo.push_back("R_HalfMass_gas");
    headerdatainfo.push_back("veldisp_xx_gas");
    headerdatainfo.push_back("veldisp_xy_gas");
    headerdatainfo.push_back("veldisp_xz_gas");
    headerdatainfo.push_back("veldisp_yx_gas");
    headerdatainfo.push_back("veldisp_yy_gas");
    headerdatainfo.push_back("veldisp_yz_gas");
    headerdatainfo.push_back("veldisp_zx_gas");
    headerdatainfo.push_back("veldisp_zy_gas");
    headerdatainfo.push_back("veldisp_zz_gas");
    headerdatainfo.push_back("Lx_gas");
    headerdatainfo.push_back("Ly_gas");
    headerdatainfo.push_back("Lz_gas");
    headerdatainfo.push_back("q_gas");
    headerdatainfo.push_back("s_gas");
    headerdatainfo.push_back("eig_xx_gas");
    headerdatainfo.push_back("eig_xy_gas");
    headerdatainfo.push_back("eig_xz_gas");
    headerdatainfo.push_back("eig_yx_gas");
    headerdatainfo.push_back("eig_yy_gas");
    headerdatainfo.push_back("eig_yz_gas");
    headerdatainfo.push_back("eig_zx_gas");
    headerdatainfo.push_back("eig_zy_gas");
    headerdatainfo.push_back("eig_zz_gas");
    headerdatainfo.push_back("Krot_gas");
    headerdatainfo.push_back("T_gas");
#ifdef STARON
    headerdatainfo.push_back("Zmet_gas");
    headerdatainfo.push_back("SFR_gas");
#endif
#endif
#ifdef STARON
    headerdatainfo.push_back("n_star");
    headerdatainfo.push_back("M_star");
    headerdatainfo.push_back("Xc_star");
    headerdatainfo.push_back("Yc_star");
    headerdatainfo.push_back("Zc_star");
    headerdatainfo.push_back("VXc_star");
    headerdatainfo.push_back("VYc_star");
    headerdatainfo.push_back("VZc_star");
    headerdatainfo.push_back("Efrac_star");
    headerdatainfo.push_back("R_HalfMass_star");
    headerdatainfo.push_back("veldisp_xx_star");
    headerdatainfo.push_back("veldisp_xy_star");
    headerdatainfo.push_back("veldisp_xz_star");
    headerdatainfo.push_back("veldisp_yx_star");
    headerdatainfo.push_back("veldisp_yy_star");
    headerdatainfo.push_back("veldisp_yz_star");
    headerdatainfo.push_back("veldisp_zx_star");
    headerdatainfo.push_back("veldisp_zy_star");
    headerdatainfo.push_back("veldisp_zz_star");
    headerdatainfo.push_back("Lx_star");
    headerdatainfo.push_back("Ly_star");
    headerdatainfo.push_back("Lz_star");
    headerdatainfo.push_back("q_star");
    headerdatainfo.push_back("s_star");
    headerdatainfo.push_back("eig_xx_star");
    headerdatainfo.push_back("eig_xy_star");
    headerdatainfo.push_back("eig_xz_star");
    headerdatainfo.push_back("eig_yx_star");
    headerdatainfo.push_back("eig_yy_star");
    headerdatainfo.push_back("eig_yz_star");
    headerdatainfo.push_back("eig_zx_star");
    headerdatainfo.push_back("eig_zy_star");
    headerdatainfo.push_back("eig_zz_star");
    headerdatainfo.push_back("Krot_star");
    headerdatainfo.push_back("tage_star");
    headerdatainfo.push_back("Zmet_star");
#endif
    ///\todo CODE ALTERED SO THAT HEADER IS ALWAYS 4 ints, this file, number of files, number of groups, total number of groups
    //Int_t noffset=0,ngtot=0;
    
#ifdef USEMPI
    sprintf(fname,"%s.properties.%d",opt.outname,ThisTask);
    for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
    for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
#else
    sprintf(fname,"%s.properties",opt.outname);
    int ThisTask=0,NProcs=1;
    ngtot=ngroups;
#endif
    cout<<"saving property data to "<<fname<<endl;

    //write header
    if (opt.ibinaryout) {
        Fout.open(fname,ios::out|ios::binary);
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&ngroups,sizeof(Int_t));
        Fout.write((char*)&ngtot,sizeof(Int_t));
        ///\todo ADD string containing information of what is in output since this will possibly change with time
    }
    else {
        Fout.open(fname,ios::out);
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ngroups<<" "<<ngtot<<endl;
        for (Int_t i=0;i<headerdatainfo.size();i++) Fout<<headerdatainfo[i]<<"("<<i+1<<") ";Fout<<endl;
        Fout<<setprecision(10);
    }
    //if need to convert from physical back to comoving
    if (opt.icomoveunit) for (Int_t i=1;i<=ngroups;i++) pdata[i].ConverttoComove(opt);
    
    long long idbound;
    //for ensuring downgrade of precision as subfind uses floats when storing values save for Mvir (??why??)
    float value,ctemp[3],mtemp[9];
    double dvalue;
    int ivalue;
    for (Int_t i=1;i<=ngroups;i++) {
        if (opt.ibinaryout) {
            idbound=pdata[i].ibound;
            Fout.write((char*)&idbound,sizeof(long long));
            dvalue=pdata[i].gMvir;
            Fout.write((char*)&dvalue,sizeof(dvalue));
            ivalue=pdata[i].num;
            Fout.write((char*)&ivalue,sizeof(ivalue));
            for (int k=0;k<3;k++) ctemp[k]=pdata[i].gcm[k];
            Fout.write((char*)ctemp,sizeof(float)*3);
            for (int k=0;k<3;k++) ctemp[k]=pdata[i].gpos[k];
            Fout.write((char*)ctemp,sizeof(float)*3);
            for (int k=0;k<3;k++) ctemp[k]=pdata[i].gcmvel[k];
            Fout.write((char*)ctemp,sizeof(float)*3);
            for (int k=0;k<3;k++) ctemp[k]=pdata[i].gvel[k];
            Fout.write((char*)ctemp,sizeof(float)*3);
            value=pdata[i].gRvir;
            Fout.write((char*)&value,sizeof(value));
            value=pdata[i].gRmbp;
            Fout.write((char*)&value,sizeof(value));
            value=pdata[i].gRmaxvel;
            Fout.write((char*)&value,sizeof(value));
            value=pdata[i].gmaxvel;
            Fout.write((char*)&value,sizeof(value));
            value=pdata[i].gsigma_v;
            Fout.write((char*)&value,sizeof(value));
            for (int k=0;k<3;k++) ctemp[k]=pdata[i].gJ[k];
            Fout.write((char*)ctemp,sizeof(float)*3);
            value=pdata[i].gq;
            Fout.write((char*)&value,sizeof(value));
            value=pdata[i].gs;
            Fout.write((char*)&value,sizeof(value));
            for (int k=0;k<3;k++) for (int n=0;n<3;n++) mtemp[k*3+n]=pdata[i].geigvec(k,n);
            Fout.write((char*)mtemp,sizeof(float)*9);
            //afterwards would be padding for 8 more floats (chars), add extra stuff like total mass, mass enclosed in Rmax, T, Pot, Efrac
            dvalue=pdata[i].gmass;
            Fout.write((char*)&dvalue,sizeof(dvalue));
            dvalue=pdata[i].gMmaxvel;
            Fout.write((char*)&dvalue,sizeof(dvalue));
            value=pdata[i].T;
            Fout.write((char*)&value,sizeof(value));
            value=pdata[i].Pot;
            Fout.write((char*)&value,sizeof(value));
            value=pdata[i].Efrac;
            Fout.write((char*)&value,sizeof(value));
            value=0;
            for (int k=0;k<1;k++) Fout.write((char*)&value,sizeof(value));//for fact that 2*2*4+3*4
        }
        else {
            Fout<<pdata[i].haloid<<" ";
            Fout<<pdata[i].ibound<<" ";
            Fout<<pdata[i].hostid<<" ";
            Fout<<pdata[i].numsubs<<" ";
            Fout<<pdata[i].num<<" ";
            Fout<<pdata[i].gMvir<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].gcm[k]<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].gpos[k]<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].gcmvel[k]<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].gvel[k]<<" ";
            Fout<<pdata[i].gmass<<" ";
            Fout<<pdata[i].gMFOF<<" ";
            Fout<<pdata[i].gM200m<<" ";
            Fout<<pdata[i].gM200c<<" ";
            Fout<<pdata[i].gMvir<<" ";
            Fout<<pdata[i].Efrac<<" ";
            Fout<<pdata[i].gRvir<<" ";
            Fout<<pdata[i].gsize<<" ";
            Fout<<pdata[i].gR200m<<" ";
            Fout<<pdata[i].gR200c<<" ";
            Fout<<pdata[i].gRvir<<" ";
            Fout<<pdata[i].gRhalfmass<<" ";
            Fout<<pdata[i].gRmaxvel<<" ";
            Fout<<pdata[i].gmaxvel<<" ";
            Fout<<pdata[i].gsigma_v<<" ";
            for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<pdata[i].gveldisp(k,n)<<" ";
            Fout<<pdata[i].glambda_B<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].gJ[k]<<" ";
            Fout<<pdata[i].gq<<" ";
            Fout<<pdata[i].gs<<" ";
            for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<pdata[i].geigvec(k,n)<<" ";
            Fout<<pdata[i].cNFW<<" ";
            Fout<<pdata[i].Krot<<" ";
            Fout<<pdata[i].T<<" ";
            Fout<<pdata[i].Pot<<" ";
#ifdef GASON
            Fout<<pdata[i].n_gas<<" ";
            Fout<<pdata[i].M_gas<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].cm_gas[k]<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].cmvel_gas[k]<<" ";
            Fout<<pdata[i].Efrac_gas<<" ";
            Fout<<pdata[i].Rhalfmass_gas<<" ";
            for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<pdata[i].veldisp_gas(k,n)<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].L_gas[k]<<" ";
            Fout<<pdata[i].q_gas<<" ";
            Fout<<pdata[i].s_gas<<" ";
            for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<pdata[i].eigvec_gas(k,n)<<" ";
            Fout<<pdata[i].Krot_gas<<" ";
            Fout<<pdata[i].Temp_gas<<" ";
#ifdef STARON
            Fout<<pdata[i].Z_gas<<" ";
            Fout<<pdata[i].SFR_gas<<" ";
#endif
#endif

#ifdef STARON
            Fout<<pdata[i].n_star<<" ";
            Fout<<pdata[i].M_star<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].cm_star[k]<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].cmvel_star[k]<<" ";
            Fout<<pdata[i].Efrac_star<<" ";
            Fout<<pdata[i].Rhalfmass_star<<" ";
            for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<pdata[i].veldisp_star(k,n)<<" ";
            for (int k=0;k<3;k++) Fout<<pdata[i].L_star[k]<<" ";
            Fout<<pdata[i].q_star<<" ";
            Fout<<pdata[i].s_star<<" ";
            for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<pdata[i].eigvec_star(k,n)<<" ";
            Fout<<pdata[i].Krot_star<<" ";
            Fout<<pdata[i].t_star<<" ";
            Fout<<pdata[i].Z_star<<" ";
#endif
            Fout<<endl;
        }
    }
    cout<<"Done"<<endl;
    Fout.close();
}
//@}

///\name Writes the hierarchy of structures
void WriteHierarchy(Options &opt, const Int_t &ngroups, const Int_t & nhierarchy, const Int_t &nfield, Int_t *nsub, Int_t *parentgid, Int_t *stype, int subflag){
    fstream Fout;
    fstream Fout2;
    char fname[500],fname2[500];
    Int_t ng=0,noffset=0;
#ifdef USEMPI
    sprintf(fname,"%s.catalog_groups.%d",opt.outname,ThisTask);
#else
    int ThisTask=0,NProcs=1;
    sprintf(fname,"%s.catalog_groups",opt.outname);
#endif
    cout<<"saving hierarchy data to "<<fname<<endl;

    if (opt.ibinaryout) Fout.open(fname,ios::out|ios::binary|ios::app);
    else Fout.open(fname,ios::out|ios::app);

    //since the hierarchy file is appended to the catalog_groups files, no header written
#ifdef USEMPI
    for (int j=0;j<NProcs;j++) ng+=mpi_ngroups[j];
    for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
#else
    ng=ngroups;
#endif

    //if subflag==0 only write number of substructures
    if (subflag==0) {
    if (opt.ibinaryout) Fout.write((char*)&nsub[1],sizeof(Int_t)*nfield);
    else for (Int_t i=1;i<=nfield;i++)Fout<<nsub[i]<<endl;
    }
    else if (subflag==1) {
    if (opt.ibinaryout) {
        Fout.write((char*)&nsub[1+nfield],sizeof(Int_t)*(ngroups-nfield));
        Fout.write((char*)&parentgid[nfield+1],sizeof(Int_t)*(ngroups-nfield));
    }
    else {
        for (Int_t i=nfield+1;i<=ngroups;i++)Fout<<nsub[i]<<endl;
        for (Int_t i=nfield+1;i<=ngroups;i++)Fout<<parentgid[i]<<endl;
    }
    }
    //write everything, no distinction made between field and substructure
    else if (subflag==-1) {
    if (opt.ibinaryout) {
        Fout.write((char*)&nsub[1],sizeof(Int_t)*ngroups);
        Fout.write((char*)&parentgid[1],sizeof(Int_t)*ngroups);
    }
    else {
        for (Int_t i=1;i<=ngroups;i++)Fout<<nsub[i]<<endl;
        for (Int_t i=1;i<=ngroups;i++)Fout<<parentgid[i]<<endl;
    }
    }
    Fout.close();

    //now write a completely separate hierarchy file which I find more intuitive to parse
#ifdef USEMPI
    sprintf(fname,"%s.hierarchy.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.hierarchy",opt.outname);
#endif
    if (opt.ibinaryout) Fout.open(fname,ios::out|ios::binary);
    else Fout.open(fname,ios::out);

    cout<<"saving hierarchy data to "<<fname<<endl;
    if (opt.ibinaryout) {
    Fout.write((char*)&ThisTask,sizeof(int));
    Fout.write((char*)&NProcs,sizeof(int));
    Fout.write((char*)&ngroups,sizeof(Int_t));
    Fout.write((char*)&ng,sizeof(Int_t));
    }
    else {
    Fout<<ThisTask<<" "<<NProcs<<endl;
    Fout<<ngroups<<" "<<ng<<endl;
    Fout<<setprecision(10);
    }
    if (subflag==0) {
    if (opt.ibinaryout) {
        Fout.write((char*)&parentgid[1],sizeof(Int_t)*nfield);
        Fout.write((char*)&nsub[1],sizeof(Int_t)*nfield);
    }
    else for (Int_t i=1;i<=nfield;i++)Fout<<parentgid[i]<<" "<<nsub[i]<<endl;
    }
    else if (subflag==1) {
    if (opt.ibinaryout) {
        Fout.write((char*)&parentgid[1+nfield],sizeof(Int_t)*(ngroups-nfield));
        Fout.write((char*)&nsub[1+nfield],sizeof(Int_t)*(ngroups-nfield));
    }
    else for (Int_t i=nfield+1;i<=ngroups;i++)Fout<<parentgid[i]<<" "<<nsub[i]<<endl;
    }
    Fout.close();
    cout<<"Done saving hierarchy"<<endl;
}
//@}

/// \name Routines that can be used to output information of a halo subvolume decomposition
//@{
///Writes cell quantites
void WriteCellValues(Options &opt, const Int_t nbodies, const Int_t ngrid, GridCell *grid, Coordinate *gvel, Matrix *gveldisp)
{
    fstream Fout;
    char fname[1000];
#ifdef USEMPI
    if(opt.gname==NULL) sprintf(fname,"%s.griddata.%d",opt.outname,ThisTask);
    else sprintf(fname,"%s.%d",opt.gname,ThisTask);
#else
    if(opt.gname==NULL) sprintf(fname,"%s.griddata",opt.outname);
    else sprintf(fname,"%s",opt.gname);
#endif
    if (opt.ibinaryout) {
    Fout.open(fname,ios::out|ios::binary);
    Fout.write((char*)&nbodies,sizeof(Int_t));
    Fout.write((char*)&ngrid,sizeof(Int_t));
    for (Int_t i=0;i<ngrid;i++){
        Fout.write((char*)&grid[i].ndim,sizeof(Int_t));
        for (int j=0;j<grid[i].ndim;j++) Fout.write((char*)&grid[i].xm[j],sizeof(Double_t));
        for (int j=0;j<grid[i].ndim;j++) Fout.write((char*)&grid[i].xbl[j],sizeof(Double_t));
        for (int j=0;j<grid[i].ndim;j++) Fout.write((char*)&grid[i].xbu[j],sizeof(Double_t));
        Fout.write((char*)&grid[i].mass,sizeof(Double_t));
        Fout.write((char*)&grid[i].rsize,sizeof(Double_t));
        Fout.write((char*)&grid[i].nparts,sizeof(Int_t));
        Fout.write((char*)grid[i].nindex,sizeof(Int_t)*grid[i].nparts);
    }
    Fout.write((char*)gvel,sizeof(Coordinate)*ngrid);
    Fout.write((char*)gveldisp,sizeof(Matrix)*ngrid);
    }
    else {
    Fout.open(fname,ios::out);
    Fout<<nbodies<<" "<<ngrid<<endl;
    Fout<<scientific<<setprecision(10);
    for (Int_t i=0;i<ngrid;i++){
        Fout<<grid[i].ndim<<" ";
        for (int j=0;j<grid[i].ndim;j++) Fout<<grid[i].xm[j]<<" ";
        for (int j=0;j<grid[i].ndim;j++) Fout<<grid[i].xbl[j]<<" ";
        for (int j=0;j<grid[i].ndim;j++) Fout<<grid[i].xbu[j]<<" ";
        Fout<<grid[i].mass<<" "<<grid[i].rsize<<" "<<grid[i].nparts<<" ";
        for (int j=0;j<grid[i].nparts;j++) Fout<<grid[i].nindex[j]<<" ";
        Fout<<endl;
    }
    for (Int_t i=0;i<ngrid;i++){
        for (int j=0;j<3;j++)Fout<<gvel[i][j]<<" ";
        for (int j=0;j<3;j++)for (int k=0;k<3;k++)Fout<<gveldisp[i](j,k)<<" ";
        Fout<<endl;
    }
    }
    Fout.close();
}
//@}


/// \name Read group ids which can be useful if {\em Halo} have already been found
//@{
Int_t ReadPFOF(Options &opt, Int_t nbodies, Int_t *pfof){
    fstream Fin;
    Int_t temp;
    char fname[400];
	Int_t ngroup=0;
    sprintf(fname,"%s.fof.grp",opt.outname);
    cout<<"reading fof data "<<fname<<endl;
    Fin.open(fname,ios::in);
    Fin>>nbodies;
    //nbodies=opt.numpart[DARKTYPE];
    //for (Int_t i=0;i<opt.numpart[GASTYPE];i++) Fin>>temp;
	for (Int_t i=0;i<nbodies;i++) {Fin>>pfof[i];if (pfof[i]>ngroup) ngroup=pfof[i];}
    Fin.close();
    cout<<"Done"<<endl;
	return ngroup;
}

//load binary group fof catalogue 
Int_t ReadFOFGroupBinary(Options &opt, Int_t nbodies, Int_t *pfof, Int_t *idtoindex, Int_t minid, Particle *p)
{//old groupcat format
  char buf[1024];
  int TotNgroups,NFiles,dummy,Ngroups,Nids;
  int *pids;
  int *numingroup;
  fstream Fin;

  //group tab contains bulk info of groups
  sprintf(buf, "%s/group_tab_%03d", opt.gname, opt.snum);
  cout<<buf<<endl;
  Fin.open(buf,ios::in|ios::binary);

  //read group header info (number of groups, number of particles in all groups, total number of groups (in case split across several files)
  //and number of files
  Fin.read((char*)&Ngroups,sizeof(int));
  Fin.read((char*)&Nids,sizeof(int));
  Fin.read((char*)&TotNgroups,sizeof(int));
  Fin.read((char*)&NFiles,sizeof(int));
  cout<<Ngroups<<" fof groups in files "<<endl;
  cout<<Nids<<" particles in groups "<<endl;

  numingroup=new int[Ngroups];
//offsets are sum of lengths starting at 0
//ids of particles order according to group with first particle beloing to group 0, read n1
  Fin.read((char*)numingroup,sizeof(int)*Ngroups);
  Fin.close();

  //group ids contains actual particle ids in the group
  sprintf(buf, "%s/group_ids_%03d", opt.gname, opt.snum);
  Fin.open(buf,ios::in|ios::binary);

  //reread header info
  Fin.read((char*)&dummy,sizeof(int));
  Fin.read((char*)&dummy,sizeof(int));
  Fin.read((char*)&dummy,sizeof(int));
  Fin.read((char*)&dummy,sizeof(int));

  //read array of particle ids belong to groups
  pids=new int[Nids];
  Fin.read((char*)pids,sizeof(int)*Nids);
  int offset=0;
  for (Int_t i=0;i<Ngroups;i++) {
    for (Int_t j=0;j<numingroup[i];j++) {
        pfof[idtoindex[pids[j+offset]-minid]]=i+1;
    }
    offset+=numingroup[i];
  }
  Fin.close();
  return Ngroups;
}

//@}
