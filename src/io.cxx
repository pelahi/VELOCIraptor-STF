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
#ifdef USEXDR
#endif
#include "nchiladaitems.h"
#ifdef USEADIOS
#include "adios.h"
#endif

///write the information stored in a unit struct as meta data into a HDF5 file
#ifdef USEHDF
inline void WriteHeaderUnitEntry(Options & opt, H5OutputFile & hfile, string datasetname, HeaderUnitInfo &u)
{
    hfile.write_attribute(datasetname, "Dimension_Mass", u.massdim);
    hfile.write_attribute(datasetname, "Dimension_Length", u.lengthdim);
    hfile.write_attribute(datasetname, "Dimension_Velocity", u.velocitydim);
    hfile.write_attribute(datasetname, "Dimension_Time", u.timedim);
    if (u.extrainfo.size()>0){
        hfile.write_attribute(datasetname, "Dimension_Extra_Info", u.extrainfo);
    }
}
#endif

///Checks if file exits by attempting to get the file attributes
///If success file obviously exists.
///If failure may mean that we don't have permission to access the folder which contains this file or doesn't exist.
///If we need to do that level of checking, lookup return values of stat which will give you more details on why stat failed.
bool FileExists(const char *fname) {
  struct stat stFileInfo;
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
        struct tipsy_dump tipsyheader;
        fstream Ftip(opt.fname, ios::in | ios::binary);
        if (!Ftip){cerr<<"ERROR: Unable to open " <<opt.fname<<endl;exit(8);}
        else cout<<"Reading tipsy format from "<<opt.fname<<endl;
        Ftip.read((char*)&tipsyheader,sizeof(tipsy_dump));
        tipsyheader.SwitchtoBigEndian();
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
    else if (opt.inputtype==IOHDF) return HDF_get_nbodies(opt.fname,opt.partsearchtype,opt);
#endif
#ifdef USEXDR
    else if (opt.inputtype==IONCHILADA) return Nchilada_get_nbodies(opt.fname,opt.partsearchtype,opt);
#endif
    return 0;
}

///Reads particle data
///To add a new interface simply alter this to include the appropriate user written call
void ReadData(Options &opt, vector<Particle> &Part, const Int_t nbodies, Particle *&Pbaryons, Int_t nbaryons)
{
    InitEndian();
#ifndef USEMPI
    int ThisTask = 0, NProcs = 1;
    Int_t Nlocal = nbodies;
#endif
    if (ThisTask==0) {
        cout<<"Reading input ... "<<endl;
#ifdef USEMPI
        cout<<"Each MPI read thread, of which there are "<<opt.nsnapread<<", will allocate ";
        cout<<opt.mpiparticlebufsize*NProcs*sizeof(Particle)/1024.0/1024.0/1024.0<<" of memory to store particle data"<<endl;
        cout<<"Sending information to non-read threads in chunks of "<<opt.mpiparticlebufsize<<" particles "<<endl;
        cout<<"This requires approximately "<<(int)(Nlocal/(double)opt.mpiparticlebufsize)<<" receives"<<endl;
#endif
    }

    if(opt.inputtype==IOTIPSY) ReadTipsy(opt,Part,nbodies, Pbaryons, nbaryons);
    else if (opt.inputtype==IOGADGET) ReadGadget(opt,Part,nbodies, Pbaryons, nbaryons);
    else if (opt.inputtype==IORAMSES) ReadRamses(opt,Part,nbodies, Pbaryons, nbaryons);
#ifdef USEHDF
    else if (opt.inputtype==IOHDF) ReadHDF(opt,Part,nbodies, Pbaryons, nbaryons);
#endif
#ifdef USEXDR
    else if (opt.inputtype==IONCHILADA) ReadNchilada(opt,Part,nbodies, Pbaryons, nbaryons);
#endif
#ifdef NOMASS
    NOMASSCheck(opt);
#endif
    AdjustHydroQuantities(opt,Part,nbodies);
    AdjustStarQuantities(opt,Part,nbodies);
    AdjustBHQuantities(opt,Part,nbodies);
#ifdef USEMPI
    MPIAdjustDomain(opt);
#endif
    if (ThisTask==0) cout<<"Done loading input data"<<endl;
    GetMemUsage(opt,__func__+string("--line--")+to_string(__LINE__), (opt.iverbose>=1));
}


//Adjust particle data to appropriate units
void AdjustHydroQuantities(Options &opt, vector<Particle> &Part, const Int_t nbodies) {
    #ifdef GASON
    for (auto &p:Part) {
        if (p.GetType()!=GASTYPE) continue;
        p.SetU(p.GetU()*opt.internalenergyinputconversion);
    }
    #ifdef STARON
    if (opt.metallicityinputconversion!=1.0) {
        for (auto &p:Part) {
            if (p.GetType()!=GASTYPE) continue;
            p.SetZmet(p.GetZmet()*opt.metallicityinputconversion);
        }
    }
    if (opt.isfrisssfr==1) {
        for (auto &p:Part) {
            if (p.GetType()!=GASTYPE) continue;
            p.SetSFR(p.GetSFR()*p.GetMass());
        }
    }
    if (opt.SFRinputconversion!=1.0) {
        for (auto &p:Part) {
            if (p.GetType()!=GASTYPE) continue;
            p.SetSFR(p.GetSFR()*opt.SFRinputconversion);
        }
    }
    #endif
    #endif
}

void AdjustStarQuantities(Options &opt, vector<Particle> &Part, const Int_t nbodies) {
    #ifdef STARON
    if (opt.metallicityinputconversion!=1.0) {
        for (auto &p:Part) {
            if (p.GetType()!=STARTYPE) continue;
            p.SetZmet(p.GetZmet()*opt.metallicityinputconversion);
        }
    }
    if (opt.istellaragescalefactor!=0 || opt.stellarageinputconversion!=1.0) {
        double tage;
        for (auto &p:Part) {
            if (p.GetType()!=STARTYPE) continue;
            //if stellar age is initially stored as scale factor of formation
            if (opt.istellaragescalefactor == 1) {
		        //make sure units are consisten with those internal to the input tables as CalcCosmicTime returns time in yrs.
                tage = CalcCosmicTime(opt,p.GetTage(),opt.a) / opt.stellaragetoyrs;
            }
            //if stellar age is initially stored as redshift of formation
            else if (opt.istellaragescalefactor == 2) {
                tage = CalcCosmicTime(opt,1.0/(p.GetTage()+1),opt.a) / opt.stellaragetoyrs;
            }
            //if stellar age is initially stored as time of formation
            else if (opt.istellaragescalefactor == 3) {
                tage = opt.a-p.GetTage();
            }
            //if stellar age is initially stored as an age
            else tage = p.GetTage();
            tage*=opt.stellarageinputconversion;
            p.SetTage(tage);
        }
    }
    #endif
}

void AdjustBHQuantities(Options &opt, vector<Particle> &Part, const Int_t nbodies) {
    #ifdef BHON
    if (opt.metallicityinputconversion!=1.0) {
        for (auto &p:Part) {
            if (p.GetType()!=BHTYPE) continue;
            p.SetZmet(p.GetZmet()*opt.metallicityinputconversion);
        }
    }
    #endif
}
//@}

///\name Read STF data files
//@{

///Read local velocity density
void ReadLocalVelocityDensity(Options &opt, const Int_t nbodies, vector<Particle> &Part){
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
    if (opt.ibinaryout==OUTBINARY) {
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
void WriteLocalVelocityDensity(Options &opt, const Int_t nbodies, vector<Particle> &Part){
    fstream Fout;
    char fname[1000];
#ifdef USEMPI
    if(opt.smname==NULL) sprintf(fname,"%s.smdata.%d",opt.outname,ThisTask);
    else sprintf(fname,"%s.%d",opt.smname,ThisTask);
#else
    if(opt.smname==NULL) sprintf(fname,"%s.smdata",opt.outname);
    else sprintf(fname,"%s",opt.smname);
#endif
    if (opt.ibinaryout==OUTBINARY) {
        Fout.open(fname,ios::out|ios::binary);
        Fout.write((char*)&nbodies,sizeof(Int_t));
        Double_t tempd;
        for(Int_t i=0;i<nbodies;i++) {tempd=Part[i].GetDensity();Fout.write((char*)&tempd,sizeof(Double_t));}
    }
    if (opt.ibinaryout==OUTBINARY) {
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

void WriteGroupCatalog(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, vector<Particle> &Part, Int_t nadditional){
    fstream Fout,Fout2,Fout3;
    string fname, fname2, fname3;
    ostringstream os;
    unsigned long long ngtot=0,nids=0,nidstot=0,nuids=0, nuidstot=0, ng=0;
#ifdef USEPARALLELHDF
    unsigned long long nwritecommtot=0, nuwritecommtot=0;
#endif
    vector<unsigned long long> groupdata;
    vector<unsigned long long> offset;
    vector<long long> partdata;
#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif
#ifdef USEHDF
    H5OutputFile Fhdf, Fhdf3;
    int itemp=0;
#if defined(USEMPI) && defined(USEPARALLELHDF)
    vector<Int_t> mpi_ngoffset(NProcs);
    Int_t ngoffset;
#endif
#endif
#ifdef USEADIOS
    int adios_err;
    uint64_t adios_groupsize , adios_totalsize ;
    int64_t adios_file_handle,adios_file_handle3;
    int64_t adios_grp_handle, adios_grp_handle3;
    int64_t adios_var_handle;
    int64_t adios_attr_handle;
#endif
#if defined(USEHDF)||defined(USEADIOS)
    DataGroupNames datagroupnames;
#endif

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    os << opt.outname << ".catalog_groups";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
#endif
    fname = os.str();

    cout<<"saving group catalog to "<<fname<<endl;
    if (opt.ibinaryout==OUTBINARY) Fout.open(fname,ios::out|ios::binary);
#ifdef USEHDF
#ifdef USEPARALLELHDF
    else if (opt.ibinaryout==OUTHDF) {
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
        }
        else{
            Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
        }
    }
#else
    //create file
    else if (opt.ibinaryout==OUTHDF) {
        Fhdf.create(string(fname));
    }
#endif
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //open an adios file
        adios_err=adios_open(&adios_file_handle, "VELOCIraptor_catalog_groups", fname, "w", MPI_COMM_WORLD);
    }
#endif
    else Fout.open(fname,ios::out);
    ng=ngroups;

#ifdef USEMPI
    if (NProcs > 1) {
        for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
    }
    else {
        ngtot = ngroups;
    }

#else
    ngtot=ngroups+nadditional;//useful if outputing field halos
#endif
    //write header
    if (opt.ibinaryout==OUTBINARY) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&ng,sizeof(unsigned long));
        Fout.write((char*)&ngtot,sizeof(unsigned long));
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            MPI_Allreduce(&ng, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            if (ThisWriteTask==0) {
                Fhdf.write_dataset(opt, datagroupnames.group[itemp++], 1, &ThisWriteComm, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.group[itemp++], 1, &NWriteComms, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.group[itemp++], 1, &nwritecommtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.group[itemp++], 1, &ngtot, -1, -1, false);
            }
            else {
                itemp=4;
            }
            Fhdf.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
        }
        else{
            itemp=0;
            Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &ThisTask);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &NProcs);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &ng);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &ngtot);
            itemp++;
        }
#else
        itemp=0;
        Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &ThisTask);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &NProcs);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &ng);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.group[itemp], 1, &ngtot);
        itemp++;
#endif
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //declare the attributes in a header group, assiging the group handle, setting the name, no time step indicator, and a flag saying yes to all statistics
        adios_err=adios_declare_group(&adios_grp_handle,"Header", "" , adios_stat_full);
        //select simple mpi method
        adios_select_method (adios_grp_handle, "MPI", "", "");
        //define some attributes
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(ThisTask).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(NProcs).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(ng).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(ngtot).c_str(),"");
        itemp++;
        ///\todo don't actually know if I should use adios attribute or var to store simple single values
    }
#endif
    else{
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ng<<" "<<ngtot<<endl;
    }

    //write group size
    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&numingroup[1],sizeof(Int_t)*ngroups);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        groupdata.resize(ng+1,0);
        for (Int_t i=1;i<=ng;i++) groupdata[i-1]=numingroup[i];
        Fhdf.write_dataset(opt, datagroupnames.group[itemp], ng, groupdata.data());
        itemp++;
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //declare a new group
        //declare the attributes in a header group, assiging the group handle, setting the name, no time step indicator, and a flag saying yes to all statistics
        adios_err=adios_declare_group(&adios_grp_handle,"Catalog_Data", "" , adios_stat_full);
        //select simple mpi method
        adios_select_method (adios_grp_handle, "MPI", "", "");
        //now stage variables (data)
        //if want to define dimensions can either create a variable that stores the dimensions or store the value as a string.
        //store local dim
        adios_err=adios_define_var(adios_grp_handle,"ng","", adios_unsigned_long,0,0,0);
        //store global dim
        adios_err=adios_define_var(adios_grp_handle,"ngtot","", adios_unsigned_long,0,0,0);
        //store mpi offset
        adios_err=adios_define_var(adios_grp_handle,"ngmpioffset","", adios_unsigned_long,0,0,0);
        //then define the group actually storing the data. Might be useful to define an offset variable as well for quick access when reading
        //offset would be the last field in the code below
        adios_err=adios_define_var(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],"ng","ngtot","ngmpioffset");
        adios_err=adios_write(adios_file_handle,"ng",&ng);
        adios_err=adios_write(adios_file_handle,"ngtot",&ngtot);
        Int_t mpioffset=0;
        for (Int_t itask=0;itask<ThisTask;itask++)mpioffset+=mpi_ngroups[itask];
        adios_err=adios_write(adios_file_handle,"ngmpioffset",&mpioffset);
        groupdata.resize(ng+1);
        for (Int_t i=1;i<=ng;i++) groupdata[i-1]=numingroup[i];
        adios_err=adios_write(adios_file_handle,datagroupnames.group[itemp].c_str(),groupdata.data());
        itemp++;
    }
#endif
    else for (Int_t i=1;i<=ngroups;i++) Fout<<numingroup[i]<<endl;

    //Write offsets for bound and unbound particles
    offset.resize(ng+2,0);
    //note before had offsets at numingroup but to account for unbound particles use value of pglist at numingroup
    if (ngroups >1) for (Int_t i=2;i<=ng;i++) offset[i]=offset[i-1]+pglist[i-1][numingroup[i-1]];
    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&(offset.data())[1],sizeof(Int_t)*ngroups);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        for (auto &x:groupdata) x=0;
        if (ng > 1) for (Int_t i=2;i<=ng;i++) groupdata[i-1]=offset[i];
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            nids = 0; for (Int_t i=1; i<=ng; i++) nids+=pglist[i][numingroup[i]];
            MPI_Allgather(&nids, 1, MPI_Int_t, mpi_ngoffset.data(), 1, MPI_Int_t, mpi_comm_write);
            ngoffset = 0;
            if (ThisWriteTask > 0)
            {
                for (auto itask = 0; itask < ThisWriteTask; itask++) ngoffset += mpi_ngoffset[itask];
                if (ng > 0) for (Int_t i=1; i<=ng; i++) groupdata[i-1] += ngoffset;
            }
        }
#endif
        Fhdf.write_dataset(opt, datagroupnames.group[itemp], ng, groupdata.data());
        itemp++;
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //don't delcare new group, just add data
        adios_err=adios_define_var(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],"ng","ngtot","ngmpioffset");
        for (Int_t i=1;i<=ng;i++) groupdata[i-1]=offset[i];
        adios_err=adios_write(adios_file_handle,datagroupnames.group[itemp].c_str(),groupdata.data());
        itemp++;
    }
#endif
    else for (Int_t i=1;i<=ng;i++) Fout<<offset[i]<<endl;

    //position of unbound particle
    for (auto &x:offset) x=0;
    if (ng >1) for (Int_t i=2;i<=ng;i++) offset[i]=offset[i-1]+numingroup[i-1]-pglist[i-1][numingroup[i-1]];
    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&(offset.data())[1],sizeof(Int_t)*ngroups);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        for (auto &x:groupdata) x=0;
        if (ng > 1) for (Int_t i=2;i<=ng;i++) groupdata[i-1]=offset[i];
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            nuids = 0; for (Int_t i=1; i<=ng; i++) nuids+=numingroup[i]-pglist[i][numingroup[i]];
            MPI_Allgather(&nuids, 1, MPI_Int_t, mpi_ngoffset.data(), 1, MPI_Int_t, mpi_comm_write);
            ngoffset = 0;
            if (ThisWriteTask > 0)
            {
                for (auto itask = 0; itask < ThisWriteTask; itask++) ngoffset += mpi_ngoffset[itask];
                if (ng > 0) for (Int_t i=1; i<=ng; i++) groupdata[i-1] += ngoffset;
            }
        }
#endif
        Fhdf.write_dataset(opt, datagroupnames.group[itemp], ng, groupdata.data());
        itemp++;
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //don't delcare new group, just add data
        adios_err=adios_define_var(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],"ng","ngtot","ngmpioffset");
        groupdata.resize(ng+1,0);
        for (Int_t i=1;i<=ng;i++) groupdata[i-1]=offset[i];
        adios_err=adios_write(adios_file_handle,datagroupnames.group[itemp].c_str(),groupdata.data());
        itemp++;
    }
#endif
    else for (Int_t i=1;i<=ng;i++) Fout<<offset[i]<<endl;
    offset.clear();

    if (opt.ibinaryout==OUTASCII || opt.ibinaryout==OUTBINARY) Fout.close();
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) Fhdf.close();
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) adios_err=adios_close(adios_file_handle);
#endif

    //now write pid files
    os.str(string());
    os <<opt.outname<< ".catalog_particles";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
#endif
    fname = os.str();
    os.str(string());
    os <<opt.outname<< ".catalog_particles.unbound";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
#endif
    fname3 = os.str();

    cout<<"saving particle catalog to "<<fname<<endl;

    if (opt.ibinaryout==OUTBINARY) {
        Fout.open(fname,ios::out|ios::binary);
        Fout3.open(fname3,ios::out|ios::binary);
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
            Fhdf3.create(string(fname3),H5F_ACC_TRUNC, 0, false);
        }
        else{
            Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
            Fhdf3.create(string(fname3),H5F_ACC_TRUNC, ThisWriteComm, false);
        }
#else
        Fhdf.create(string(fname));
        Fhdf3.create(string(fname3));
#endif
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        adios_err=adios_open(&adios_file_handle, "VELOCIraptor_catalog_particles", fname, "w", MPI_COMM_WORLD);
        adios_err=adios_open(&adios_file_handle3, "VELOCIraptor_catalog_particles.unbound", fname3, "w", MPI_COMM_WORLD);
    }
#endif
    else {
        Fout.open(fname,ios::out);
        Fout3.open(fname3,ios::out);
    }

    //see above regarding unbound particle
    //for (Int_t i=1;i<=ngroups;i++) nids+=numingroup[i];
    nids = nuids = 0;
    for (Int_t i=1;i<=ngroups;i++) {nids+=pglist[i][numingroup[i]];nuids+=numingroup[i]-pglist[i][numingroup[i]];}
#ifdef USEMPI
    if (NProcs > 1) {
        MPI_Allreduce(&nids, &nidstot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nuids, &nuidstot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
        nidstot=nids;
        nuidstot=nuids;
    }
#else
    nidstot=nids;
    nuidstot=nuids;
#endif

    //write header
    if (opt.ibinaryout==OUTBINARY) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&nids,sizeof(unsigned long));
        Fout.write((char*)&nidstot,sizeof(unsigned long));

        Fout3.write((char*)&ThisTask,sizeof(int));
        Fout3.write((char*)&NProcs,sizeof(int));
        Fout3.write((char*)&nuids,sizeof(unsigned long));
        Fout3.write((char*)&nuidstot,sizeof(unsigned long));
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            MPI_Allreduce(&nids, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            MPI_Allreduce(&nuids, &nuwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            if (ThisWriteTask == 0) {
                itemp=0;
                Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &ThisWriteComm, -1, -1, false);
                Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &ThisWriteComm, -1, -1, false);
                itemp++;
                Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &NWriteComms, -1, -1, false);
                Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &NWriteComms, -1, -1, false);
                itemp++;
                Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &nwritecommtot, -1, -1, false);
                Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &nuwritecommtot, -1, -1, false);
                itemp++;
                Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &nidstot, -1, -1, false);
                Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &nuidstot, -1, -1, false);
                itemp++;
            }
            else {
                itemp=4;
            }
            Fhdf.close();
            Fhdf3.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
            Fhdf3.append(string(fname3));
        }
        else{
            itemp=0;
            Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &ThisTask);
            Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &ThisTask);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &NProcs);
            Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &NProcs);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &nids);
            Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &nuids);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &nidstot);
            Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &nuidstot);
            itemp++;
        }
#else
        itemp=0;
        Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &ThisTask);
        Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &ThisTask);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &NProcs);
        Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &NProcs);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &nids);
        Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &nuids);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.part[itemp], 1, &nidstot);
        Fhdf3.write_dataset(opt, datagroupnames.part[itemp], 1, &nuidstot);
        itemp++;
#endif

    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //declare the attributes in a header group, assiging the group handle, setting the name, no time step indicator, and a flag saying yes to all statistics
        adios_err=adios_declare_group(&adios_grp_handle,"Header", "" , adios_stat_full);
        //select simple mpi method
        adios_select_method (adios_grp_handle, "MPI", "", "");
        //define some attributes
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(ThisTask).c_str(),"");
        adios_err=adios_define_attribute(adios_grp_handle3,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(ThisTask).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(NProcs).c_str(),"");
        adios_err=adios_define_attribute(adios_grp_handle3,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(NProcs).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(nids).c_str(),"");
        adios_err=adios_define_attribute(adios_grp_handle3,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(nuids).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(nidstot).c_str(),"");
        adios_err=adios_define_attribute(adios_grp_handle3,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],to_string(nuidstot).c_str(),"");
        itemp++;
        ///\todo don't actually know if I should use adios attribute or var to store simple single values
    }
#endif
    else {
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<nids<<" "<<nidstot<<endl;

        Fout3<<ThisTask<<" "<<NProcs<<endl;
        Fout3<<nuids<<" "<<nuidstot<<endl;
    }

    Int_t *idval;
    idval=new Int_t[nids+1];
    nids=0;
    for (Int_t i=1;i<=ngroups;i++)
        for (Int_t j=0;j<pglist[i][numingroup[i]];j++)
            idval[nids++]=Part[pglist[i][j]].GetPID();
    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)idval,sizeof(Int_t)*nids);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        partdata.resize(nids+1);
        for (Int_t i=0;i<nids;i++) partdata[i]=idval[i];
        Fhdf.write_dataset(opt, datagroupnames.part[itemp], nids, partdata.data());
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        adios_err=adios_declare_group(&adios_grp_handle,"Catalog_Data", "" , adios_stat_full);
        adios_select_method (adios_grp_handle, "MPI", "", "");
        //store local dim
        adios_err=adios_define_var(adios_grp_handle,"nids","", adios_unsigned_long,0,0,0);
        //store global dim
        adios_err=adios_define_var(adios_grp_handle,"nidstot","", adios_unsigned_long,0,0,0);
        //store mpi offset
        adios_err=adios_define_var(adios_grp_handle,"nidsmpioffset","", adios_unsigned_long,0,0,0);
        adios_err=adios_define_var(adios_grp_handle,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],"nids","nidstot","nidsmpioffset");
        adios_err=adios_write(adios_file_handle,"nids",&nids);
        adios_err=adios_write(adios_file_handle,"nidstot",&nidstot);
        Int_t mpioffset=0;
        //for (Int_t itask=0;itask<ThisTask;itask++)mpioffset+=mpi_ngroups[itask];
        adios_err=adios_write(adios_file_handle,"nidsmpioffset",&mpioffset);
        if (nids > 0) {
            partdata.resize(nids+1);
            for (Int_t i=0;i<nids;i++) partdata[i-1]=idval[i];
            adios_err=adios_write(adios_file_handle,datagroupnames.group[itemp].c_str(),partdata.data());
        }
    }
#endif
    else for (Int_t i=0;i<nids;i++) Fout<<idval[i]<<endl;
    delete[] idval;
    if (opt.ibinaryout==OUTASCII || opt.ibinaryout==OUTBINARY) Fout.close();
#ifdef USEHDF
    if (opt.ibinaryout==OUTHDF) Fhdf.close();
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) adios_err=adios_close(adios_file_handle);
#endif

    idval=new Int_t[nuids+1];
    nuids=0;
    for (Int_t i=1;i<=ngroups;i++)
        for (Int_t j=pglist[i][numingroup[i]];j<numingroup[i];j++)
            idval[nuids++]=Part[pglist[i][j]].GetPID();
    if (opt.ibinaryout==OUTBINARY) Fout3.write((char*)idval,sizeof(Int_t)*nuids);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        partdata.resize(nuids+1);
        for (Int_t i=0;i<nuids;i++) partdata[i]=idval[i];
        Fhdf3.write_dataset(opt, datagroupnames.part[itemp], nuids, partdata.data());
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        adios_err=adios_declare_group(&adios_grp_handle3,"Catalog_Data", "" , adios_stat_full);
        adios_select_method (adios_grp_handle3, "MPI", "", "");
        //store local dim
        adios_err=adios_define_var(adios_grp_handle3,"nuids","", adios_unsigned_long,0,0,0);
        //store global dim
        adios_err=adios_define_var(adios_grp_handle3,"nuidstot","", adios_unsigned_long,0,0,0);
        //store mpi offset
        adios_err=adios_define_var(adios_grp_handle3,"nuidsmpioffset","", adios_unsigned_long,0,0,0);
        adios_err=adios_define_var(adios_grp_handle3,datagroupnames.group[itemp].c_str(),"",datagroupnames.adiosgroupdatatype[itemp],"nuids","nuidstot","nuidsmpioffset");
        adios_err=adios_write(adios_file_handle3,"nuids",&nuids);
        adios_err=adios_write(adios_file_handle3,"nuidstot",&nuidstot);
        Int_t mpioffset=0;
        //for (Int_t itask=0;itask<ThisTask;itask++)mpioffset+=mpi_ngroups[itask];
        adios_err=adios_write(adios_file_handle3,"nidsmpioffset",&mpioffset);
        if (nuids > 0) {
            partdata.resize(nuids+1);
            for (Int_t i=0;i<nuids;i++) partdata[i-1]=idval[i];
            adios_err=adios_write(adios_file_handle3,datagroupnames.group[itemp].c_str(),partdata.data());
            delete[] data;
        }
    }
#endif
    else for (Int_t i=0;i<nuids;i++) Fout3<<idval[i]<<endl;
    delete[] idval;

    if (opt.ibinaryout==OUTASCII || opt.ibinaryout==OUTBINARY) Fout3.close();
#ifdef USEHDF
    if (opt.ibinaryout==OUTHDF) Fhdf3.close();
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) adios_err=adios_close(adios_file_handle3);
#endif

#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif
}

///if particles are separately searched (i.e. \ref Options.iBaryonSearch is set) then produce list of particle types
void WriteGroupPartType(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, vector<Particle> &Part){
    fstream Fout,Fout2;
    string fname, fname2;
    ostringstream os, os2;
    unsigned long long nids=0, nidstot=0, nuids=0, nuidstot=0;
#ifdef USEPARALLELHDF
    unsigned long long nwritecommtot=0, nuwritecommtot=0;
#endif
    int *typeval;

#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif
#ifdef USEHDF
    H5OutputFile Fhdf,Fhdf2;
    int itemp;
#endif
#if defined(USEHDF)||defined(USEADIOS)
    DataGroupNames datagroupnames;
#endif
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    os << opt.outname << ".catalog_parttypes";
    os2 << opt.outname << ".catalog_parttypes.unbound";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
        os2 << "." << ThisWriteComm;
#else
        os<<"."<<ThisTask;
        os2 << "." << ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
        os2 << "." << ThisTask;
    }
#endif
    fname = os.str();
    fname2 = os2.str();
    cout<<"saving particle type info to "<<fname<<endl;


    if (opt.ibinaryout==OUTBINARY) {
        Fout.open(fname,ios::out|ios::binary);
        Fout2.open(fname2,ios::out|ios::binary);
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
            Fhdf2.create(string(fname2),H5F_ACC_TRUNC, 0, false);
        }
        else{
            Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
            Fhdf2.create(string(fname2),H5F_ACC_TRUNC, ThisWriteComm, false);
        }
#else
        //create file
        Fhdf.create(string(fname));
        Fhdf2.create(string(fname2));
#endif
    }
#endif
    else {
        Fout.open(fname,ios::out);
        Fout2.open(fname2,ios::out);
    }

    for (Int_t i=1;i<=ngroups;i++) {nids+=pglist[i][numingroup[i]];nuids+=numingroup[i]-pglist[i][numingroup[i]];}
#ifdef USEMPI
    if(NProcs > 1) {
        MPI_Allreduce(&nids, &nidstot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&nuids, &nuidstot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
        nidstot=nids;
        nuidstot=nuids;
    }
#else
    nidstot=nids;
    nuidstot=nuids;
#endif

    //write header
    if (opt.ibinaryout==OUTBINARY) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&nids,sizeof(Int_t));
        Fout.write((char*)&nidstot,sizeof(Int_t));

        Fout2.write((char*)&ThisTask,sizeof(int));
        Fout2.write((char*)&NProcs,sizeof(int));
        Fout2.write((char*)&nuids,sizeof(Int_t));
        Fout2.write((char*)&nuidstot,sizeof(Int_t));
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            MPI_Allreduce(&nids, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            MPI_Allreduce(&nuids, &nuwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            if (ThisWriteTask == 0) {
                itemp=0;
                Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &ThisWriteComm, -1, -1, false);
                Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &ThisWriteComm, -1, -1, false);
                itemp++;
                Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &NWriteComms, -1, -1, false);
                Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &NWriteComms, -1, -1, false);
                itemp++;
                Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &nwritecommtot, -1, -1, false);
                Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &nuwritecommtot, -1, -1, false);
                itemp++;
                Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &nidstot, -1, -1, false);
                Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &nuidstot, -1, -1, false);
                itemp++;
            }
            else {
                itemp=4;
            }
            Fhdf.close();
            Fhdf2.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
            Fhdf2.append(string(fname2));
        }
        else{
            itemp=0;
            Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &ThisTask);
            Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &ThisTask);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &NProcs);
            Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &NProcs);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &nids);
            Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &nuids);
            itemp++;
            Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &nidstot);
            Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &nuidstot);
            itemp++;
        }
#else
        itemp=0;
        Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &ThisTask);
        Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &ThisTask);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &NProcs);
        Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &NProcs);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &nids);
        Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &nuids);
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.types[itemp], 1, &nidstot);
        Fhdf2.write_dataset(opt, datagroupnames.types[itemp], 1, &nuidstot);
        itemp++;
#endif
    }
#endif
    else {
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<nids<<" "<<nidstot<<endl;

        Fout2<<ThisTask<<" "<<NProcs<<endl;
        Fout2<<nuids<<" "<<nuidstot<<endl;
    }

    typeval=new int[nids+1];
    nids=0;
    for (Int_t i=1;i<=ngroups;i++)
        for (Int_t j=0;j<pglist[i][numingroup[i]];j++)
            typeval[nids++]=Part[pglist[i][j]].GetType();
    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)typeval,sizeof(int)*nids);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        unsigned short *data=new unsigned short[nids];
        for (Int_t i=0;i<nids;i++) data[i]=typeval[i];
        Fhdf.write_dataset(opt, datagroupnames.types[itemp], nids, data);
        delete[] data;
    }
#endif
    else for (Int_t i=0;i<nids;i++) Fout<<typeval[i]<<endl;
    delete[] typeval;
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

    typeval=new int[nuids];
    nuids=0;
    for (Int_t i=1;i<=ngroups;i++)
        for (Int_t j=pglist[i][numingroup[i]];j<numingroup[i];j++)
            typeval[nuids++]=Part[pglist[i][j]].GetType();
    if (opt.ibinaryout==OUTBINARY) Fout2.write((char*)typeval,sizeof(int)*nuids);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        unsigned short *data=new unsigned short[nuids];
        for (Int_t i=0;i<nuids;i++) data[i]=typeval[i];
        Fhdf2.write_dataset(opt, datagroupnames.types[itemp], nuids, data);
        delete[] data;
    }
#endif
    else for (Int_t i=0;i<nuids;i++) Fout2<<typeval[i]<<endl;
    delete[] typeval;
    if (opt.ibinaryout!=OUTHDF) Fout2.close();
#ifdef USEHDF
    else Fhdf2.close();
#endif

#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif
}

///Write the particles in each SO region
///Note that this particle list will not be exclusive
///\todo optimisation memory wise can be implemented by not creating an array
///to store all ids and then copying info from the array of vectors into it.
void WriteSOCatalog(Options &opt, const Int_t ngroups, vector<Int_t> *SOpids, vector<int> *SOtypes){
    fstream Fout;
    string fname;
    ostringstream os;
    unsigned long ng=0,ngtot=0,nSOids=0,nSOidstot=0;
#ifdef USEPARALLELHDF
    unsigned long long nwritecommtot=0, nSOwritecommtot=0;
#endif
    vector<unsigned long> offset;
    vector<long long> idval;
    vector<int> typeval;
    vector<Int_t> numingroup;

#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif
#ifdef USEHDF
    H5OutputFile Fhdf;
    int itemp=0;
#if defined(USEMPI) && defined(USEPARALLELHDF)
    vector<Int_t> mpi_offset(NProcs);
    Int_t nSOidoffset;
#endif
#endif
#ifdef USEADIOS
    int adios_err;
    uint64_t adios_groupsize , adios_totalsize ;
    int64_t adios_file_handle;
    int64_t adios_grp_handle;
    int64_t adios_var_handle;
    int64_t adios_attr_handle;
#endif
#if defined(USEHDF)||defined(USEADIOS)
    DataGroupNames datagroupnames;
#endif

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    ng=ngroups;
#ifdef USEMPI
    if (NProcs >1) {
        MPI_Allreduce(&ng, &ngtot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
        ngtot = ng;
    }
#else
    ngtot=ng;
#endif
    for (Int_t i=1;i<=ngroups;i++) nSOids+=SOpids[i].size();
#ifdef USEMPI
    if (NProcs > 1) {
        MPI_Allreduce(&nSOids, &nSOidstot, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
        nSOidstot = nSOids;
    }
#else
    nSOidstot = nSOids;
#endif

    os << opt.outname <<".catalog_SOlist";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
#endif
    fname = os.str();

    if (opt.iverbose) cout<<"saving SO particle lists to "<<fname<<endl;
    if (opt.ibinaryout==OUTBINARY) Fout.open(fname,ios::out|ios::binary);
#ifdef USEHDF
    //create file
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
        }
        else{
            Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
        }
#else
        Fhdf.create(string(fname));
#endif
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //open an adios file
        adios_err=adios_open(&adios_file_handle, "VELOCIraptor_SOlist", fname, "w", MPI_COMM_WORLD);
    }
#endif
    else Fout.open(fname,ios::out);

    //write header
    if (opt.ibinaryout==OUTBINARY) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&ng,sizeof(unsigned long));
        Fout.write((char*)&ngtot,sizeof(unsigned long));
        Fout.write((char*)&nSOids,sizeof(unsigned long));
        Fout.write((char*)&nSOidstot,sizeof(unsigned long));
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            MPI_Allreduce(&ng, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            MPI_Allreduce(&nSOids, &nSOwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            if (ThisWriteTask == 0) {
                itemp=0;
                Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ThisWriteComm, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &NWriteComms, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &nwritecommtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ngtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &nSOwritecommtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &nSOidstot, -1, -1, false);
            }
            else {
                itemp=6;
            }
            Fhdf.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
        }
        else{
            itemp=0;
            Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ThisTask);
            Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &NProcs);
            Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ng);
            Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ngtot);
            Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &nSOids);
            Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &nSOidstot);
        }
#else
        itemp=0;
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ThisTask);
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &NProcs);
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ng);
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &ngtot);
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &nSOids);
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp++], 1, &nSOidstot);
#endif
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS)
    {
        //declare the attributes in a header group, assiging the group handle, setting the name, no time step indicator, and a flag saying yes to all statistics
        adios_err=adios_declare_group(&adios_grp_handle,"Header", "" , adios_stat_full);
        //select simple mpi method
        adios_select_method (adios_grp_handle, "MPI", "", "");
        //define some attributes
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],to_string(ThisTask).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],to_string(NProcs).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],to_string(ng).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],to_string(ngtot).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],to_string(nSOids).c_str(),"");
        itemp++;
        adios_err=adios_define_attribute(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],to_string(nSOidstot).c_str(),"");
        itemp++;
        ///\todo don't actually know if I should use adios attribute or var to store simple single values
    }
#endif
    else
    {
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ng<<" "<<ngtot<<endl;
        Fout<<nSOids<<" "<<nSOidstot<<endl;
    }

    //write group size
    if (opt.ibinaryout==OUTBINARY) {
        numingroup.resize(ngroups+1,0);
        for (auto i=1;i<=ngroups;i++) numingroup[i] = SOpids[i].size();
        if (ngroups > 0 ) Fout.write((char*)&(numingroup.data())[1],sizeof(Int_t)*ngroups);
        numingroup.clear();
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        unsigned int *data=NULL;
        if (ng > 0) {
	        data=new unsigned int[ng];
        	for (Int_t i=1;i<=ng;i++) data[i-1]=SOpids[i].size();
        }
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp], ng, data);
        itemp++;
        delete[] data;
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //declare a new group
        //declare the attributes in a header group, assiging the group handle, setting the name, no time step indicator, and a flag saying yes to all statistics
        adios_err=adios_declare_group(&adios_grp_handle,"Catalog_Data", "" , adios_stat_full);
        //select simple mpi method
        adios_select_method (adios_grp_handle, "MPI", "", "");
        //now stage variables (data)
        //if want to define dimensions can either create a variable that stores the dimensions or store the value as a string.
        //store local dim
        adios_err=adios_define_var(adios_grp_handle,"ng","", adios_unsigned_long,0,0,0);
        //store global dim
        adios_err=adios_define_var(adios_grp_handle,"ngtot","", adios_unsigned_long,0,0,0);
        //store mpi offset
        adios_err=adios_define_var(adios_grp_handle,"ngmpioffset","", adios_unsigned_long,0,0,0);
        //then define the group actually storing the data. Might be useful to define an offset variable as well for quick access when reading
        //offset would be the last field in the code below
        adios_err=adios_define_var(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],"ng","ngtot","ngmpioffset");
        adios_err=adios_write(adios_file_handle,"ng",&ng);
        adios_err=adios_write(adios_file_handle,"ngtot",&ngtot);
        Int_t mpioffset=0;
        for (Int_t itask=0;itask<ThisTask;itask++)mpioffset+=mpi_ngroups[itask];
        adios_err=adios_write(adios_file_handle,"ngmpioffset",&mpioffset);
        unsigned int *data = NULL;
        if (ng > 0){
            data = new unsigned int[ng];
            for (Int_t i=1;i<=ng;i++) data[i-1]=SOpids[i].size();
        }
        adios_err=adios_write(adios_file_handle,datagroupnames.SO[itemp].c_str(),data);
        delete[] data;
        itemp++;
    }
#endif
    else {
        for (Int_t i=1;i<=ngroups;i++) Fout<<SOpids[i].size()<<endl;
    }


    //Write offsets
    offset.resize(ngroups+1,0);
    if (ngroups > 0) {
       offset[0] = offset[1] = 0;
       for (Int_t i=2;i<=ngroups;i++) offset[i]=offset[i-1]+SOpids[i].size();
    }

    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&(offset.data())[1],sizeof(Int_t)*ngroups);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        unsigned long *data = NULL;
        if (ng > 0) {
            data = new unsigned long[ng];
            for (Int_t i=1;i<=ng;i++) data[i-1]=offset[i];
        }
#ifdef USEPARALLELHDF
    if(opt.mpinprocswritesize>1){
            MPI_Allgather(&nSOids, 1, MPI_Int_t, mpi_offset.data(), 1, MPI_Int_t, mpi_comm_write);
            if (ThisWriteTask > 0)
            {
                nSOidoffset = 0; for (auto itask = 0; itask < ThisWriteTask; itask++) nSOidoffset += mpi_offset[itask];
                for (Int_t i=1; i<=ng; i++) data[i-1] += nSOidoffset;
            }
    }
#endif
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp], ng, data);
        itemp++;
        delete[] data;
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        //don't delcare new group, just add data
        adios_err=adios_define_var(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],"ng","ngtot","ngmpioffset");
        unsigned long *data = NULL;
        if (ng > 0) {
            data = new unsigned long[ng];
            for (Int_t i=1;i<=ng;i++) data[i-1]=offset[i];
        }
        adios_err=adios_write(adios_file_handle,datagroupnames.SO[itemp].c_str(),data);
        delete[] data;
        itemp++;
    }
#endif
    else {
        for (Int_t i=1;i<=ngroups;i++) Fout<<offset[i]<<endl;
    }
    offset.clear();

    if (nSOids>0) {
        idval.resize(nSOids,0);
        nSOids=0;
        for (Int_t i=1;i<=ngroups;i++) {
            for (Int_t j=0;j<SOpids[i].size();j++) {
                idval[nSOids++]=SOpids[i][j];
            }
            SOpids[i].clear();
        }
#if defined(GASON) || defined(STARON) || defined(BHON)
        typeval.resize(nSOids,0);
        nSOids=0;
        for (Int_t i=1;i<=ngroups;i++) {
            for (Int_t j=0;j<SOtypes[i].size();j++) {
                typeval[nSOids++]=SOtypes[i][j];
            }
            SOtypes[i].clear();
        }
#endif
    }
    if (opt.ibinaryout==OUTBINARY) {
        if (nSOids>0) Fout.write((char*)idval.data(),sizeof(Int_t)*nSOids);
#if defined(GASON) || defined(STARON) || defined(BHON)
        if (nSOids>0) Fout.write((char*)typeval.data(),sizeof(int)*nSOids);
#endif
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp], nSOids, idval.data());
        #if defined(GASON) || defined(STARON) || defined(BHON)
        itemp++;
        Fhdf.write_dataset(opt, datagroupnames.SO[itemp], nSOids, typeval.data());
        #endif
    }
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) {
        adios_err=adios_declare_group(&adios_grp_handle,"Particle_Data", "" , adios_stat_full);
        adios_select_method (adios_grp_handle, "MPI", "", "");
        //store local dim
        adios_err=adios_define_var(adios_grp_handle,"nSOids","", adios_unsigned_long,0,0,0);
        //store global dim
        adios_err=adios_define_var(adios_grp_handle,"nSOidstot","", adios_unsigned_long,0,0,0);
        //store mpi offset
        adios_err=adios_define_var(adios_grp_handle,"nSOidsmpioffset","", adios_unsigned_long,0,0,0);
        adios_err=adios_define_var(adios_grp_handle,datagroupnames.SO[itemp].c_str(),"",datagroupnames.adiosSOdatatype[itemp],"nSOids","nSOidstot","nSOidsmpioffset");
        adios_err=adios_write(adios_file_handle,"nSOids",&nSOids);
        adios_err=adios_write(adios_file_handle,"nSOidstot",&nSOidstot);
        Int_t mpioffset=0;
        adios_err=adios_write(adios_file_handle,"nSOidsmpioffset",&mpioffset);
        adios_err=adios_write(adios_file_handle,datagroupnames.SO[itemp].c_str(),idval);
#if defined(GASON) || defined(STARON) || defined(BHON)
        itemp++;
        adios_err=adios_write(adios_file_handle,datagroupnames.SO[itemp].c_str(),typeval);
#endif
    }
#endif
    else {
        for (Int_t i=0;i<nSOids;i++) Fout<<idval[i]<<endl;
#if defined(GASON) || defined(STARON) || defined(BHON)
        for (Int_t i=0;i<nSOids;i++) Fout<<typeval[i]<<endl;
#endif
    }
    if (nSOids>0) idval.clear();
#if defined(GASON) || defined(STARON) || defined(BHON)
    if (nSOids>0) typeval.clear();
#endif

    if (opt.ibinaryout==OUTASCII || opt.ibinaryout==OUTBINARY) Fout.close();
#ifdef USEHDF
    if (opt.ibinaryout==OUTHDF) Fhdf.close();
#endif
#ifdef USEADIOS
    else if (opt.ibinaryout==OUTADIOS) adios_err=adios_close(adios_file_handle);
#endif

#ifdef USEMPI
    MPIFreeWriteComm();
#endif
}

//@}

///\name Final outputs such as properties and output that can be used to construct merger trees and substructure hierarchy
//@{
///Writes the bulk properties of the substructures
///\todo need to add in 500crit mass and radial output in here and in \ref allvars.h
void WriteProperties(Options &opt, const Int_t ngroups, PropData *pdata){
    fstream Fout;
    string fname;
    ostringstream os;
    char buf[40];
    unsigned long long ngtot=0, noffset=0, ng=ngroups;
#ifdef USEPARALLELHDF
    unsigned long long nwritecommtot=0;
#endif
    //if need to convert from physical back to comoving
    if (opt.icomoveunit) {
        opt.p*=opt.h/opt.a;
        for (Int_t i=1;i<=ngroups;i++) pdata[i].ConverttoComove(opt);
    }
#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif

#ifdef USEHDF
    H5OutputFile Fhdf;
    int itemp=0;
#endif
#if defined(USEHDF)||defined(USEADIOS)
    DataGroupNames datagroupnames;
#endif

    PropDataHeader head(opt);

    os << opt.outname <<".properties";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
    if (NProcs > 1) {
        for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
        for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
    }
    else {
        ngtot = ngroups;
        noffset = 0;
    }
#else
    int ThisTask=0,NProcs=1;
    ngtot=ngroups;
#endif
    fname = os.str();
    cout<<"saving property data to "<<fname<<endl;

    //write header
    if (opt.ibinaryout==OUTBINARY) {
        Fout.open(fname,ios::out|ios::binary);
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&ng,sizeof(long unsigned));
        Fout.write((char*)&ngtot,sizeof(long unsigned));
        int hsize=head.headerdatainfo.size();
        Fout.write((char*)&hsize,sizeof(int));
        ///\todo ADD string containing information of what is in output since this will possibly change with time
        for (Int_t i=0;i<head.headerdatainfo.size();i++) {
            strcpy(buf,head.headerdatainfo[i].c_str());
            Fout.write(buf,sizeof(char)*40);
        }
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
        }
        else{
             Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
        }
#else
        Fhdf.create(string(fname));
#endif
        itemp=0;
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            MPI_Allreduce(&ng, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            //if parallel HDF then only
            if (ThisWriteTask==0) {
                Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ThisWriteComm, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &NWriteComms, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &nwritecommtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ngtot, -1, -1, false);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icosmologicalin);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icomoveunit);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.p);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.a);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.lengthtokpc);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.velocitytokms);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.masstosolarmass);
#if defined(GASON) || defined(STARON) || defined(BHON)
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.metallicitytosolar);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.SFRtosolarmassperyear);
                Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.stellaragetoyrs);
#endif
                WriteVELOCIraptorConfigToHDF(opt,Fhdf);
                WriteSimulationInfoToHDF(opt,Fhdf);
                WriteUnitInfoToHDF(opt,Fhdf);
            }
            Fhdf.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
        }
        else{
            Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ThisTask);
            Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &NProcs);
            Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ng);
            Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ngtot);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icosmologicalin);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icomoveunit);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.p);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.a);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.lengthtokpc);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.velocitytokms);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.masstosolarmass);
#if defined(GASON) || defined(STARON) || defined(BHON)
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.metallicitytosolar);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.SFRtosolarmassperyear);
            Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.stellaragetoyrs);
#endif
            WriteVELOCIraptorConfigToHDF(opt,Fhdf);
            WriteSimulationInfoToHDF(opt,Fhdf);
            WriteUnitInfoToHDF(opt,Fhdf);
        }
#else
        Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ThisTask);
        Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &NProcs);
        Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ng);
        Fhdf.write_dataset(opt, datagroupnames.prop[itemp++], 1, &ngtot);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icosmologicalin);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.icomoveunit);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.p);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.a);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.lengthtokpc);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.velocitytokms);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.masstosolarmass);
#if defined(GASON) || defined(STARON) || defined(BHON)
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.metallicitytosolar);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.SFRtosolarmassperyear);
        Fhdf.write_attribute(string("/"), datagroupnames.prop[itemp++], opt.stellaragetoyrs);
#endif
        WriteVELOCIraptorConfigToHDF(opt,Fhdf);
        WriteSimulationInfoToHDF(opt,Fhdf);
        WriteUnitInfoToHDF(opt,Fhdf);
#endif
    }
#endif
    else {
        Fout.open(fname,ios::out);
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ngroups<<" "<<ngtot<<endl;
        for (Int_t i=0;i<head.headerdatainfo.size();i++) Fout<<head.headerdatainfo[i]<<"("<<i+1<<") ";Fout<<endl;
        Fout<<setprecision(10);
    }

    // long long idbound;
    //for ensuring downgrade of precision as subfind uses floats when storing values save for Mvir (??why??)
    // float value,ctemp[3],mtemp[9];
    // double dvalue;
    // int ivalue;
    for (Int_t i=1;i<=ngroups;i++) {
        if (opt.ibinaryout==OUTBINARY) {
            pdata[i].WriteBinary(Fout,opt);
        }
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            //pdata[i].WriteHDF(Fhdf);
            //for hdf may be more useful to produce an array of the appropriate size and write each data set in one go
            //requires allocating memory
        }
#endif
        else if (opt.ibinaryout==OUTASCII){
            pdata[i].WriteAscii(Fout,opt);
        }
    }
#ifdef USEHDF
    if (opt.ibinaryout==OUTHDF) {
        //for hdf may be more useful to produce an array of the appropriate size and write each data set in one go
        //requires allocating memory
        int itemp;
        //void pointer to hold data
        void *data;
        //allocate enough memory to store largest data type
        data= ::operator new(sizeof(long long)*(ng+1));
        itemp=0;

        //first is halo ids, then id of most bound particle, host halo id, number of direct subhaloes, number of particles
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].haloid;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((long long*)data)[i]=pdata[i+1].ibound;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((long long*)data)[i]=pdata[i+1].iminpot;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((long long*)data)[i]=pdata[i+1].hostid;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].numsubs;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].num;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((int*)data)[i]=pdata[i+1].stype;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        if (opt.iKeepFOF==1){
            for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].directhostid;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].hostfofid;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        //now halo properties that are doubles
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gMvir;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gcm[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gposmbp[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gposminpot[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gcmvel[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gvelmbp[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gvelminpot[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gmass;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gMFOF;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gM200m;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gM200c;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gMBN98;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Efrac;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRvir;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gsize;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gR200m;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gR200c;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRBN98;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRhalfmass;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRmaxvel;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRhalf200m;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRhalf200c;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRhalfBN98;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gmaxvel;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gsigma_v;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gveldisp(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].glambda_B;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gJ[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gq;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gs;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].geigvec(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cNFW;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cNFW200c;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cNFW200m;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cNFWBN98;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Krot;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].T;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Pot;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;


        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].RV_sigma_v;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].RV_veldisp(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].RV_lambda_B;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].RV_J[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].RV_q;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].RV_s;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].RV_eigvec(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        if (opt.iextrahalooutput) {
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gJ200m[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gJ200c[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gJBN98[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            if (opt.iInclusiveHalo>0) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gM200m_excl;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gM200c_excl;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gMBN98_excl;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gR200m_excl;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gR200c_excl;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gRBN98_excl;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

                for (int k=0;k<3;k++){
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gJ200m_excl[k];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
                for (int k=0;k<3;k++){
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gJ200c_excl[k];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
                for (int k=0;k<3;k++){
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].gJBN98_excl[k];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
            }
        }

#ifdef GASON
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].n_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_gas_rvmax;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_gas_30kpc;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_gas_500c;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cm_gas[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cmvel_gas[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Efrac_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Rhalfmass_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].veldisp_gas(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_gas[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].q_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].s_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].eigvec_gas(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Krot_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Temp_mean_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
#ifdef STARON
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Z_mean_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SFR_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
#endif
    if (opt.iextragasoutput) {
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_gas;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_gas[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_gas[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_gas[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        if (opt.iInclusiveHalo>0) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_excl_gas;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_excl_gas;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_excl_gas;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_excl_gas[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_excl_gas[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_excl_gas[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
        }
    }
#endif

#ifdef STARON
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].n_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_star_rvmax;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_star_30kpc;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_star_500c;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cm_star[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].cmvel_star[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Efrac_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Rhalfmass_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].veldisp_star(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_star[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].q_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].s_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].eigvec_star(k,n);
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Krot_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].t_mean_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Z_mean_star;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        if (opt.iextrastaroutput) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_star;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_star;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_star;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_star[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_star[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_star[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }

            if (opt.iInclusiveHalo>0) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_excl_star;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_excl_star;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_excl_star;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

                for (int k=0;k<3;k++){
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_excl_star[k];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
                for (int k=0;k<3;k++){
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_excl_star[k];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
                for (int k=0;k<3;k++){
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_excl_star[k];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
            }
        }
#endif
#ifdef BHON
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].n_bh;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_bh;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
#endif
#ifdef HIGHRES
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].n_interloper;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_interloper;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        if (opt.iextrainterloperoutput) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_interloper;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_interloper;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_interloper;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            if (opt.iInclusiveHalo>0) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_excl_interloper;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_excl_interloper;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_excl_interloper;
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
        }
#endif

#if defined(GASON) && defined(STARON)
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_gas_sf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Rhalfmass_gas_sf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].sigV_gas_sf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_gas_sf[k];
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    }

    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Krot_gas_sf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Temp_mean_gas_sf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Z_mean_gas_sf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    if (opt.iextragasoutput) {
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_gas_sf;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_gas_sf;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_gas_sf;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_gas_sf[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_gas_sf[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_gas_sf[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        if (opt.iInclusiveHalo>0) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_excl_gas_sf;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_excl_gas_sf;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_excl_gas_sf;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_excl_gas_sf[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_excl_gas_sf[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_excl_gas_sf[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
        }
    }
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_gas_nsf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Rhalfmass_gas_nsf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].sigV_gas_nsf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_gas_nsf[k];
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    }

    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Krot_gas_nsf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Temp_mean_gas_nsf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].Z_mean_gas_nsf;
    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
    if (opt.iextragasoutput) {
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_gas_nsf;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_gas_nsf;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_gas_nsf;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_gas_nsf[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_gas_nsf[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
        for (int k=0;k<3;k++){
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_gas_nsf[k];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }

        if (opt.iInclusiveHalo>0) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200mean_excl_gas_nsf;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_200crit_excl_gas_nsf;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].M_BN98_excl_gas_nsf;
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;

            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200mean_excl_gas_nsf[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_200crit_excl_gas_nsf[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (int k=0;k<3;k++){
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].L_BN98_excl_gas_nsf[k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
        }
    }
#endif

        //output extra hydro/star/bh props
#ifdef GASON
        if (opt.gas_internalprop_names.size() + opt.gas_chem_names.size() + opt.gas_chemproduction_names.size()>0) {
            for (auto &extrafield:opt.gas_internalprop_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].hydroprop.GetInternalProperties(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
            for (auto &extrafield:opt.gas_chem_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].hydroprop.GetChemistry(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
            for (auto &extrafield:opt.gas_chemproduction_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].hydroprop.GetChemistryProduction(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
        }
#endif
#ifdef STARON
        if (opt.star_internalprop_names.size() + opt.star_chem_names.size() + opt.star_chemproduction_names.size()>0) {
            for (auto &extrafield:opt.star_internalprop_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].starprop.GetInternalProperties(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
            for (auto &extrafield:opt.star_chem_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].starprop.GetChemistry(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
            for (auto &extrafield:opt.star_chemproduction_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].starprop.GetChemistryProduction(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
        }
#endif
#ifdef BHON
        if (opt.bh_internalprop_names.size() + opt.bh_chem_names.size() + opt.bh_chemproduction_names.size()>0) {
            for (auto &extrafield:opt.bh_internalprop_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].bhprop.GetInternalProperties(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
            for (auto &extrafield:opt.bh_chem_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].bhprop.GetChemistry(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
            for (auto &extrafield:opt.bh_chemproduction_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].bhprop.GetChemistryProduction(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
        }
#endif
#ifdef EXTRADMON
        if (opt.extra_dm_internalprop_names.size()>0) {
            for (auto &extrafield:opt.extra_dm_internalprop_output_names)
            {
                for (Int_t i=0;i<ngroups;i++)
                    ((Double_t*)data)[i]=pdata[i+1].extradmprop.GetExtraProperties(extrafield);
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);
                itemp++;
            }
        }
#endif


        //output apertures
        if (opt.iaperturecalc && opt.aperturenum>0){
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((unsigned int*)data)[i]=pdata[i+1].aperture_npart[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef GASON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((unsigned int*)data)[i]=pdata[i+1].aperture_npart_gas[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((unsigned int*)data)[i]=pdata[i+1].aperture_npart_gas_sf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((unsigned int*)data)[i]=pdata[i+1].aperture_npart_gas_nsf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((unsigned int*)data)[i]=pdata[i+1].aperture_npart_star[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#ifdef BHON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((unsigned int*)data)[i]=pdata[i+1].aperture_npart_bh[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#ifdef HIGHRES
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((unsigned int*)data)[i]=pdata[i+1].aperture_npart_interloper[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif

            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef GASON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_gas[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_gas_sf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_gas_nsf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_star[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#ifdef BHON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_bh[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#ifdef HIGHRES
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_interloper[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef GASON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_gas[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_gas_sf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_gas_nsf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_star[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_veldisp[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef GASON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_veldisp_gas[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_veldisp_gas_sf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_veldisp_gas_nsf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_veldisp_star[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#if defined(GASON) && defined(STARON)
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_SFR_gas[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_gas[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_gas_sf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_gas_nsf[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_star[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#ifdef GASON
            if (opt.gas_extraprop_aperture_calc) {
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.gas_internalprop_output_names_aperture) {
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_gas[j].GetInternalProperties(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_chem_output_names_aperture){
                    for (Int_t i=0;i<ngroups;i++)
                        ((Double_t*)data)[i]=pdata[i+1].aperture_properties_gas[j].GetChemistry(x);
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_chemproduction_output_names_aperture){
                    for (Int_t i=0;i<ngroups;i++)
                        ((Double_t*)data)[i]=pdata[i+1].aperture_properties_gas[j].GetChemistryProduction(x);
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
            }
        }
#endif
#ifdef STARON
            if (opt.star_extraprop_aperture_calc) {
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.star_internalprop_output_names_aperture) {
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_star[j].GetInternalProperties(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                    }
                }
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.star_chem_output_names_aperture){
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_star[j].GetChemistry(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                    }
                }
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.star_chemproduction_output_names_aperture){
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_star[j].GetChemistryProduction(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                    }
                }
            }
#endif
#ifdef BHON
            if (opt.bh_extraprop_aperture_calc) {
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.bh_internalprop_output_names_aperture) {
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_bh[j].GetInternalProperties(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                    }
                }
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.bh_chem_output_names_aperture){
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_bh[j].GetChemistry(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                    }
                }
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.bh_chemproduction_output_names_aperture){
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_bh[j].GetChemistryProduction(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                    }
                }
            }
#endif
#ifdef EXTTRADMON
            if (opt.extra_dm_extraprop_aperture_calc) {
                for (auto j=0;j<opt.aperturenum;j++) {
                    for (auto &x:opt.extra_dm_internalprop_output_names_aperture) {
                        for (Int_t i=0;i<ngroups;i++)
                            ((Double_t*)data)[i]=pdata[i+1].aperture_properties_extra_dm[j].GetExtraProperties(x);
                        Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                    }
                }
            }
#endif

        }
        //output apertures
        if (opt.iaperturecalc && opt.apertureprojnum>0){
            for (auto k=0;k<3;k++) {
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_proj[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef GASON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_proj_gas[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_proj_gas_sf[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_proj_gas_nsf[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_mass_proj_star[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_proj[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef GASON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_proj_gas[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_proj_gas_sf[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_proj_gas_nsf[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_rhalfmass_proj_star[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
#if defined(GASON) && defined(STARON)
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_SFR_proj_gas[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_proj_gas[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_proj_gas_sf[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_proj_gas_nsf[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].aperture_Z_proj_star[j][k];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#endif
            }
        }
        if (opt.SOnum>0) {
            for (auto j=0;j<opt.SOnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_mass[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
            for (auto j=0;j<opt.SOnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_radius[j];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef GASON
            if (opt.iextragasoutput && opt.iextrahalooutput) {
                for (auto j=0;j<opt.SOnum;j++) {
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_mass_gas[j];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
#ifdef STARON
#endif
            }
#endif
#ifdef STARON
            if (opt.iextrastaroutput && opt.iextrahalooutput) {
                for (auto j=0;j<opt.SOnum;j++) {
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_mass_star[j];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
            }
#endif
#ifdef HIGHRES
            if (opt.iextrainterloperoutput && opt.iextrahalooutput) {
                for (auto j=0;j<opt.SOnum;j++) {
                    for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_mass_interloper[j];
                    Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                }
            }
#endif
        }
        if (opt.SOnum>0 && opt.iextrahalooutput) {
        for (auto j=0;j<opt.SOnum;j++) {
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum[j][0];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum[j][1];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum[j][2];
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
        }
#ifdef GASON
        if (opt.iextragasoutput) {
            for (auto j=0;j<opt.SOnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum_gas[j][0];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum_gas[j][1];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum_gas[j][2];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
#ifdef STARON
#endif
        }
#endif
#ifdef STARON
        if (opt.iextrastaroutput) {
            for (auto j=0;j<opt.SOnum;j++) {
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum_star[j][0];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum_star[j][1];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
                for (Int_t i=0;i<ngroups;i++) ((Double_t*)data)[i]=pdata[i+1].SO_angularmomentum_star[j][2];
                Fhdf.write_dataset(opt, head.headerdatainfo[itemp],ng,data,head.hdfpredtypeinfo[itemp]);itemp++;
            }
        }
#endif
        }
        //delete memory associated with void pointer
        ::operator delete(data);
        // delete[] propdataspace;
        // delete[] propdataset;
    }
#endif
    cout<<"Done"<<endl;
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

    //write the units as metadata for each data set
#ifdef USEHDF
    Fhdf.append(string(fname), H5F_ACC_RDWR, 0, false);
#ifdef USEPARALLELHDF
    if (ThisWriteTask==0) {
#endif
        for (auto ientry=0;ientry<head.headerdatainfo.size();ientry++) {
            WriteHeaderUnitEntry(opt, Fhdf, head.headerdatainfo[ientry], head.unitdatainfo[ientry]);
        }
#ifdef USEPARALLELHDF
    }
#endif
    Fhdf.close();
#endif

#ifdef USEMPI
    MPIFreeWriteComm();
#endif
}

void WriteProfiles(Options &opt, const Int_t ngroups, PropData *pdata){
    fstream Fout;
    string fname;
    ostringstream os;
    char buf[40];
    unsigned long long  ngtot=0, noffset=0, ng=ngroups, nhalos=0, nhalostot=0;
#ifdef USEPARALLELHDF
    unsigned long long nwritecommtot=0, nhalowritecommtot=0;
#endif
    vector<unsigned long long> indices(ngroups), haloindices;

    //void pointer to hold data
    void *data;
    int itemp=0, nbinsedges = opt.profilenbins+1;

    //if need to convert from physical back to comoving
    if (opt.icomoveunit) {
        opt.p*=opt.h/opt.a;
        for (Int_t i=1;i<=ngroups;i++) {
            if (pdata[i].gNFOF >= opt.profileminFOFsize && pdata[i].num >= opt.profileminsize) {
                pdata[i].ConvertProfilestoComove(opt);
            }
        }
    }
    if (opt.iInclusiveHalo>0) {
        nhalos = 0;
        haloindices.resize(ngroups);
        for (auto i=1;i<=ngroups;i++){
            if (pdata[i].gNFOF >= opt.profileminFOFsize && pdata[i].num >= opt.profileminsize && pdata[i].hostid == -1) {
                haloindices[nhalos++] = i;
            }
        }
#ifdef USEMPI
        MPI_Allgather(&nhalos, 1, MPI_Int_t, mpi_nhalos, 1, MPI_Int_t, MPI_COMM_WORLD);
#endif
    }
#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif
#ifdef USEHDF
    H5OutputFile Fhdf;
    vector<hsize_t> dims;
#endif
#if defined(USEHDF)||defined(USEADIOS)
    DataGroupNames datagroupnames;
#endif
    ProfileDataHeader head(opt);

    //since profiles can be called for a subset of objects, get the total number to be written
    if (opt.profileminsize > 0 || opt.profileminFOFsize > 0) {
        ng = 0;
        for (auto i=1;i<=ngroups;i++) if (pdata[i].gNFOF >= opt.profileminFOFsize && pdata[i].num >= opt.profileminsize) {
            indices[ng++] = i;
        }
#ifdef USEMPI
        MPI_Allgather(&ng, 1, MPI_Int_t, mpi_ngroups, 1, MPI_Int_t, MPI_COMM_WORLD);
#endif
    }
    else {
        for (auto i=1;i<=ngroups;i++) indices[i-1] = i;
    }

    os << opt.outname << ".profiles";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
    if (NProcs > 1) {
        for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
        for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
        for (int j=0;j<NProcs;j++) nhalostot+=mpi_nhalos[j];
    }
    else {
        ngtot = ngroups;
        noffset = 0;
        nhalostot = nhalos;
    }
#else
    int ThisTask=0,NProcs=1;
    ngtot=ngroups;
    nhalostot=nhalos;
#endif
    fname = os.str();
    cout<<"saving profiles "<<fname<<endl;
    //allocate enough memory to store largest data type
    data= ::operator new(sizeof(long long)*(opt.profilenbins+1));
    ((Double_t*)data)[0]=0.0;for (auto i=0;i<opt.profilenbins;i++) ((Double_t*)data)[i+1]=opt.profile_bin_edges[i];
    //write header
    if (opt.ibinaryout==OUTBINARY) {
        Fout.open(fname,ios::out|ios::binary);
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&ng,sizeof(long unsigned));
        Fout.write((char*)&ngtot,sizeof(long unsigned));
        Fout.write((char*)&nhalos,sizeof(long unsigned));
        Fout.write((char*)&nhalostot,sizeof(long unsigned));
        int hsize=head.headerdatainfo.size();
        Fout.write((char*)&hsize,sizeof(int));
        strcpy(buf,"Radial_norm");
        Fout.write(buf,sizeof(char)*40);
        strcpy(buf,opt.profileradnormstring.c_str());
        Fout.write(buf,sizeof(char)*40);
        strcpy(buf,"Num_radial_bin_edges");
        Fout.write((char*)&nbinsedges,sizeof(int));
        strcpy(buf,"Radial_bin_edges");
        Fout.write(buf,sizeof(char)*40);
        Fout.write((char*)data,sizeof(Double_t)*(opt.profilenbins+1));

        strcpy(buf,"ID");
        Fout.write(buf,sizeof(char)*40);

    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
            MPI_Allreduce(&ng, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            MPI_Allreduce(&nhalos, &nhalowritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            if (ThisTask == 0) {
                itemp=0;
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ThisWriteComm, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &NWriteComms, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nwritecommtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ngtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nhalowritecommtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nhalostot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, opt.profileradnormstring, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &opt.iInclusiveHalo, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nbinsedges, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], nbinsedges, data, H5T_NATIVE_DOUBLE, -1, false);
            }
            else {
                itemp=10;
            }
            Fhdf.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
        }
        else{
            Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
            itemp=0;
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ThisTask);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &NProcs);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ng);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ngtot);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nhalos);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nhalostot);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, opt.profileradnormstring);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &opt.iInclusiveHalo);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nbinsedges);
            Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], nbinsedges, data, H5T_NATIVE_DOUBLE);
        }

#else
        Fhdf.create(string(fname));
        itemp=0;
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ThisTask);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &NProcs);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ng);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &ngtot);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nhalos);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nhalostot);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, opt.profileradnormstring);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &opt.iInclusiveHalo);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], 1, &nbinsedges);
        Fhdf.write_dataset(opt, datagroupnames.profile[itemp++], nbinsedges, data, H5T_NATIVE_DOUBLE);
#endif

    }
#endif
    else {
        Fout.open(fname,ios::out);
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ngroups<<" "<<ngtot<<endl;
        Fout<<nhalos<<" "<<nhalostot<<endl;
        Fout<<"Radial_norm="<<opt.profileradnormstring.c_str()<<endl;
        Fout<<"Inclusive_profiles_flag="<<opt.iInclusiveHalo<<endl;
        Fout<<nbinsedges<<endl;
        Fout<<"Radial_bin_edges=0 ";for (auto i=0;i<opt.profilenbins;i++) Fout<<opt.profile_bin_edges[i]<<" ";Fout<<endl;
        Fout<<"ID "<<opt.profileradnormstring<<" ";
        //for (auto i=0;i<)
        Fout<<setprecision(10);
        Fout<<endl;
    }
    ::operator delete(data);
    if (opt.ibinaryout==OUTBINARY) {
        //for (Int_t i=1;i<=ngroups;i++) pdata[i].WriteProfileBinary(Fout,opt);
    }
    else if (opt.ibinaryout==OUTASCII){
        //for (Int_t i=1;i<=ngroups;i++) pdata[i].WriteProfileAscii(Fout,opt);
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        itemp=0;
        data= ::operator new(sizeof(long long)*(ng));
        //first is halo ids, then normalisation
        for (auto i=0;i<ng;i++) ((unsigned long*)data)[i]=pdata[indices[i]].haloid;
        Fhdf.write_dataset(opt, head.headerdatainfo[itemp], ng, data, head.hdfpredtypeinfo[itemp]);itemp++;
        if (opt.iprofilenorm == PROFILERNORMR200CRIT) {
            for (Int_t i=0;i<ng;i++) {
                if (opt.iInclusiveHalo >0){
                    if (pdata[indices[i]].hostid == -1) ((Double_t*)data)[i]=pdata[indices[i]].gR200c_excl;
                    else ((Double_t*)data)[i]=pdata[indices[i]].gR200c;
                }
                else ((Double_t*)data)[i]=pdata[indices[i]].gR200c;
            }
            Fhdf.write_dataset(opt, head.headerdatainfo[itemp], ng, data, head.hdfpredtypeinfo[itemp]);itemp++;
        }
        //otherwise no normalisation and don't need to write data block
        ::operator delete(data);

        //now move onto 2d arrays;
        data= ::operator new(sizeof(int)*(ng)*(opt.profilenbins));
        dims.resize(2);dims[0]=ng;dims[1]=opt.profilenbins;
        //write all the npart arrays
        for (Int_t i=0;i<ng;i++) {
            for (auto j=0;j<opt.profilenbins;j++) {
                ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_npart[j];
            }
        }
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;

#ifdef GASON
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_npart_gas[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
#ifdef STARON
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_npart_gas_sf[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_npart_gas_nsf[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
#endif
#endif
#ifdef STARON
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_npart_star[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
#endif

        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_mass[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
#ifdef GASON
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_mass_gas[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
#ifdef STARON
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_mass_gas_sf[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_mass_gas_nsf[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
#endif
#endif
#ifdef STARON
        for (Int_t i=0;i<ng;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[indices[i]].profile_mass_star[j];
        Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
#endif
        ::operator delete(data);

        //write all the npart arrays for halos only if inclusive masses calculated
        if (opt.iInclusiveHalo >0) {
            data= ::operator new(sizeof(int)*(nhalos)*(opt.profilenbins));
            dims.resize(2);dims[0]=nhalos;dims[1]=opt.profilenbins;

            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_npart_inclusive[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #ifdef GASON
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_npart_inclusive_gas[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #ifdef STARON
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_npart_inclusive_gas_sf[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_npart_inclusive_gas_nsf[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #endif
    #endif
    #ifdef STARON
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((unsigned int*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_npart_inclusive_star[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #endif

            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_mass_inclusive[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #ifdef GASON
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_mass_inclusive_gas[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #ifdef STARON
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_mass_inclusive_gas_sf[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_mass_inclusive_gas_nsf[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #endif
    #endif
    #ifdef STARON
            for (Int_t i=0;i<nhalos;i++) for (auto j=0;j<opt.profilenbins;j++) ((float*)data)[i*opt.profilenbins+j]=pdata[haloindices[i]].profile_mass_inclusive_star[j];
            Fhdf.write_dataset_nd(opt, head.headerdatainfo[itemp], 2, dims.data(), data, head.hdfpredtypeinfo[itemp]);itemp++;
    #endif
            ::operator delete(data);
        }
        //delete memory associated with void pointer
        // delete[] profiledataspace;
        // delete[] profiledataset;
    }
#endif
    cout<<"Done"<<endl;
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

#ifdef USEMPI
    MPIFreeWriteComm();
#endif
}

//@}

///\name Writes the hierarchy of structures
void WriteHierarchy(Options &opt, const Int_t &ngroups, const Int_t & nhierarchy, const Int_t &nfield, Int_t *nsub, Int_t *parentgid, Int_t *stype, int subflag){
    fstream Fout;
    fstream Fout2;
    string fname;
    ostringstream os;
    unsigned long long ng=ngroups,ngtot=0,noffset=0;
#ifdef USEPARALLELHDF
     unsigned long long nwritecommtot = 0;
#endif
#ifdef USEMPI
    MPIBuildWriteComm(opt);
#endif
#ifdef USEHDF
    H5OutputFile Fhdf;
    int itemp=0;
#endif
#if defined(USEHDF)||defined(USEADIOS)
    DataGroupNames datagroupnames;
#endif
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

    os << opt.outname << ".catalog_groups";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
    #endif
    fname = os.str();
    cout<<"saving hierarchy data to "<<fname<<endl;

    if (opt.ibinaryout==OUTBINARY) Fout.open(fname,ios::out|ios::binary|ios::app);
#ifdef USEHDF
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            Fhdf.append(string(fname),H5F_ACC_RDWR);
        }
        else{
            Fhdf.append(string(fname),H5F_ACC_RDWR,ThisWriteComm,false);
        }
#else
        Fhdf.append(string(fname),H5F_ACC_RDWR);
#endif
    }
#endif
    else Fout.open(fname,ios::out|ios::app);

    //since the hierarchy file is appended to the catalog_groups files, no header written
#ifdef USEMPI
    if (NProcs > 1) {
        for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
        for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
    }
    else {
        ngtot = ngroups;
        noffset = 0;
    }
#else
    ngtot=ngroups;
#endif

    //if subflag==0 only write number of substructures
    if (subflag==0) {
        if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&nsub[1],sizeof(Int_t)*nfield);
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            itemp=4;
            unsigned int *data=new unsigned int[nfield];
            for (Int_t i=1;i<=nfield;i++) data[i-1]=nsub[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp], nfield, data);
            delete[] data;
        }
#endif
        else for (Int_t i=1;i<=nfield;i++)Fout<<nsub[i]<<endl;
    }
    else if (subflag==1) {
        if (opt.ibinaryout==OUTBINARY) {
            Fout.write((char*)&nsub[1+nfield],sizeof(Int_t)*(ngroups-nfield));
            Fout.write((char*)&parentgid[nfield+1],sizeof(Int_t)*(ngroups-nfield));
        }
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            itemp=4;
            unsigned int *data=new unsigned int[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data[i-nfield-1]=nsub[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups-nfield, data);
            delete[] data;
            long long *data2=new long long[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data2[i-nfield-1]=parentgid[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups-nfield, data2);
            delete[] data2;
        }
#endif
        else {
            for (Int_t i=nfield+1;i<=ngroups;i++)Fout<<nsub[i]<<endl;
            for (Int_t i=nfield+1;i<=ngroups;i++)Fout<<parentgid[i]<<endl;
        }
    }
    //write everything, no distinction made between field and substructure
    else if (subflag==-1) {
        if (opt.ibinaryout==OUTBINARY) {
            Fout.write((char*)&nsub[1],sizeof(Int_t)*ngroups);
            Fout.write((char*)&parentgid[1],sizeof(Int_t)*ngroups);
        }
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            itemp=4;
            unsigned int *data=new unsigned int[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data[i-1]=nsub[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups, data);
            delete[] data;
            long long *data2=new long long[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data2[i-1]=parentgid[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups, data2);
            delete[] data2;
        }
#endif
        else {
            for (Int_t i=1;i<=ngroups;i++)Fout<<nsub[i]<<endl;
            for (Int_t i=1;i<=ngroups;i++)Fout<<parentgid[i]<<endl;
        }
    }
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

    //now write a completely separate hierarchy file which I find more intuitive to parse
    os.str(string());
    os << opt.outname << ".hierarchy";
#ifdef USEMPI
    if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        os<<"."<<ThisWriteComm;
#else
        os<<"."<<ThisTask;
#endif
    }
    else {
        os<<"."<<ThisTask;
    }
#endif
    fname = os.str();
    if (opt.ibinaryout==OUTBINARY) Fout.open(fname,ios::out|ios::binary);
    else Fout.open(fname,ios::out);

    cout<<"saving hierarchy data to "<<fname<<endl;
    if (opt.ibinaryout==OUTBINARY) {
        Fout.write((char*)&ThisTask,sizeof(int));
        Fout.write((char*)&NProcs,sizeof(int));
        Fout.write((char*)&ng,sizeof(unsigned long));
        Fout.write((char*)&ngtot,sizeof(unsigned long));
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
#ifdef USEPARALLELHDF
        if(opt.mpinprocswritesize>1){
            //if parallel then open file in serial so task 0 writes header
            Fhdf.create(string(fname),H5F_ACC_TRUNC, 0, false);
            MPI_Allreduce(&ng, &nwritecommtot, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, mpi_comm_write);
            if (ThisWriteTask == 0) {
                itemp=0;
                Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ThisWriteComm, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &NWriteComms, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &nwritecommtot, -1, -1, false);
                Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ngtot, -1, -1, false);
            }
            else {
                itemp = 4;
            }
            Fhdf.close();
            MPI_Barrier(MPI_COMM_WORLD);
            //reopen for parallel write
            Fhdf.append(string(fname));
        }
        else{
            Fhdf.create(string(fname),H5F_ACC_TRUNC, ThisWriteComm, false);
            itemp=0;
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ThisTask);
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &NProcs);
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ng);
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ngtot);
        }
#else
        Fhdf.create(string(fname));
        itemp=0;
        Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ThisTask);
        Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &NProcs);
        Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ng);
        Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], 1, &ngtot);
#endif
    }
#endif
    else {
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ng<<" "<<ngtot<<endl;
        Fout<<setprecision(10);
    }

    if (subflag==0) {
        if (opt.ibinaryout==OUTBINARY) {
            Fout.write((char*)&parentgid[1],sizeof(Int_t)*nfield);
            Fout.write((char*)&nsub[1],sizeof(Int_t)*nfield);
        }
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            itemp=4;
            unsigned int *data=new unsigned int[nfield];
            for (Int_t i=1;i<=nfield;i++) data[i-1]=nsub[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], nfield, data);
            delete[] data;
            long long *data2=new long long[nfield];
            for (Int_t i=1;i<=nfield;i++) data2[i-1]=parentgid[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], nfield, data2);
            delete[] data2;
        }
#endif
        else for (Int_t i=1;i<=nfield;i++)Fout<<parentgid[i]<<" "<<nsub[i]<<endl;
    }
    else if (subflag==1) {
        if (opt.ibinaryout==OUTBINARY) {
            Fout.write((char*)&parentgid[1+nfield],sizeof(Int_t)*(ngroups-nfield));
            Fout.write((char*)&nsub[1+nfield],sizeof(Int_t)*(ngroups-nfield));
        }
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            unsigned int *data=new unsigned int[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data[i-nfield-1]=nsub[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups-nfield, data);
            delete[] data;
            long long *data2=new long long[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data2[i-nfield-1]=parentgid[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups-nfield, data2);
            delete[] data2;
        }
#endif
        else for (Int_t i=nfield+1;i<=ngroups;i++)Fout<<parentgid[i]<<" "<<nsub[i]<<endl;
    }
    //write everything, no distinction made between field and substructure
    else if (subflag==-1) {
        if (opt.ibinaryout==OUTBINARY) {
            Fout.write((char*)&nsub[1],sizeof(Int_t)*ngroups);
            Fout.write((char*)&parentgid[1],sizeof(Int_t)*ngroups);
        }
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            unsigned int *data=new unsigned int[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data[i-1]=nsub[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups, data);
            delete[] data;
            long long *data2=new long long[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data2[i-1]=parentgid[i];
            Fhdf.write_dataset(opt, datagroupnames.hierarchy[itemp++], ngroups, data2);
            delete[] data2;
        }
#endif
        else {
            for (Int_t i=1;i<=ngroups;i++)Fout<<nsub[i]<<endl;
            for (Int_t i=1;i<=ngroups;i++)Fout<<parentgid[i]<<endl;
        }
    }
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif
#ifdef USEMPI
    MPIFreeWriteComm();
#endif
    cout<<"Done saving hierarchy"<<endl;
}

///Write subfind style format of properties, where selection of properties
///are written for FOF objects and a larger selection of properties are written
///for each object
void WriteSUBFINDProperties(Options &opt, const Int_t ngroups, PropData *pdata){
    fstream Fout;
    char fname[1000];
    // char buf[40];
    unsigned long long ngtot=0, noffset=0, ng=ngroups;

    //if need to convert from physical back to comoving
    if (opt.icomoveunit) {
        opt.p*=opt.h/opt.a;
        for (Int_t i=1;i<=ngroups;i++) pdata[i].ConverttoComove(opt);
    }
#ifdef USEHDF
    PropDataHeader head(opt);

#ifdef USEMPI
    sprintf(fname,"%s.subfindproperties.%d",opt.outname,ThisTask);
    for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
    for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
#else
    sprintf(fname,"%s.subproperties",opt.outname);
    int ThisTask=0,NProcs=1;
    ngtot=ngroups;
#endif
    cout<<"saving property data to "<<fname<<endl;
#endif
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
    if (opt.ibinaryout==OUTBINARY) {
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
    // Int_t temp;
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

///\name Write configuration/simulation info
//@{
void WriteVELOCIraptorConfig(Options &opt){
    fstream Fout;
    char fname[1000];
#ifndef USEMPI
    int ThisTask=0;
#endif

    if (ThisTask==0) {
        ConfigInfo config(opt);
        sprintf(fname,"%s.configuration",opt.outname);
        Fout.open(fname,ios::out);
#ifndef OLDCCOMPILER
        for (Int_t i=0;i<config.nameinfo.size();i++) {
            Fout<<config.nameinfo[i]<<"=";
            Fout<<config.datainfo[i]<<" # ";
            Fout<<config.datatype[i]<<" ";
            Fout<<endl;
        }
#else
        Fout<<"C compiler is too old and config file output relies on std 11 implentation to write info. UPDATE YOUR COMPILER "<<endl;
#endif
        Fout.close();
    }
}


#ifdef USEHDF
void WriteVELOCIraptorConfigToHDF(Options &opt, H5OutputFile &Fhdf){
    hid_t group_id = Fhdf.create_group("Configuration");
    ConfigInfo config(opt);
    for (auto i=0;i<config.nameinfo.size();i++) {
        Fhdf.write_attribute(string("/Configuration"), config.nameinfo[i], config.datainfo[i]);
    }
    Fhdf.close_group(group_id);
}
#endif


void WriteSimulationInfo(Options &opt){
    fstream Fout;
    char fname[1000];
#ifndef USEMPI
    int ThisTask=0;
#endif

    if (ThisTask==0) {
        SimInfo siminfo(opt);
        sprintf(fname,"%s.siminfo",opt.outname);
        Fout.open(fname,ios::out);
#ifndef OLDCCOMPILER
        for (Int_t i=0;i<siminfo.nameinfo.size();i++) {
            Fout<<siminfo.nameinfo[i]<<" : ";
            Fout<<siminfo.datainfo[i]<<" # ";
            Fout<<siminfo.datatype[i]<<" ";
            Fout<<endl;
        }
#else
        Fout<<"C compiler is too old and config file output relies on std 11 implentation to write info. UPDATE YOUR COMPILER "<<endl;
#endif
        Fout.close();
    }
}

#ifdef USEHDF
void WriteSimulationInfoToHDF(Options &opt, H5OutputFile &Fhdf){
    hid_t group_id = Fhdf.create_group("SimulationInfo");
    SimInfo siminfo(opt);
    for (auto i=0;i<siminfo.nameinfo.size();i++) {
        Fhdf.write_attribute(string("/SimulationInfo"), siminfo.nameinfo[i], siminfo.datainfo[i]);
    }
    Fhdf.close_group(group_id);
}
#endif

void WriteUnitInfo(Options &opt){
    fstream Fout;
    char fname[1000];
#ifndef USEMPI
    int ThisTask=0;
#endif

    if (ThisTask==0) {
        UnitInfo unitinfo(opt);
        sprintf(fname,"%s.units",opt.outname);
        Fout.open(fname,ios::out);
#ifndef OLDCCOMPILER
        for (Int_t i=0;i<unitinfo.nameinfo.size();i++) {
            Fout<<unitinfo.nameinfo[i]<<" : ";
            Fout<<unitinfo.datainfo[i]<<" # ";
            Fout<<unitinfo.datatype[i]<<" ";
            Fout<<endl;
        }
#else
        Fout<<"C compiler is too old and config file output relies on std 11 implentation to write info. UPDATE YOUR COMPILER "<<endl;
#endif
        Fout.close();
    }
}

#ifdef USEHDF
void WriteUnitInfoToHDF(Options &opt, H5OutputFile &Fhdf){
    hid_t group_id = Fhdf.create_group("UnitInfo");
    UnitInfo unitinfo(opt);
    for (auto i=0;i<unitinfo.nameinfo.size();i++) {
        Fhdf.write_attribute(string("/UnitInfo"), unitinfo.nameinfo[i], unitinfo.datainfo[i]);
    }
    Fhdf.close_group(group_id);
}
#endif

//@}

///\name output simulation state
//@{
void PrintCosmology(Options &opt){
    if (opt.iverbose) {
        cout<<"Cosmology (h, Omega_m, Omega_cdm, Omega_b, Omega_L, Omega_r, Omega_nu, Omega_k, Omega_de, w_de) =";
        cout<<"("<<opt.h<<", ";
        cout<<opt.Omega_m<<", ";
        cout<<opt.Omega_cdm<<", ";
        cout<<opt.Omega_b<<", ";
        cout<<opt.Omega_Lambda<<", ";
        cout<<opt.Omega_r<<", ";
        cout<<opt.Omega_nu<<", ";
        cout<<opt.Omega_k<<", ";
        cout<<opt.Omega_de<<", ";
        cout<<opt.w_de<<", ";
        cout<<")"<<endl;
    }
}

void PrintSimulationState(Options &opt){
    if (opt.iverbose) {
        cout<<"Current simulation state "<<endl;
        cout<<"Scale factor :"<<opt.a<<endl;
        cout<<"Period :"<<opt.p<<endl;
        if (opt.icosmologicalin) {
            double Hubble=opt.h*opt.H*sqrt(opt.Omega_k*pow(opt.a,-2.0)+opt.Omega_m*pow(opt.a,-3.0)
            +opt.Omega_r*pow(opt.a,-4.0)+opt.Omega_Lambda+opt.Omega_de*pow(opt.a,-3.0*(1+opt.w_de)));
            cout<<"Cosmological simulation with "<<endl;
            cout<<"Hubble expansion :"<<Hubble<<endl;
            cout<<"Critical Density :"<<opt.rhobg/opt.Omega_m<<endl;
            cout<<"Matter density :"<<opt.rhobg<<endl;
        }
    }
}
//@}

#ifdef SWIFTINTERFACE
///write an HDF file that stores where particles are written.
void WriteSwiftExtendedOutput(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, vector<Particle> &Part)
{
    return;
}
#endif

#ifdef EXTENDEDHALOOUTPUT
/// \name Routines that can be used to output information of a halo subvolume decomposition
//@{
///Writes cell quantites
void WriteExtendedOutput (Options &opt, Int_t numgroups, Int_t nbodies, PropData *pdata, Particle *p, Int_t * pfof)
{
    fstream Fout;
    char fname[1000];
    int ngtot = 0;
    int noffset = 0;

#ifdef USEMPI
    for (int j = 0; j < NProcs; j++) ngtot += mpi_ngroups[j];
    for (int j = 0; j < ThisTask; j++) noffset += mpi_ngroups[j];
#else
    int ThisTask = 0;
    int NProcs = 1;
    ngtot = numgroups;
#endif
    cout << "numgroups  " << numgroups << endl;
    numgroups++;

    Int_t *  nfilespergroup = new Int_t   [numgroups];
    Int_t ** filesofgroup   = new Int_t * [numgroups];

    Int_t *  ntaskspergroup = new Int_t   [numgroups];
    Int_t ** tasksofgroup   = new Int_t * [numgroups];

    // First send to Task 0, Groups and number of files over which the group is distributed
    // from original reading

    // Groups id have already the offset
    Int_t ** npartsofgroupinfile = new Int_t * [numgroups];
    Int_t ** npartsofgroupintask = new Int_t * [numgroups];

#ifdef USEMPI
    Int_t ntosendtotask      [NProcs];
    Int_t ntorecievefromtask [NProcs];
    for (Int_t i = 0; i < NProcs; i++)
    {
        ntosendtotask[i] = 0;
        ntorecievefromtask[i] = 0;
    }
#endif

    // Initialize arrays
    for (Int_t i = 0; i < numgroups; i++)
    {
        npartsofgroupinfile[i] = new Int_t [opt.num_files];

#ifdef USEMPI
        npartsofgroupintask[i] = new Int_t [NProcs];
        ntaskspergroup[i] = 0;
        for(Int_t j = 0; j < NProcs; j++)
            npartsofgroupintask[i][j] = 0;
#endif

        nfilespergroup[i] = 0;
        for(Int_t j = 0; j < opt.num_files; j++)
            npartsofgroupinfile[i][j] = 0;
    }

    // Set IdStruct IdTopHost and fill arrays
    for (Int_t i = 0; i < nbodies; i++)
    {
        if (pfof[i] == 0)
        {
            p[i].SetIdTopHost(0);
            p[i].SetIdHost(0);
            p[i].SetIdStruct(0);
        }
        else
        {
            p[i].SetIdStruct (pdata[pfof[i]].haloid);
            if (pdata[pfof[i]].hostfofid == 0)
                p[i].SetIdTopHost (pfof[i] + noffset);
            else
                p[i].SetIdTopHost (pdata[pfof[i]].hostfofid);
            if (pdata[pfof[i]].hostid < 0)
                p[i].SetIdHost (pfof[i] + noffset);
            else
                p[i].SetIdHost (pdata[pfof[i]].hostid);
        }

        npartsofgroupinfile[pfof[i]][p[i].GetOFile()]++;
#ifdef USEMPI
        npartsofgroupintask[pfof[i]][p[i].GetOTask()]++;
#endif
    }

    for (Int_t i = 0; i < numgroups; i++)
    {
#ifdef USEMPI
        for(Int_t j = 0; j < NProcs; j++)
            if(npartsofgroupintask[i][j] > 0)
                ntaskspergroup[i]++;
#endif
        for(Int_t j = 0; j < opt.num_files; j++)
            if(npartsofgroupinfile[i][j] > 0)
                nfilespergroup[i]++;
    }

    // Shorten arrays
    for (Int_t i = 0; i < numgroups; i++)
    {
        Int_t k = 0;
#ifdef USEMPI
        tasksofgroup[i] = new Int_t [ntaskspergroup[i]];
        for(Int_t j = 0; j < NProcs; j++)
        {
            if(npartsofgroupintask[i][j] > 0)
            {
                tasksofgroup[i][k] = j;
                ntosendtotask[j] += npartsofgroupintask[i][j];
                k++;
            }
        }
        k = 0;
#endif
        filesofgroup[i] = new Int_t [nfilespergroup[i]];
        for(Int_t j = 0; j < opt.num_files; j++)
        {
            if(npartsofgroupinfile[i][j] > 0)
            {
                filesofgroup[i][k] = j;
                k++;
            }
        }
        delete [] npartsofgroupinfile[i];
#ifdef USEMPI
        delete [] npartsofgroupintask[i];
#endif
    }
    delete [] npartsofgroupinfile;
#ifdef USEMPI
    delete [] npartsofgroupintask;
#endif
    // Now nfilespergroup has the number of files over which the group is distributed and
    // filesofgroup has the id of the files
    //
    // Write FilesOfGroup File
    int myturn = 0;

#ifdef USEMPI
    MPI_Status status;
    MPI_Request rqst;

    if (ThisTask == 0)
    {
#endif
        myturn = 1;
        char fog [1000];
        sprintf (fog, "%s.filesofgroup", opt.outname);
        Fout.open (fog, ios::out);
        for (Int_t i = 1; i < numgroups; i++)
        {
            Fout << pdata[i].haloid << "  " << nfilespergroup[i] << endl;
            for (Int_t j = 0; j < nfilespergroup[i]; j++)
            Fout << filesofgroup[i][j] << " ";
            Fout << endl;
        }
        Fout.close();

#ifdef USEMPI
        if (NProcs > 1)
            MPI_Isend (&myturn, 1, MPI_INT, ThisTask+1, ThisTask, MPI_COMM_WORLD, &rqst);
    }
    else
    {
        MPI_Recv (&myturn, 1, MPI_INT, ThisTask-1, ThisTask-1, MPI_COMM_WORLD, &status);

        ofstream fout;
        char fog[1000];
        sprintf (fog, "%s.filesofgroup", opt.outname);
        fout.open (fog, ios::app);
        for (Int_t i = 1; i < numgroups; i++)
        {
            fout << pdata[i].haloid << " " << nfilespergroup[i] << endl;
            for (Int_t j = 0; j < nfilespergroup[i]; j++)fout << filesofgroup[i][j] << " ";
            fout << endl;
        }
        fout.close();

        if (ThisTask != NProcs-1)
            MPI_Isend (&myturn, 1, MPI_INT, ThisTask+1, ThisTask, MPI_COMM_WORLD, &rqst);
    }
    delete [] filesofgroup;
    delete [] tasksofgroup;
    delete [] ntaskspergroup;
    delete [] nfilespergroup;

    if (opt.iverbose) cout << ThisTask << "filesofgroup written" << endl;
    // Send and Particles before writing Extended Files
    // Communicate to all other processors how many particles are going to be sent
    ntosendtotask[ThisTask] = 0;
    for (Int_t i = 1; i < NProcs; i++)
    {
        int src = (ThisTask + NProcs - i) % NProcs;
        int dst = (ThisTask + i) % NProcs;
        MPI_Isend (&ntosendtotask[dst], 1, MPI_INT, dst, ThisTask, MPI_COMM_WORLD, &rqst);
        MPI_Recv (&ntorecievefromtask[src], 1, MPI_INT, src, src, MPI_COMM_WORLD, &status);
    }

    // Declare and allocate Particle arrays for sending and receiving
    Particle ** PartsToSend = new Particle * [NProcs];
    Particle ** PartsToRecv = new Particle * [NProcs];

    int * count = new int [NProcs];
    for (Int_t i = 0; i < NProcs; i++)
    {
        count[i] = 0;
        PartsToSend[i] = new Particle [ntosendtotask[i]+1];
        PartsToRecv[i] = new Particle [ntorecievefromtask[i]+1];
    }

    // Copy Particles to send
    for (Int_t i = 0; i < nbodies; i++)
        if (p[i].GetOTask() != ThisTask)
            PartsToSend[p[i].GetOTask()][count[p[i].GetOTask()]++] = Particle(p[i]);

    //determine if number of particles can fit into a single send
    int bufferFlag = 1;
    long int  maxNumPart = LOCAL_MAX_MSGSIZE / (long int) sizeof(Particle);
    int localMax = 0;
    int globalMax = 0;
    //find max local send and global send
    for (Int_t i = 0; i < NProcs; i++) if (ntosendtotask[i] > localMax) localMax = ntosendtotask[i];
    MPI_Allreduce (&localMax, &globalMax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if (globalMax >= maxNumPart) bufferFlag = 1;

    // Send and Receive Particles
    //if splitting sends into chunks
    if (bufferFlag)
    {
        int numBuffersToSend [NProcs];
        int numBuffersToRecv [NProcs];
        int numPartInBuffer = maxNumPart;

        for (int jj = 0; jj < NProcs; jj++)
        {
            numBuffersToSend[jj] = 0;
            numBuffersToRecv[jj] = 0;
            if (ntosendtotask[jj] > 0) numBuffersToSend[jj] = (ntosendtotask[jj]/numPartInBuffer) + 1;
        }
        //broadcast numbers sent
        for (int i = 1; i < NProcs; i++)
        {
            int src = (ThisTask + NProcs - i) % NProcs;
            int dst = (ThisTask + i) % NProcs;
            MPI_Isend (&numBuffersToSend[dst], 1, MPI_INT, dst, 0, MPI_COMM_WORLD, &rqst);
            MPI_Recv  (&numBuffersToRecv[src], 1, MPI_INT, src, 0, MPI_COMM_WORLD, &status);
        }
        MPI_Barrier (MPI_COMM_WORLD);

        //for each mpi thread send info as necessary
        for (int i = 1; i < NProcs; i++)
        {
            int src = (ThisTask + NProcs - i) % NProcs;
            int dst = (ThisTask + i) % NProcs;
            Int_t size = numPartInBuffer;
            int buffOffset = 0;
            //send buffers
            for (int jj = 0; jj < numBuffersToSend[dst]-1; jj++)
            {
              MPI_Isend (&size, 1, MPI_Int_t, dst, (int)(jj+1), MPI_COMM_WORLD, &rqst);
              MPI_Isend (&PartsToSend[dst][buffOffset], sizeof(Particle)*size, MPI_BYTE,
                         dst, (int)(10000+jj+1), MPI_COMM_WORLD, &rqst);
            }
            //and if anything is remaining
            size = ntosendtotask[dst] % numPartInBuffer;
            if (size > 0 && numBuffersToSend[dst] > 0)
            {
              MPI_Isend (&size, 1, MPI_Int_t, dst, (int)(numBuffersToSend[dst]), MPI_COMM_WORLD, &rqst);
              MPI_Isend (&PartsToSend[dst][buffOffset], sizeof(Particle)*size, MPI_BYTE,
                         dst, (int)(10000+numBuffersToSend[dst]), MPI_COMM_WORLD, &rqst);
            }

            // Receive Buffers
            buffOffset = 0;
            for (int jj = 0; jj < numBuffersToRecv[src]; jj++)
            {
              Int_t numInBuffer = 0;
              MPI_Recv (&numInBuffer, 1, MPI_Int_t, src, (int)(jj+1), MPI_COMM_WORLD, &status);
              MPI_Recv (&PartsToRecv[src][buffOffset], sizeof(Particle)*numInBuffer,
                        MPI_BYTE, src, (int)(10000+jj+1), MPI_COMM_WORLD, &status);
              buffOffset += numInBuffer;
            }
        }
    }
    else
    {
        for (Int_t i = 1; i < NProcs; i++)
        {
            int src = (ThisTask + NProcs - i) % NProcs;
            int dst = (ThisTask + i) % NProcs;
            MPI_Isend (PartsToSend[dst], ntosendtotask[dst]*sizeof(Particle), MPI_BYTE, dst, ThisTask, MPI_COMM_WORLD, &rqst);
            MPI_Recv (PartsToRecv[src], ntorecievefromtask[src]*sizeof(Particle), MPI_BYTE, src, src, MPI_COMM_WORLD, &status);
        }
    }
#endif

    // Organize particles in files for the ExtendedOutput
    // First determine which files ThisTask should write
    int npartthistask = 0;
    int npartperfile[opt.num_files];

    for (Int_t i = 0; i < opt.num_files; i++)
        npartperfile[i] = 0;

    // Loop over particles to get how many go to each file
    for (Int_t i = 0; i < nbodies; i++) if (p[i].GetOTask() == ThisTask)
    {
        npartthistask++;
        npartperfile[p[i].GetOFile()]++;
    }

#ifdef USEMPI
    for (Int_t i = 0; i < NProcs; i++)
        for (Int_t j = 0; j < ntorecievefromtask[i]; j++)
            if (PartsToRecv[i][j].GetOTask() == ThisTask)
            {
                npartthistask++;
                npartperfile[PartsToRecv[i][j].GetOFile()]++;
            }
#endif

    // Allocate and fill arrays of particle, halo, host, and igm ids
    // Here I am assuming that we have n particles with indexes rangin from 0 to n-1
    int ** Id, ** IdStruct, ** IdTopHost, ** IdHost;

    Id       = new int * [opt.num_files];
    IdStruct = new int * [opt.num_files];
    IdTopHost    = new int * [opt.num_files];
    IdHost   = new int * [opt.num_files];

    for (Int_t i = 0; i < opt.num_files; i++)
    if (npartperfile[i] > 0)
    {
        Id[i]       = new int [npartperfile[i]];
        IdStruct[i] = new int [npartperfile[i]];
        IdTopHost[i]    = new int [npartperfile[i]];
        IdHost[i]   = new int [npartperfile[i]];
    }

    for (Int_t i = 0; i < nbodies; i++)
    {
        if (p[i].GetOTask() == ThisTask)
        {
            Id       [p[i].GetOFile()][p[i].GetOIndex()] = p[i].GetPID();
            IdStruct [p[i].GetOFile()][p[i].GetOIndex()] = p[i].GetIdStruct();
            IdHost   [p[i].GetOFile()][p[i].GetOIndex()] = p[i].GetIdHost();
            IdTopHost    [p[i].GetOFile()][p[i].GetOIndex()] = p[i].GetIdTopHost();
        }
    }
#ifdef USEMPI
    for (Int_t i = 0; i < NProcs; i++)
        for (Int_t j = 0; j < ntorecievefromtask[i]; j++)
            if (PartsToRecv[i][j].GetOTask() == ThisTask)
            {
                Id       [PartsToRecv[i][j].GetOFile()][PartsToRecv[i][j].GetOIndex()] = PartsToRecv[i][j].GetPID();
                IdStruct [PartsToRecv[i][j].GetOFile()][PartsToRecv[i][j].GetOIndex()] = PartsToRecv[i][j].GetIdStruct();
                IdHost   [PartsToRecv[i][j].GetOFile()][PartsToRecv[i][j].GetOIndex()] = PartsToRecv[i][j].GetIdHost();
                IdTopHost    [PartsToRecv[i][j].GetOFile()][PartsToRecv[i][j].GetOIndex()] = PartsToRecv[i][j].GetIdTopHost();
            }
#endif

    // Write ExtendedFiles
    for (Int_t i = 0; i < opt.num_files; i++)
        if (npartperfile[i] > 0)
        {
            sprintf (fname,"%s.extended.%d",opt.outname,i);
            Fout.open (fname,ios::out);
            for (Int_t j = 0; j < npartperfile[i]; j++)
            {
                Fout << setw(12) << Id[i][j]       << "  ";
                Fout << setw(7)  << IdStruct[i][j] << "  ";
                Fout << setw(7)  << IdHost[i][j]   << "  ";
                Fout << setw(7)  << IdTopHost[i][j]    << "  ";
                Fout << endl;
            }
            Fout.close();
        }
#ifdef USEMPI
    MPI_Barrier (MPI_COMM_WORLD);
#endif
    cout << ThisTask << " Finished writing extended output" << endl;
}
//@}
#endif
