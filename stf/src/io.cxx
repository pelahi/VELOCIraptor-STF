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

void WriteGroupCatalog(Options &opt, const Int_t ngroups, Int_t *numingroup, Int_t **pglist, Particle *Part, Int_t nadditional){
    fstream Fout,Fout2,Fout3;
    char fname[500];
    char fname2[500];
    char fname3[500];
    unsigned long noffset=0,ngtot=0,nids=0,nidstot,nuids=0,nuidstot,ng=0;
    Int_t *offset;
#ifdef USEHDF
    H5File Fhdf,Fhdf3;
    H5std_string datasetname;
    DataSpace dataspace;
    DataSet dataset;
    hsize_t *dims;
    hsize_t rank;
    HDFCatalogNames hdfnames;
    int itemp=0;
#endif

#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif

#ifdef USEMPI
    sprintf(fname,"%s.catalog_groups.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.catalog_groups",opt.outname);
#endif

    cout<<"saving group catalog to "<<fname<<endl;
    if (opt.ibinaryout==OUTBINARY) Fout.open(fname,ios::out|ios::binary);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        //create file 
        Fhdf=H5File(fname,H5F_ACC_TRUNC);
        //Fhdf.H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    }
#endif    
    else Fout.open(fname,ios::out);
    ng=ngroups;

#ifdef USEMPI
    for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
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
        //set file info
        dims=new hsize_t[1];
        dims[0]=1;
        rank=1;
        itemp=0;
        //datasetname=H5std_string("File_id");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&ThisTask,hdfnames.groupdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Num_of_files");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&NProcs,hdfnames.groupdatatype[itemp]);
        itemp++;
        
        //datasetname=H5std_string("Num_of_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&ng,hdfnames.groupdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Total_num_of_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&ngtot,hdfnames.groupdatatype[itemp]);
        itemp++;
        delete[] dims;
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
        dims=new hsize_t[1];
        dims[0]=ng;
        rank=1;
        //datasetname=H5std_string("Group_size");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        unsigned int *data=new unsigned int[ng];
        for (Int_t i=1;i<=ng;i++) data[i-1]=numingroup[i];
        dataset.write(data,hdfnames.groupdatatype[itemp]);
        itemp++;
        delete[] data;
        delete[] dims;
    }
#endif
    else for (Int_t i=1;i<=ngroups;i++) Fout<<numingroup[i]<<endl;


    //Write offsets for bound and unbound particles
    offset=new Int_t[ngroups+1];
    offset[1]=0;
    //note before had offsets at numingroup but to account for unbound particles use value of pglist at numingroup 
    for (Int_t i=2;i<=ngroups;i++) offset[i]=offset[i-1]+pglist[i-1][numingroup[i-1]];

    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&offset[1],sizeof(Int_t)*ngroups);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        dims=new hsize_t[1];
        dims[0]=ng;
        rank=1;
        //datasetname=H5std_string("Offset");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        unsigned long *data=new unsigned long[ng];
        for (Int_t i=1;i<=ng;i++) data[i-1]=offset[i];
        dataset.write(data,hdfnames.groupdatatype[itemp]);
        itemp++;
        delete[] data;
        delete[] dims;
    }
#endif
    else for (Int_t i=1;i<=ngroups;i++) Fout<<offset[i]<<endl;

    //position of unbound particle
    for (Int_t i=2;i<=ngroups;i++) offset[i]=offset[i-1]+numingroup[i-1]-pglist[i-1][numingroup[i-1]];
    if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&offset[1],sizeof(Int_t)*ngroups);
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        dims=new hsize_t[1];
        dims[0]=ng;
        rank=1;
        //datasetname=H5std_string("Offset_unbound");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        unsigned long *data=new unsigned long[ng];
        for (Int_t i=1;i<=ng;i++) data[i-1]=offset[i];
        dataset.write(data,hdfnames.groupdatatype[itemp]);
        itemp++;
        delete[] data;
        delete[] dims;
    }
#endif
    else for (Int_t i=1;i<=ngroups;i++) Fout<<offset[i]<<endl;

    delete[] offset;
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

    //now write pid files
#ifdef USEMPI
    sprintf(fname,"%s.catalog_particles.%d",opt.outname,ThisTask);
    sprintf(fname3,"%s.catalog_particles.unbound.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.catalog_particles",opt.outname);
    sprintf(fname3,"%s.catalog_particles.unbound",opt.outname);
#endif
    cout<<"saving particle catalog to "<<fname<<endl;

    if (opt.ibinaryout==OUTBINARY) {
        Fout.open(fname,ios::out|ios::binary);
        Fout3.open(fname3,ios::out|ios::binary);
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        Fhdf=H5File(fname,H5F_ACC_TRUNC);
        Fhdf3=H5File(fname3,H5F_ACC_TRUNC);
    }
#endif
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
        //set file info
        dims=new hsize_t[1];
        dims[0]=1;
        rank=1;
        itemp=0;
        //datasetname=H5std_string("File_id");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&ThisTask,hdfnames.partdatatype[itemp]);
        dataset = Fhdf3.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&ThisTask,hdfnames.partdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Num_of_files");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&NProcs,hdfnames.partdatatype[itemp]);
        dataset = Fhdf3.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&NProcs,hdfnames.partdatatype[itemp]);
        itemp++;
        
        //datasetname=H5std_string("Num_of_particles_in_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&nids,hdfnames.partdatatype[itemp]);
        dataset = Fhdf3.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&nuids,hdfnames.partdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Total_num_of_particles_in_all_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&nidstot,hdfnames.partdatatype[itemp]);
        dataset = Fhdf3.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
        dataset.write(&nuidstot,hdfnames.partdatatype[itemp]);
        itemp++;
        delete[] dims;
    }
#endif
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
        if (opt.ibinaryout==OUTBINARY) Fout.write((char*)idval,sizeof(Int_t)*nids);
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            dims=new hsize_t[1];
            dims[0]=nids;
            rank=1;
            //datasetname=H5std_string("Particle_IDs");
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
            long long *data=new long long[nids];
            for (Int_t i=0;i<nids;i++) data[i]=idval[i];
            dataset.write(data,hdfnames.partdatatype[itemp]);
            delete[] data;
            delete[] dims;
        }
#endif
        else for (Int_t i=0;i<nids;i++) Fout<<idval[i]<<endl;
        delete[] idval;
    }
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

    if (nuids>0) {
        idval=new Int_t[nuids];
        nuids=0;
        for (Int_t i=1;i<=ngroups;i++) 
            for (Int_t j=pglist[i][numingroup[i]];j<numingroup[i];j++) 
                idval[nuids++]=Part[pglist[i][j]].GetPID();
        if (opt.ibinaryout==OUTBINARY) Fout3.write((char*)idval,sizeof(Int_t)*nuids);
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            dims=new hsize_t[1];
            dims[0]=nuids;
            rank=1;
            //datasetname=H5std_string("Particle_IDs");
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf3.createDataSet(hdfnames.part[itemp], hdfnames.partdatatype[itemp], dataspace);
            long long *data=new long long[nuids];
            for (Int_t i=0;i<nuids;i++) data[i]=idval[i];
            dataset.write(data,hdfnames.partdatatype[itemp]);
            delete[] data;
            delete[] dims;
        }
#endif
        else for (Int_t i=0;i<nuids;i++) Fout3<<idval[i]<<endl;
        delete[] idval;
    }
    if (opt.ibinaryout!=OUTHDF) Fout3.close();
#ifdef USEHDF
    else Fhdf3.close();
#endif

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

#ifdef USEHDF
    H5File Fhdf,Fhdf2;
    H5std_string datasetname;
    DataSpace dataspace;
    DataSet dataset;
    hsize_t *dims;
    hsize_t rank;
    HDFCatalogNames hdfnames;
    int itemp;
#endif

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


    if (opt.ibinaryout==OUTBINARY) {
        Fout.open(fname,ios::out|ios::binary);
        Fout2.open(fname2,ios::out|ios::binary);
    }
#ifdef USEHDF
    else if (opt.ibinaryout==OUTHDF) {
        //create file 
        Fhdf=H5File(fname,H5F_ACC_TRUNC);
        //Fhdf.H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        Fhdf2=H5File(fname2,H5F_ACC_TRUNC);
    }
#endif    
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
        //set file info
        dims=new hsize_t[1];
        dims[0]=1;
        rank=1;
        itemp=0;
        //datasetname=H5std_string("File_id");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&ThisTask,hdfnames.typesdatatype[itemp]);
        dataset = Fhdf2.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&ThisTask,hdfnames.typesdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Num_of_files");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&NProcs,hdfnames.typesdatatype[itemp]);
        dataset = Fhdf2.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&NProcs,hdfnames.typesdatatype[itemp]);
        itemp++;
        
        //datasetname=H5std_string("Num_of_particles_in_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&nids,hdfnames.typesdatatype[itemp]);
        dataset = Fhdf2.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&nuids,hdfnames.typesdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Total_num_of_particles_in_all_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&nidstot,hdfnames.typesdatatype[itemp]);
        dataset = Fhdf2.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
        dataset.write(&nuidstot,hdfnames.typesdatatype[itemp]);
        itemp++;
        delete[] dims;
    }
#endif
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
        if (opt.ibinaryout==OUTBINARY) Fout.write((char*)typeval,sizeof(int)*nids);
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            dims=new hsize_t[1];
            dims[0]=nids;
            rank=1;
            //datasetname=H5std_string("Particle_types");
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
            unsigned short *data=new unsigned short[nids];
            for (Int_t i=0;i<nids;i++) data[i]=typeval[i];
            dataset.write(data,hdfnames.typesdatatype[itemp]);
            delete[] data;
            delete[] dims;
        }
#endif
        else for (Int_t i=0;i<nids;i++) Fout<<typeval[i]<<endl;
        delete[] typeval;
    }
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

    if (nuids>0) {
        typeval=new int[nuids];
        nuids=0;
        for (Int_t i=1;i<=ngroups;i++) 
            for (Int_t j=pglist[i][numingroup[i]];j<numingroup[i];j++) 
                typeval[nuids++]=Part[pglist[i][j]].GetType();
        if (opt.ibinaryout==OUTBINARY) Fout2.write((char*)typeval,sizeof(int)*nuids);
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            dims=new hsize_t[1];
            dims[0]=nuids;
            rank=1;
            datasetname=H5std_string("Particle_types");
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf2.createDataSet(hdfnames.types[itemp], hdfnames.typesdatatype[itemp], dataspace);
            unsigned short *data=new unsigned short[nuids];
            for (Int_t i=0;i<nuids;i++) data[i]=typeval[i];
            dataset.write(data,hdfnames.typesdatatype[itemp]);
            delete[] data;
            delete[] dims;
        }
#endif
        else for (Int_t i=0;i<nuids;i++) Fout2<<typeval[i]<<endl;
        delete[] typeval;
    }
    if (opt.ibinaryout!=OUTHDF) Fout2.close();
#ifdef USEHDF
    else Fhdf2.close();
#endif

#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//@}


///\name Final outputs such as properties and output that can be used to construct merger trees and substructure hierarchy
//@{
///Writes the bulk properties of the substructures
///\todo need to add in 500crit mass and radial output in here and in \ref allvars.h
void WriteProperties(Options &opt, const Int_t ngroups, PropData *pdata){
    fstream Fout;
    char fname[1000];
    char buf[40];
    long unsigned ngtot=0, noffset=0, ng=ngroups;

#ifdef USEHDF
    H5File Fhdf;
    H5std_string datasetname;
    DataSpace dataspace;
    DataSet dataset;
    hsize_t *dims;
    int rank;
    DataSpace *propdataspace;
    DataSet *propdataset;
    HDFCatalogNames hdfnames;
    int itemp=0;
#endif

    PropDataHeader head;
    
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
        Fhdf=H5File(fname,H5F_ACC_TRUNC);
        //set file info
        dims=new hsize_t[1];
        dims[0]=1;
        rank=1;
        //since header info is the same as the group catalog files, write to hdf file using same interface
        itemp=0;
        //datasetname=H5std_string("File_id");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&ThisTask,hdfnames.groupdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Num_of_files");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&NProcs,hdfnames.groupdatatype[itemp]);
        itemp++;
        
        //datasetname=H5std_string("Num_of_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&ng,hdfnames.groupdatatype[itemp]);
        itemp++;

        //datasetname=H5std_string("Total_num_of_groups");
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.group[itemp], hdfnames.groupdatatype[itemp], dataspace);
        dataset.write(&ngtot,hdfnames.groupdatatype[itemp]);

        //load data spaces
        propdataspace=new DataSpace[head.headerdatainfo.size()];
        propdataset=new DataSet[head.headerdatainfo.size()];
        dims[0]=ng;
        rank=1;
        dataspace=DataSpace(rank,dims);
        for (Int_t i=0;i<head.headerdatainfo.size();i++) {
            datasetname=H5std_string(head.headerdatainfo[i]);
            propdataspace[i]=DataSpace(rank,dims);
            propdataset[i] = Fhdf.createDataSet(datasetname, head.predtypeinfo[i], propdataspace[i]);
        }
        delete[] dims;
    }
#endif
    else {
        Fout.open(fname,ios::out);
        Fout<<ThisTask<<" "<<NProcs<<endl;
        Fout<<ngroups<<" "<<ngtot<<endl;
        for (Int_t i=0;i<head.headerdatainfo.size();i++) Fout<<head.headerdatainfo[i]<<"("<<i+1<<") ";Fout<<endl;
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
        if (opt.ibinaryout==OUTBINARY) {
            pdata[i].WriteBinary(Fout);
        }
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            //pdata[i].WriteHDF(Fhdf);
            //for hdf may be more useful to produce an array of the appropriate size and write each data set in one go
            //requires allocating memory
        }
#endif
        else if (opt.ibinaryout==OUTASCII){
            pdata[i].WriteAscii(Fout);
        }
    }
#ifdef USEHDF
    if (opt.ibinaryout==OUTHDF) {
        //for hdf may be more useful to produce an array of the appropriate size and write each data set in one go
        //requires allocating memory
        int *iarray,itemp;
        unsigned int *uiarray;
        long long *larray;
        unsigned long *ularray;
        double *darray;
        //void pointer to hold data
        void *data;
        //allocate enough memory to store largest data type
        data= ::operator new(sizeof(long long)*(ng+1));
        itemp=0;

        //first is halo ids, then id of most bound particle, host halo id, number of direct subhaloes, number of particles
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].haloid;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((long long*)data)[i]=pdata[i+1].ibound;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((long long*)data)[i]=pdata[i+1].hostid;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].numsubs;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].num;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        
        //now halo properties that are doubles
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gMvir;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gcm[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gpos[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gcmvel[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gvel[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gmass;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gMFOF;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gM200m;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gM200c;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gMvir;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Efrac;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gRvir;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gsize;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gR200m;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gR200c;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gRvir;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gRhalfmass;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gRmaxvel;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gmaxvel;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gsigma_v;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gveldisp(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].glambda_B;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gJ[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gq;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].gs;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].geigvec(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].cNFW;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Krot;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].T;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Pot;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;


        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].RV_sigma_v;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].RV_veldisp(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].RV_lambda_B;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].RV_J[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].RV_q;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].RV_s;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].RV_eigvec(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

#ifdef GASON
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].n_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_gas_rvmax;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_gas_30kpc;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_gas_500c;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;

        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].cm_gas[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].cmvel_gas[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Efrac_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Rhalfmass_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].veldisp_gas(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].L_gas[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].q_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].s_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].eigvec_gas(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Krot_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Temp_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
#ifdef STARON
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Z_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].SFR_gas;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
#endif
#endif

#ifdef STARON
        for (Int_t i=0;i<ngroups;i++) ((unsigned long*)data)[i]=pdata[i+1].n_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_star_rvmax;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_star_30kpc;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].M_star_500c;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;

        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].cm_star[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].cmvel_star[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Efrac_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Rhalfmass_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].veldisp_star(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }
        for (int k=0;k<3;k++){
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].L_star[k];
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].q_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].s_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (int k=0;k<3;k++) for (int n=0;n<3;n++) {
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].eigvec_star(k,n);
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        }

        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Krot_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].t_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
        for (Int_t i=0;i<ngroups;i++) ((double*)data)[i]=pdata[i+1].Z_star;
        propdataset[itemp].write(data,head.predtypeinfo[itemp]);
        itemp++;
#endif
        //delete memory associated with void pointer
        ::operator delete(data);
        delete[] propdataspace;
        delete[] propdataset;
    }
#endif
    cout<<"Done"<<endl;
    if (opt.ibinaryout!=OUTHDF) Fout.close();
#ifdef USEHDF
    else Fhdf.close();
#endif

}
//@}

///\name Writes the hierarchy of structures
void WriteHierarchy(Options &opt, const Int_t &ngroups, const Int_t & nhierarchy, const Int_t &nfield, Int_t *nsub, Int_t *parentgid, Int_t *stype, int subflag){
    fstream Fout;
    fstream Fout2;
    char fname[500],fname2[500];
    unsigned long ng=ngroups,ngtot=0,noffset=0;
#ifdef USEHDF
    H5File Fhdf;
    H5std_string datasetname;
    DataSpace dataspace;
    DataSet dataset;
    hsize_t *dims;
    int rank;
    HDFCatalogNames hdfnames;
    int itemp=0;
#endif

    #ifdef USEMPI
    sprintf(fname,"%s.catalog_groups.%d",opt.outname,ThisTask);
#else
    int ThisTask=0,NProcs=1;
    sprintf(fname,"%s.catalog_groups",opt.outname);
#endif
    cout<<"saving hierarchy data to "<<fname<<endl;

    if (opt.ibinaryout==OUTBINARY) Fout.open(fname,ios::out|ios::binary|ios::app);
#ifdef USEHDF
    if (opt.ibinaryout==OUTHDF) {
        Fhdf=H5File(fname,H5F_ACC_RDWR);       
    }
#endif
    else Fout.open(fname,ios::out|ios::app);

    //since the hierarchy file is appended to the catalog_groups files, no header written
#ifdef USEMPI
    for (int j=0;j<NProcs;j++) ngtot+=mpi_ngroups[j];
    for (int j=0;j<ThisTask;j++)noffset+=mpi_ngroups[j];
#else
    ngtot=ngroups;
#endif

    //if subflag==0 only write number of substructures
    if (subflag==0) {
        if (opt.ibinaryout==OUTBINARY) Fout.write((char*)&nsub[1],sizeof(Int_t)*nfield);
#ifdef USEHDF
        else if (opt.ibinaryout==OUTHDF) {
            dims=new hsize_t[1];
            dims[0]=nfield;
            rank=1;
            itemp=4;
            //datasetname=H5std_string("Number_of_substructures_in_halo");
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            unsigned int *data=new unsigned int[nfield];
            for (Int_t i=1;i<=nfield;i++) data[i-1]=nsub[i];
            dataset.write(data,hdfnames.hierarchydatatype[itemp]);
            delete[] data;
            delete[] dims;
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
            dims=new hsize_t[1];
            dims[0]=ngroups-nfield;
            rank=1;
            itemp=4;
            //datasetname=H5std_string("Number_of_substructures_in_subhalo");
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            unsigned int *data=new unsigned int[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data[i-nfield-1]=nsub[i];
            dataset.write(data,hdfnames.hierarchydatatype[itemp]);
            delete[] data;
            itemp++;
            //datasetname=H5std_string("Parent_halo_ID");
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            long long *data2=new long long[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data2[i-nfield-1]=parentgid[i];
            dataset.write(data2,hdfnames.hierarchydatatype[itemp]);
            delete[] data2;
            delete[] dims;
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
            dims=new hsize_t[1];
            dims[0]=ngroups;
            rank=1;
            itemp=4;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            unsigned int *data=new unsigned int[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data[i-1]=nsub[i];
            dataset.write(data,hdfnames.hierarchydatatype[itemp]);
            delete[] data;
            itemp++;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            long long *data2=new long long[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data2[i-1]=parentgid[i];
            dataset.write(data2,hdfnames.hierarchydatatype[itemp]);
            delete[] data2;
            delete[] dims;
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
#ifdef USEMPI
    sprintf(fname,"%s.hierarchy.%d",opt.outname,ThisTask);
#else
    sprintf(fname,"%s.hierarchy",opt.outname);
#endif
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
        Fhdf=H5File(fname,H5F_ACC_TRUNC);
        //set file info
        dims=new hsize_t[1];
        dims[0]=1;
        rank=1;
        itemp=0;
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
        dataset.write(&ThisTask,hdfnames.hierarchydatatype[itemp]);
        itemp++;

        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
        dataset.write(&NProcs,hdfnames.hierarchydatatype[itemp]);
        itemp++;
        
        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
        dataset.write(&ng,hdfnames.hierarchydatatype[itemp]);
        itemp++;

        dataspace=DataSpace(rank,dims);
        dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
        dataset.write(&ngtot,hdfnames.hierarchydatatype[itemp]);
        itemp++;
        delete[] dims;
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
            dims=new hsize_t[1];
            dims[0]=nfield;
            rank=1;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            unsigned int *data=new unsigned int[nfield];
            for (Int_t i=1;i<=nfield;i++) data[i-1]=nsub[i];
            dataset.write(data,hdfnames.hierarchydatatype[itemp]);
            delete[] data;
            itemp++;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            long long *data2=new long long[nfield];
            for (Int_t i=1;i<=nfield;i++) data2[i-1]=parentgid[i];
            dataset.write(data2,hdfnames.hierarchydatatype[itemp]);
            delete[] data2;
            delete[] dims;
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
            dims=new hsize_t[1];
            dims[0]=ngroups-nfield;
            rank=1;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            unsigned int *data=new unsigned int[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data[i-nfield-1]=nsub[i];
            dataset.write(data,hdfnames.hierarchydatatype[itemp]);
            delete[] data;
            itemp++;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            long long *data2=new long long[ngroups-nfield];
            for (Int_t i=nfield+1;i<=ngroups;i++) data2[i-nfield-1]=parentgid[i];
            dataset.write(data2,hdfnames.hierarchydatatype[itemp]);
            delete[] data2;
            delete[] dims;
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
            dims=new hsize_t[1];
            dims[0]=ngroups;
            rank=1;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            unsigned int *data=new unsigned int[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data[i-1]=nsub[i];
            dataset.write(data,hdfnames.hierarchydatatype[itemp]);
            delete[] data;
            itemp++;
            dataspace=DataSpace(rank,dims);
            dataset = Fhdf.createDataSet(hdfnames.hierarchy[itemp], hdfnames.hierarchydatatype[itemp], dataspace);
            long long *data2=new long long[ngroups];
            for (Int_t i=1;i<=ngroups;i++) data2[i-1]=parentgid[i];
            dataset.write(data2,hdfnames.hierarchydatatype[itemp]);
            delete[] data2;
            delete[] dims;
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

    // Send and Receive Particles
    for (Int_t i = 1; i < NProcs; i++)
    {
        int src = (ThisTask + NProcs - i) % NProcs;
        int dst = (ThisTask + i) % NProcs;
        MPI_Isend (PartsToSend[dst], ntosendtotask[dst]*sizeof(Particle), MPI_BYTE, dst, ThisTask, MPI_COMM_WORLD, &rqst);
        MPI_Recv (PartsToRecv[src], ntorecievefromtask[src]*sizeof(Particle), MPI_BYTE, src, src, MPI_COMM_WORLD, &status);      
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
