/*! \file stfio.cxx
 *  \brief this file contains routines for reading stf style group catalog data
 */

#include "TreeFrog.h"

#ifdef USEHDF
#include "H5Cpp.h"
using namespace H5;
#endif

///\name routines for opening closing the set of files produced by velociraptor
//@{
void OpenBinaryorAsciiFiles(string &infile, int ibinary, int numfiletypes, int k, int mpi_ninput, int ifieldhalos, int itypematch,
    fstream &Fgroup, fstream &Fpart, fstream &Fupart,
    fstream &Fsgroup, fstream &Fspart, fstream &Fsupart,
    fstream &Fparttype, fstream &Fuparttype,
    fstream &Fsparttype, fstream &Fsuparttype,
    fstream &Fhaloinfo, fstream &Fshaloinfo,
    int iverbose)
{
    string fname1,fname2,fname3;
    string fname4,fname5,fname6;
    string fname7,fname8;
    string fname9,fname10;
    string fname11,fname12;
    string fnamearray[12];
    fstream *Farray[12];
    int itemp;
    VELOCIraptorFileTypeNames vftn;
    vftn.UpdateName(infile,fnamearray,k,mpi_ninput);
    fname1=fnamearray[0];
    fname2=fnamearray[1];
    fname3=fnamearray[2];
    fname4=fnamearray[3];
    fname5=fnamearray[4];
    fname6=fnamearray[5];
    fname7=fnamearray[6];
    fname8=fnamearray[7];
    fname9=fnamearray[8];
    fname10=fnamearray[9];
    fname11=fnamearray[10];
    fname12=fnamearray[11];

    //set names of files
    itemp=0;
    fnamearray[itemp++]=fname1;fnamearray[itemp++]=fname2;fnamearray[itemp++]=fname3;
    if (ifieldhalos) {
        fnamearray[itemp++]=fname4;fnamearray[itemp++]=fname5;fnamearray[itemp++]=fname6;
    }
    if (itypematch!=ALLTYPEMATCH) {
        fnamearray[itemp++]=fname7;fnamearray[itemp++]=fname8;
        if (ifieldhalos) {
            fnamearray[itemp++]=fname9;fnamearray[itemp++]=fname10;
        }
    }
    fnamearray[itemp++]=fname11;
    if (ibinary==INBINARY) {
        Fgroup.open(fname1.c_str(),ios::in|ios::binary);
        Fpart.open(fname2.c_str(),ios::in|ios::binary);
        Fupart.open(fname3.c_str(),ios::in|ios::binary);
        if (ifieldhalos) {
            Fsgroup.open(fname4.c_str(),ios::in|ios::binary);
            Fspart.open(fname5.c_str(),ios::in|ios::binary);
            Fsupart.open(fname6.c_str(),ios::in|ios::binary);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7.c_str(),ios::in|ios::binary);
            Fuparttype.open(fname8.c_str(),ios::in|ios::binary);
            if (ifieldhalos) {
                Fsparttype.open(fname9.c_str(),ios::in|ios::binary);
                Fsuparttype.open(fname10.c_str(),ios::in|ios::binary);
            }
        }
        Fhaloinfo.open(fname11.c_str(),ios::in|ios::binary);
        if (ifieldhalos) Fshaloinfo.open(fname12.c_str(),ios::in|ios::binary);
    }
    else {
        Fgroup.open(fname1.c_str(),ios::in);
        Fpart.open(fname2.c_str(),ios::in);
        Fupart.open(fname3.c_str(),ios::in);
        if (ifieldhalos) {
            Fsgroup.open(fname4.c_str(),ios::in);
            Fspart.open(fname5.c_str(),ios::in);
            Fsupart.open(fname6.c_str(),ios::in);
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.open(fname7.c_str(),ios::in);
            Fuparttype.open(fname8.c_str(),ios::in);
            if (ifieldhalos) {
                Fsparttype.open(fname9.c_str(),ios::in);
                Fsuparttype.open(fname10.c_str(),ios::in);
            }
        }
        Fhaloinfo.open(fname11.c_str(),ios::in);
        if (ifieldhalos) Fshaloinfo.open(fname12.c_str(),ios::in);
    }

    itemp=0;
    Farray[itemp++]=&Fgroup;Farray[itemp++]=&Fpart;Farray[itemp++]=&Fupart;
    if (ifieldhalos) {
        Farray[itemp++]=&Fsgroup;Farray[itemp++]=&Fspart;Farray[itemp++]=&Fsupart;
    }
    if (itypematch!=ALLTYPEMATCH) {
        Farray[itemp++]=&Fparttype;Farray[itemp++]=&Fuparttype;
        if (ifieldhalos) {
            Farray[itemp++]=&Fsparttype;Farray[itemp++]=&Fsuparttype;
        }
    }
    Farray[itemp++]=&Fhaloinfo;
    if (ifieldhalos) Farray[itemp++]=&Fshaloinfo;
    for (int i=0;i<numfiletypes;i++) {
        if(!Farray[i]->is_open()){
            cerr<<"can't open "<<fnamearray[i]<<endl;
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,9);
#else
            exit(9);
#endif
        }
        else if (iverbose>1) cout<<"open "<<fnamearray[i]<<endl;
    }
}
void CloseBinaryorAsciiFiles(
    fstream &Fgroup, fstream &Fpart, fstream &Fupart,
    fstream &Fsgroup, fstream &Fspart, fstream &Fsupart,
    fstream &Fparttype, fstream &Fuparttype,
    fstream &Fsparttype, fstream &Fsuparttype,
    fstream &Fhaloinfo,fstream &Fshaloinfo,
    int ifieldhalos, int itypematch
)
{
    Fgroup.close();
    Fpart.close();
    Fupart.close();
    if (ifieldhalos) {
        Fsgroup.close();
        Fspart.close();
        Fsupart.close();
    }
    if (itypematch!=ALLTYPEMATCH) {
        Fparttype.close();
        Fuparttype.close();
        if (ifieldhalos) {
            Fsparttype.close();
            Fsuparttype.close();
        }
    }
    Fhaloinfo.close();
    if (ifieldhalos) Fshaloinfo.close();
}

#ifdef USEHDF
void OpenHDFFiles(string &infile, int numfiletypes, int k, int mpi_ninput, int ifieldhalos, int itypematch,
    H5File &Fgroup, H5File &Fpart, H5File &Fupart,
    H5File &Fsgroup, H5File &Fspart, H5File &Fsupart,
    H5File &Fparttype, H5File &Fuparttype,
    H5File &Fsparttype, H5File &Fsuparttype,
    H5File &Fhaloinfo, H5File &Fshaloinfo,
    int iverbose)
{
    string fname1,fname2,fname3;
    string fname4,fname5,fname6;
    string fname7,fname8;
    string fname9,fname10;
    string fname11,fname12;
    string fnamearray[12];
    H5File *Farray[12];
    VELOCIraptorFileTypeNames vftn;
    int itemp;
    vftn.UpdateName(infile,fnamearray,k,mpi_ninput);
    fname1=fnamearray[0];
    fname2=fnamearray[1];
    fname3=fnamearray[2];
    fname4=fnamearray[3];
    fname5=fnamearray[4];
    fname6=fnamearray[5];
    fname7=fnamearray[6];
    fname8=fnamearray[7];
    fname9=fnamearray[8];
    fname10=fnamearray[9];
    fname11=fnamearray[10];
    fname12=fnamearray[11];

    //set names of files
    itemp=0;
    fnamearray[itemp++]=fname1;fnamearray[itemp++]=fname2;fnamearray[itemp++]=fname3;
    if (ifieldhalos) {
        fnamearray[itemp++]=fname4;fnamearray[itemp++]=fname5;fnamearray[itemp++]=fname6;
    }
    if (itypematch!=ALLTYPEMATCH) {
        fnamearray[itemp++]=fname7;fnamearray[itemp++]=fname8;
        if (ifieldhalos) {
            fnamearray[itemp++]=fname9;fnamearray[itemp++]=fname10;
        }
    }
    fnamearray[itemp++]=fname11;
    if (ifieldhalos) fnamearray[itemp++]=fname12;
    itemp=0;
    Farray[itemp++]=&Fgroup;Farray[itemp++]=&Fpart;Farray[itemp++]=&Fupart;
    if (ifieldhalos) {
        Farray[itemp++]=&Fsgroup;Farray[itemp++]=&Fspart;Farray[itemp++]=&Fsupart;
    }
    if (itypematch!=ALLTYPEMATCH) {
        Farray[itemp++]=&Fparttype;Farray[itemp++]=&Fuparttype;
        if (ifieldhalos) {
            Farray[itemp++]=&Fsparttype;Farray[itemp++]=&Fsuparttype;
        }
    }
    Farray[itemp++]=&Fhaloinfo;
    if (ifieldhalos) Farray[itemp++]=&Fshaloinfo;
    for (int i=0;i<numfiletypes;i++) {
        try {
            Farray[i]->openFile(fnamearray[i].c_str(),H5F_ACC_RDONLY);
            if (iverbose) cout<<"open "<<fnamearray[i]<<" "<<endl;
        }
        catch(const H5::FileIException&) {
            cerr<<"can't open "<<fnamearray[i]<<endl;
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD,9);
#else
            exit(9);
#endif
        }
    }
}
void CloseHDFFiles(
    H5File &Fgroup, H5File &Fpart, H5File &Fupart,
    H5File &Fsgroup, H5File &Fspart, H5File &Fsupart,
    H5File &Fparttype, H5File &Fuparttype,
    H5File &Fsparttype, H5File &Fsuparttype,
    H5File &Fhaloinfo, H5File &Fshaloinfo,
    int ifieldhalos, int itypematch
)
{
    Fgroup.close();
    Fpart.close();
    Fupart.close();
    if (ifieldhalos) {
        Fsgroup.close();
        Fspart.close();
        Fsupart.close();
    }
    if (itypematch!=ALLTYPEMATCH) {
        Fparttype.close();
        Fuparttype.close();
        if (ifieldhalos) {
            Fsparttype.close();
            Fsuparttype.close();
        }
    }
    Fhaloinfo.close();
    if (ifieldhalos) Fshaloinfo.close();
}
#endif
//@}

///\name helper routines to read information from velociraptor output files
//@{
///read information from the group catalog file and correct the number of files if necessary
inline void STFReadNumFileInfoAndCorrectNumFile(int &itask, int &nprocs, int &nmpicount, int &mpi_ninput, fstream &Fgroup, fstream &Fsgroup,
#ifdef USEHDF
    H5File &Fhdfgroup, DataSet &dataset, DataSpace &dataspace,
#endif
    int ibinary, int ifieldhalos)
{
#ifdef USEHDF
    HDFCatalogNames hdfnames;
#endif
    int itemp;

    if (ibinary==INBINARY) {
        Fgroup.read((char*)&itask,sizeof(int));
        Fgroup.read((char*)&nprocs,sizeof(int));
        if (ifieldhalos) {
            Fsgroup.read((char*)&itask,sizeof(int));
            Fsgroup.read((char*)&nprocs,sizeof(int));
        }
    }
#ifdef USEHDF
    else if (ibinary==INHDF) {
        itemp=0;
        dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&itask,hdfnames.groupdatatype[itemp],dataspace);
        itemp++;
        dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&nprocs,hdfnames.groupdatatype[itemp],dataspace);
    }
#endif
    else {
        Fgroup>>itask>>nprocs;
        if (ifieldhalos) {
            Fsgroup>>itask>>nprocs;
        }
    }
    if (nprocs!=mpi_ninput&&mpi_ninput>0) {
        nmpicount=mpi_ninput=nprocs;
    }
    else if ((mpi_ninput==0&&nprocs!=1)) {
        cout<<"Error, number of mpi outputs was set to zero but file indicates there are more than one mpi output"<<endl;
        cout<<"Terminating"<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,9);
#else
        exit(9);
#endif
    }
}


///read information from the group catalog file and correct the number of files if necessary
inline void STFReadNumGroups(unsigned long &nglocal, unsigned long &TotalNumberofHalos, unsigned long &nsglocal, fstream &Fgroup, fstream &Fsgroup,
#ifdef USEHDF
    H5File &Fhdfgroup, H5File &Fhdfsgroup, DataSet &dataset, DataSpace &dataspace,
#endif
    int ibinary, int ifieldhalos)
{
#ifdef USEHDF
    HDFCatalogNames hdfnames;
#endif
    int itemp;
    unsigned long ltemp;

    //now read data
    if (ibinary==INBINARY) {
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        nglocal=ltemp;
        Fgroup.read((char*)&ltemp,sizeof(unsigned long));
        TotalNumberofHalos=ltemp;
        if (ifieldhalos) {
            Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
            nsglocal=ltemp;
            Fsgroup.read((char*)&ltemp,sizeof(unsigned long));
        }
    }
#ifdef USEHDF
    else if (ibinary==INHDF) {
        //read from group file number of groups
        itemp=2;
        dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&ltemp,hdfnames.groupdatatype[itemp],dataspace);
        nglocal=ltemp;
        if (ifieldhalos) {
            dataset=Fhdfsgroup.openDataSet(hdfnames.group[itemp]);
            dataspace=dataset.getSpace();
            dataset.read(&ltemp,hdfnames.groupdatatype[itemp],dataspace);
            nsglocal=ltemp;
        }
        itemp++;

        dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&ltemp,hdfnames.groupdatatype[itemp],dataspace);
        TotalNumberofHalos=ltemp;
    }
#endif
    else {
        Fgroup>>nglocal;
        Fgroup>>TotalNumberofHalos;
        if (ifieldhalos) {
            Fsgroup>>nsglocal;
            Fsgroup>>itemp;
        }
    }
}

inline void STFReadNumData(unsigned long &nids, unsigned long &nsids, unsigned long &nuids, unsigned long &nsuids,
    unsigned long &nidstot, unsigned long &nsidstot, unsigned long &nuidstot, unsigned long &nsuidstot,
    fstream &Fpart, fstream &Fupart, fstream &Fspart, fstream &Fsupart,
    fstream &Fparttype, fstream &Fuparttype, fstream &Fsparttype, fstream &Fsuparttype,
#ifdef USEHDF
    H5File &Fhdfpart, H5File &Fhdfspart, H5File &Fhdfupart, H5File &Fhdfsupart,
    H5File &Fhdfparttype, H5File &Fhdfsparttype, H5File &Fhdfuparttype, H5File &Fhdfsuparttype,
    DataSet &dataset, DataSpace &dataspace,
#endif
    int ibinary, int ifieldhalos, int itypematch)
{

#ifdef USEHDF
    HDFCatalogNames hdfnames;
#endif
    int itemp;
    unsigned long ltemp;
    int itask,nprocs;
    nids=nsids=nuids=nsuids=nidstot=nsidstot=nuidstot=nsuidstot=0;
    //read the header info before reading lengths of the halos
    if (ibinary==INBINARY) {
        Fpart.read((char*)&itask,sizeof(int));
        Fpart.read((char*)&nprocs,sizeof(int));
        Fupart.read((char*)&itask,sizeof(int));
        Fupart.read((char*)&nprocs,sizeof(int));

        Fpart.read((char*)&nids,sizeof(unsigned long));
        Fpart.read((char*)&nidstot,sizeof(unsigned long));//total ids
        Fupart.read((char*)&nuids,sizeof(unsigned long));
        Fupart.read((char*)&nuidstot,sizeof(unsigned long));//total ids

        if (ifieldhalos) {
            Fspart.read((char*)&itask,sizeof(int));
            Fspart.read((char*)&nprocs,sizeof(int));
            Fsupart.read((char*)&itask,sizeof(int));
            Fsupart.read((char*)&nprocs,sizeof(int));

            Fspart.read((char*)&nsids,sizeof(unsigned long));
            Fspart.read((char*)&nsidstot,sizeof(unsigned long));//total ids
            Fsupart.read((char*)&nsuids,sizeof(unsigned long));
            Fsupart.read((char*)&nsuidstot,sizeof(unsigned long));//total ids
        }

        if (itypematch!=ALLTYPEMATCH) {
            Fparttype.read((char*)&itask,sizeof(int));
            Fparttype.read((char*)&nprocs,sizeof(int));
            Fuparttype.read((char*)&itask,sizeof(int));
            Fuparttype.read((char*)&nprocs,sizeof(int));

            Fparttype.read((char*)&nids,sizeof(unsigned long));
            Fparttype.read((char*)&nidstot,sizeof(unsigned long));//total ids
            Fuparttype.read((char*)&nuids,sizeof(unsigned long));
            Fuparttype.read((char*)&nuidstot,sizeof(unsigned long));//total ids

            if (ifieldhalos) {
                Fsparttype.read((char*)&itask,sizeof(int));
                Fsparttype.read((char*)&nprocs,sizeof(int));
                Fsuparttype.read((char*)&itask,sizeof(int));
                Fsuparttype.read((char*)&nprocs,sizeof(int));

                Fsparttype.read((char*)&nsids,sizeof(unsigned long));
                Fsparttype.read((char*)&nsidstot,sizeof(unsigned long));//total ids
                Fsuparttype.read((char*)&nsuids,sizeof(unsigned long));
                Fsuparttype.read((char*)&nsuidstot,sizeof(unsigned long));//total ids
            }
        }
    }
#ifdef USEHDF
    else if (ibinary==INHDF) {
        //if hdf do not have to offset headers at all
        //but still need to read number of ids, etc
        itemp=2;
        dataset=Fhdfpart.openDataSet(hdfnames.part[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
        nids=ltemp;
        dataset=Fhdfupart.openDataSet(hdfnames.part[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
        nuids=ltemp;
        itemp++;

        dataset=Fhdfpart.openDataSet(hdfnames.part[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
        nidstot=ltemp;
        dataset=Fhdfupart.openDataSet(hdfnames.part[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
        nuidstot=ltemp;

        if (ifieldhalos) {
            itemp=2;
            dataset=Fhdfspart.openDataSet(hdfnames.part[itemp]);
            dataspace=dataset.getSpace();
            dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
            nsids=ltemp;
            dataset=Fhdfsupart.openDataSet(hdfnames.part[itemp]);
            dataspace=dataset.getSpace();
            dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
            nsuids=ltemp;
            itemp++;

            dataset=Fhdfspart.openDataSet(hdfnames.part[itemp]);
            dataspace=dataset.getSpace();
            dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
            nsidstot=ltemp;
            dataset=Fhdfsupart.openDataSet(hdfnames.part[itemp]);
            dataspace=dataset.getSpace();
            dataset.read(&ltemp,hdfnames.partdatatype[itemp],dataspace);
            nsuidstot=ltemp;
        }
    }
#endif
    else {
        Fpart>>itask>>nprocs;
        Fupart>>itask>>nprocs;
        Fpart>>nids>>nidstot;
        Fupart>>nuids>>nuidstot;
        if (ifieldhalos) {
            Fspart>>itask>>nprocs;
            Fsupart>>itask>>nprocs;
            Fspart>>nsids>>nsidstot;
            Fsupart>>nsuids>>nsuidstot;
        }
        if (itypematch!=ALLTYPEMATCH) {
            Fparttype>>itask>>nprocs;
            Fuparttype>>itask>>nprocs;
            Fparttype>>nids>>nidstot;
            Fuparttype>>nuids>>nuidstot;
            if (ifieldhalos) {
                Fsparttype>>itask>>nprocs;
                Fsuparttype>>itask>>nprocs;
                Fsparttype>>nsids>>nsidstot;
                Fsuparttype>>nsuids>>nsuidstot;
            }
        }
    }
}


///read information from the group catalog file and correct the number of files if necessary
inline void STFReadHaloProperties(const unsigned long nglocal, HaloData *Halo, fstream &Fhaloinfo,
#ifdef USEHDF
    H5File &Fhdfhaloinfo, DataSet &dataset, DataSpace &dataspace,
#endif
    int ibinary)
{
#ifdef USEHDF
    HDFCatalogNames hdfnames;
#endif
    int itemp;
    unsigned long ltemp;
    void *data;
    //now read data
    if (ibinary==INBINARY) {
        //offset by header
        Fhaloinfo.read((char*)&itemp,sizeof(int));
        Fhaloinfo.read((char*)&itemp,sizeof(int));
        Fhaloinfo.read((char*)&itemp,sizeof(unsigned long));
        Fhaloinfo.read((char*)&ltemp,sizeof(unsigned long));
    }
#ifdef USEHDF
    else if (ibinary==INHDF) {
        data=::operator new(sizeof(double)*(nglocal+1));
        //read from halo properties, the positions, velocities, rmax and vmax of halos
        itemp=4;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Xcm[0]=((double*)data)[i];
        itemp++;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Xcm[1]=((double*)data)[i];
        itemp++;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Xcm[2]=((double*)data)[i];
        itemp++;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Vcm[0]=((double*)data)[i];
        itemp++;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Vcm[1]=((double*)data)[i];
        itemp++;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Vcm[2]=((double*)data)[i];
        itemp++;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Rmax=((double*)data)[i];
        itemp++;
        dataset=Fhdfhaloinfo.openDataSet(hdfnames.properties[itemp]);
        dataspace=dataset.getSpace();
        dataset.read(data,hdfnames.propertiesdatatype[itemp],dataspace);
        for (auto i=0;i<nglocal;i++) Halo[i].Vmax=((double*)data)[i];
        ::operator delete(data);
    }
#endif
    else {
        Fhaloinfo>>itemp;
        Fhaloinfo>>itemp;
        Fhaloinfo>>ltemp;
        Fhaloinfo>>ltemp;
    }
}
//@}

int CheckType(unsigned int t, int tmatch){
    if (tmatch>ALLTYPEMATCH) return (t==tmatch);
    else if (tmatch==DMGASTYPEMATCH) return (t==GASTYPEMATCH || t==DMTYPEMATCH);
    else return 0;
}

#ifdef USEMPI
///\name mpi related routines that open files for load balancing calculations
//@{
///Read halo data from an idividual snapshot;
unsigned long MPIReadHaloGroupCatalogDataNum(string &infile, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch)
{
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    fstream Fhaloinfo, Fshaloinfo; //halo properties
#ifdef USEHDF
    H5File Fhdfgroup,Fhdfpart,Fhdfupart; //field objects
    H5File Fhdfsgroup,Fhdfspart,Fhdfsupart; //sublevels
    H5File Fhdfparttype,Fhdfuparttype; //field objects
    H5File Fhdfsparttype,Fhdfsuparttype; //sublevels
    H5File Fhdfhaloinfo, Fhdfshaloinfo; //halo properties
    DataSpace dataspace;
    DataSet dataset;
    HDFCatalogNames hdfnames;
    void * data;
#endif
    //void pointer to
    void * Farray;
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }
    //halo properties
    numfiletypes+=1;
    if (ifieldhalos) numfiletypes+=1;

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;

    nglocal=nsglocal=0;
    //load the group catalogues and the related particle files
    //check if files exist and open
    if (ibinary!=INHDF) OpenBinaryorAsciiFiles(infile, ibinary, numfiletypes, 0, mpi_ninput, ifieldhalos, itypematch,
        Fgroup, Fpart, Fupart,
        Fsgroup, Fspart, Fsupart,
        Fparttype, Fuparttype, Fsparttype, Fsuparttype,
        Fhaloinfo, Fshaloinfo
    );
#ifdef USEHDF
    else OpenHDFFiles(infile, numfiletypes, 0, mpi_ninput, ifieldhalos, itypematch,
        Fhdfgroup, Fhdfpart, Fhdfupart,
        Fhdfsgroup, Fhdfspart, Fhdfsupart,
        Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
        Fhdfhaloinfo, Fhdfshaloinfo
    );
#endif

    //read number of output files
    STFReadNumFileInfoAndCorrectNumFile(itask, nprocs, nmpicount, mpi_ninput, Fgroup, Fsgroup,
#ifdef USEHDF
        Fhdfgroup, dataset, dataspace,
#endif
        ibinary, ifieldhalos);

    STFReadNumGroups(nglocal, TotalNumberofHalos, nsglocal, Fgroup, Fsgroup,
#ifdef USEHDF
        Fhdfgroup, Fhdfsgroup, dataset, dataspace,
#endif
        ibinary, ifieldhalos);
    if (ibinary!=INHDF) CloseBinaryorAsciiFiles(Fgroup, Fpart, Fupart,
        Fsgroup, Fspart, Fsupart,
        Fparttype, Fuparttype, Fsparttype, Fsuparttype,
        Fhaloinfo, Fshaloinfo,
        ifieldhalos, itypematch);
#ifdef USEHDF
    else CloseHDFFiles(Fhdfgroup, Fhdfpart, Fhdfupart,
        Fhdfsgroup, Fhdfspart, Fhdfsupart,
        Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
        Fhdfhaloinfo, Fhdfshaloinfo,
        ifieldhalos, itypematch);
#endif

    return TotalNumberofHalos;
}

///Read the number of particles in halos from an idividual snapshot, useful for splitting snaphots across threads in load balanced fashion
///\todo need to update the particle number to be just read from .catalog_particles and .catalog_particles.unbound directly
///instead of looping over files and loading arrays
unsigned long MPIReadHaloGroupCatalogDataParticleNum(string &infile, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch)
{

    HaloData *Halo;
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofParticles=0;
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    fstream Fhaloinfo, Fshaloinfo; //halo properties
#ifdef USEHDF
    H5File Fhdfgroup,Fhdfpart,Fhdfupart; //field objects
    H5File Fhdfsgroup,Fhdfspart,Fhdfsupart; //sublevels
    H5File Fhdfparttype,Fhdfuparttype; //field objects
    H5File Fhdfsparttype,Fhdfsuparttype; //sublevels
    H5File Fhdfhaloinfo, Fhdfshaloinfo; //halo properties
    DataSpace dataspace;
    DataSet dataset;
    HDFCatalogNames hdfnames;
    void * data;
#endif
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal,ngtot;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }
    numfiletypes+=1;
    if (ifieldhalos) numfiletypes+=1;

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    //load the group catalogues and the related particle files
    //check if files exist and open
    if (ibinary!=INHDF) OpenBinaryorAsciiFiles(infile, ibinary, numfiletypes, 0, mpi_ninput, ifieldhalos, itypematch,
        Fgroup, Fpart, Fupart,
        Fsgroup, Fspart, Fsupart,
        Fparttype, Fuparttype, Fsparttype, Fsuparttype,
        Fhaloinfo, Fshaloinfo
    );
#ifdef USEHDF
    else OpenHDFFiles(infile, numfiletypes, 0, mpi_ninput, ifieldhalos, itypematch,
        Fhdfgroup, Fhdfpart, Fhdfupart,
        Fhdfsgroup, Fhdfspart, Fhdfsupart,
        Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
        Fhdfhaloinfo, Fhdfshaloinfo
    );
#endif

    //read number of output files
    STFReadNumFileInfoAndCorrectNumFile(itask, nprocs, nmpicount, mpi_ninput, Fgroup, Fsgroup,
#ifdef USEHDF
        Fhdfgroup, dataset, dataspace,
#endif
        ibinary, ifieldhalos);

    STFReadNumData(nids, nsids, nuids, nsuids,
        nidstot, nsidstot, nuidstot, nsuidstot,
        Fpart, Fupart, Fspart, Fsupart, Fparttype, Fuparttype, Fsparttype, Fsuparttype,
#ifdef USEHDF
        Fhdfpart, Fhdfspart, Fhdfupart, Fhdfsupart,
        Fhdfparttype, Fhdfsparttype, Fhdfuparttype, Fhdfsuparttype,
        dataset, dataspace,
#endif
        ibinary, ifieldhalos, itypematch);
    TotalNumberofParticles=nidstot+nsidstot+nuidstot+nsuidstot;

    if (ibinary!=INHDF) CloseBinaryorAsciiFiles(Fgroup, Fpart, Fupart,
        Fsgroup, Fspart, Fsupart,
        Fparttype, Fuparttype, Fsparttype, Fsuparttype,
        Fhaloinfo, Fshaloinfo,
        ifieldhalos, itypematch);
#ifdef USEHDF
    else CloseHDFFiles(Fhdfgroup, Fhdfpart, Fhdfupart,
        Fhdfsgroup, Fhdfspart, Fhdfsupart,
        Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
        Fhdfhaloinfo, Fhdfshaloinfo,
        ifieldhalos, itypematch);
#endif

    return TotalNumberofParticles;
}

///Allocate memory to store halo information
HaloData *MPIReadHaloGroupCatalogDataAllocation(string &infile, Int_t &numhalos, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch)
{
    HaloData *Halo;
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    fstream Fhaloinfo, Fshaloinfo; //halo properties
#ifdef USEHDF
    H5File Fhdfgroup,Fhdfpart,Fhdfupart; //field objects
    H5File Fhdfsgroup,Fhdfspart,Fhdfsupart; //sublevels
    H5File Fhdfparttype,Fhdfuparttype; //field objects
    H5File Fhdfsparttype,Fhdfsuparttype; //sublevels
    H5File Fhdfhaloinfo, Fhdfshaloinfo; //halo properties
    DataSpace dataspace;
    DataSet dataset;
    HDFCatalogNames hdfnames;
    void * data;
#endif
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long ngtot;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }
    //halo properties
    numfiletypes+=1;
    if (ifieldhalos) numfiletypes+=1;

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (int k=0;k<nmpicount;k++) {
        nglocal=nsglocal=0;
        //load the group catalogues and the related particle files
        //check if files exist and open
        if (ibinary!=INHDF) OpenBinaryorAsciiFiles(infile, ibinary, numfiletypes, k, mpi_ninput, ifieldhalos, itypematch,
            Fgroup, Fpart, Fupart,
            Fsgroup, Fspart, Fsupart,
            Fparttype, Fuparttype, Fsparttype, Fsuparttype,
            Fhaloinfo, Fshaloinfo
        );
#ifdef USEHDF
        else OpenHDFFiles(infile, numfiletypes, k, mpi_ninput, ifieldhalos, itypematch,
            Fhdfgroup, Fhdfpart, Fhdfupart,
            Fhdfsgroup, Fhdfspart, Fhdfsupart,
            Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
            Fhdfhaloinfo, Fhdfshaloinfo
        );
#endif

        //read number of output files
        STFReadNumFileInfoAndCorrectNumFile(itask, nprocs, nmpicount, mpi_ninput, Fgroup, Fsgroup,
#ifdef USEHDF
            Fhdfgroup, dataset, dataspace,
#endif
            ibinary, ifieldhalos);

        STFReadNumGroups(nglocal, TotalNumberofHalos, nsglocal, Fgroup, Fsgroup,
#ifdef USEHDF
            Fhdfgroup, Fhdfsgroup, dataset, dataspace,
#endif
            ibinary, ifieldhalos);

        //allocate memory for halos
        if (k==0) {
            Halo=new HaloData[TotalNumberofHalos];
            numhalos=TotalNumberofHalos;
        }
        if (nglocal>0) {
            numingroup=new Int_t[nglocal+1];

            if (ibinary==INBINARY) {
                Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
                for (Int_t i=0;i<nglocal;i++) {
                    Halo[i+noffset].Alloc(numingroup[i]);
                }
            }
            else {
                for (Int_t i=0;i<nglocal;i++) {
                    Fgroup>>numingroup[i];
                    Halo[i+noffset].Alloc(numingroup[i]);
                }
            }

            delete[] numingroup;
        }

        if (ifieldhalos) {
            if (nsglocal>0) {
                //now read substructure data
                numingroup=new Int_t[nsglocal];
                if (ibinary==INBINARY) {
                    Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
                    for (Int_t i=0;i<nsglocal;i++) {
                        Halo[i+nglocal+noffset].Alloc(numingroup[i]);
                    }
                }
                else {
                    for (Int_t i=0;i<nsglocal;i++) {
                        Fsgroup>>numingroup[i];
                        Halo[i+nglocal+noffset].Alloc(numingroup[i]);
                    }
                }
                delete[] numingroup;
            }
        }
        noffset+=nglocal+nsglocal;

        if (ibinary!=INHDF) CloseBinaryorAsciiFiles(Fgroup, Fpart, Fupart,
            Fsgroup, Fspart, Fsupart,
            Fparttype, Fuparttype, Fsparttype, Fsuparttype,
            Fhaloinfo, Fshaloinfo,
            ifieldhalos, itypematch);
#ifdef USEHDF
        else CloseHDFFiles(Fhdfgroup, Fhdfpart, Fhdfupart,
            Fhdfsgroup, Fhdfspart, Fhdfsupart,
            Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
            Fhdfhaloinfo, Fhdfshaloinfo,
            ifieldhalos, itypematch);
#endif
    }

    return Halo;
}

///Read halo data from an individual snapshot;
void MPIReadHaloGroupCatalogData(string &infile, Int_t &numhalos, HaloData *&Halo, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch, int iverbose)
{
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    fstream Fhaloinfo, Fshaloinfo; //halo properties
#ifdef USEHDF
    H5File Fhdfgroup,Fhdfpart,Fhdfupart; //field objects
    H5File Fhdfsgroup,Fhdfspart,Fhdfsupart; //sublevels
    H5File Fhdfparttype,Fhdfuparttype; //field objects
    H5File Fhdfsparttype,Fhdfsuparttype; //sublevels
    H5File Fhdfhaloinfo, Fhdfshaloinfo; //halo properties
    DataSpace dataspace;
    DataSet dataset;
    HDFCatalogNames hdfnames;
    void * data;
#endif
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    float *X,*Y, *Z, *VX, *VY, *VZ, *Rmax, *Vmax;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }
    //halo properties
    numfiletypes+=1;
    if (ifieldhalos) numfiletypes+=1;

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (int k=0;k<nmpicount;k++) {
        nglocal=nsglocal=0;
        //load the group catalogues and the related particle files
        //check if files exist and open
        if (ibinary!=INHDF) OpenBinaryorAsciiFiles(infile, ibinary, numfiletypes, k, mpi_ninput, ifieldhalos, itypematch,
            Fgroup, Fpart, Fupart,
            Fsgroup, Fspart, Fsupart,
            Fparttype, Fuparttype, Fsparttype, Fsuparttype,
            Fhaloinfo, Fshaloinfo
        );
#ifdef USEHDF
        else OpenHDFFiles(infile, numfiletypes, k, mpi_ninput, ifieldhalos, itypematch,
            Fhdfgroup, Fhdfpart, Fhdfupart,
            Fhdfsgroup, Fhdfspart, Fhdfsupart,
            Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
            Fhdfhaloinfo, Fhdfshaloinfo
        );
#endif

        //read number of output files
        STFReadNumFileInfoAndCorrectNumFile(itask, nprocs, nmpicount, mpi_ninput, Fgroup, Fsgroup,
#ifdef USEHDF
            Fhdfgroup, dataset, dataspace,
#endif
            ibinary, ifieldhalos);

        STFReadNumGroups(nglocal, TotalNumberofHalos, nsglocal, Fgroup, Fsgroup,
#ifdef USEHDF
            Fhdfgroup, Fhdfsgroup, dataset, dataspace,
#endif
            ibinary, ifieldhalos);

        //allocate memory for halos
        if (k==0) {
            if (iverbose) cout<<infile<<" has "<<nglocal<<endl;
            if (iverbose) cout<<" and the total number of halos in all files is "<<TotalNumberofHalos<<endl;
            numhalos=TotalNumberofHalos;
        }

        STFReadNumData(nids, nsids, nuids, nsuids,
            nidstot, nsidstot, nuidstot, nsuidstot,
            Fpart, Fupart, Fspart, Fsupart,
            Fparttype, Fuparttype, Fsparttype, Fsuparttype,
#ifdef USEHDF
            Fhdfpart, Fhdfspart, Fhdfupart, Fhdfsupart,
            Fhdfparttype, Fhdfsparttype, Fhdfuparttype, Fhdfsuparttype,
            dataset, dataspace,
#endif
            ibinary, ifieldhalos, itypematch);

        if (nglocal>0) {
            numingroup=new Int_t[nglocal+1];
            offset=new Int_t[nglocal+1];
            uoffset=new Int_t[nglocal+1];
            numingroupbound=new Int_t[nglocal+1];
            idval=new Int_t[nids+nuids+1];
            if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];

            if (ibinary==INBINARY) {
                Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
                Fgroup.read((char*)offset,sizeof(Int_t)*nglocal);
                Fgroup.read((char*)uoffset,sizeof(Int_t)*nglocal);
            }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from group file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nglocal+1));
                dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                //group size is a unsigned int
                ///\todo should generalize type casting
                for (Int_t i=0;i<nglocal;i++) {numingroup[i]=((unsigned int*)data)[i];Halo[i+noffset].Alloc(numingroup[i]);}
                itemp++;

                dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nglocal;i++) offset[i]=((unsigned long long*)data)[i];
                itemp++;
                dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nglocal;i++) uoffset[i]=((unsigned long long*)data)[i];
                itemp++;
                ::operator delete(data);
            }
#endif
            else {
                for (Int_t i=0;i<nglocal;i++) Fgroup>>numingroup[i];
                for (Int_t i=0;i<nglocal;i++) Fgroup>>offset[i];
                for (Int_t i=0;i<nglocal;i++) Fgroup>>uoffset[i];
            }

            //now read bound particle list
            if (ibinary==INBINARY) {
                Fpart.read((char*)idval,sizeof(Int_t)*nids);
                Fupart.read((char*)&idval[nids],sizeof(Int_t)*nuids);
                if (itypematch!=ALLTYPEMATCH) {
                    Fparttype.read((char*)typeval,sizeof(UInt_t)*nids);
                    Fuparttype.read((char*)&typeval[nids],sizeof(UInt_t)*nuids);
                }
            }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from particle file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nids+nuids+1));

                dataset=Fhdfpart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nids;i++) idval[i]=((long long*)data)[i];
                dataset=Fhdfupart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nuids;i++) idval[i+nids]=((long long*)data)[i];
                if (itypematch!=ALLTYPEMATCH) {
                    dataset=Fhdfparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nids;i++) typeval[i]=((unsigned short*)data)[i];
                    dataset=Fhdfuparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nuids;i++) typeval[i+nids]=((unsigned short*)data)[i];
                }

                ::operator delete(data);
            }
#endif
            else {
                for (Int_t i=0;i<nids;i++) Fpart>>idval[i];
                for (Int_t i=0;i<nuids;i++) Fupart>>idval[i+nids];
                if (itypematch!=ALLTYPEMATCH) {
                    for (Int_t i=0;i<nids;i++) Fparttype>>typeval[i];
                    for (Int_t i=0;i<nuids;i++) Fuparttype>>typeval[i+nids];
                }
            }

            for (Int_t i=0;i<nglocal-1;i++) {
                numingroupbound[i]=Halo[i+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
            }
            numingroupbound[nglocal-1]=Halo[nglocal-1+noffset].NumberofParticles-(nuids-uoffset[nglocal-1]);
            //if don't care about particle type, simply store info
            if (itypematch==ALLTYPEMATCH) {
                for (Int_t i=0;i<nglocal;i++) {
                    for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+noffset].ParticleID[j]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) Halo[i+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nids];
                }
            }
            else {
                for (Int_t i=0;i<nglocal;i++) {
                    counter=0;
                    for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                    Halo[i+noffset].NumberofParticles=counter;
                }
            }
            delete[] idval;
            delete[] offset;
            delete[] uoffset;
            delete[] numingroupbound;
            delete[] numingroup;
            if (itypematch!=ALLTYPEMATCH) delete[] typeval;
            //now read positions if desired
            //alter interface so that it is passed a file and reads data
            STFReadHaloProperties(nglocal, &Halo[noffset], Fhaloinfo,
            #ifdef USEHDF
                Fhdfhaloinfo, dataset, dataspace,
            #endif
                ibinary);
        }

        for (Int_t i=0;i<nglocal+nsglocal;i++) Halo[i+noffset].haloID=i+1+noffset;
        if (ifieldhalos) {
            //cout<<infile<<" has sublevels with "<<nsglocal<<" "<<nsids<<" "<<nsuids<<endl;
            if (nsglocal>0) {
                //now read substructure data
                numingroup=new Int_t[nsglocal];
                offset=new Int_t[nsglocal];
                uoffset=new Int_t[nsglocal];
                numingroupbound=new Int_t[nsglocal];
                idval=new Int_t[nsids+nsuids];
                if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];
                if (ibinary==INBINARY) {
                    Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
                    Fsgroup.read((char*)offset,sizeof(Int_t)*nsglocal);
                    Fsgroup.read((char*)uoffset,sizeof(Int_t)*nsglocal);
                }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from group file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nsglocal+1));
                dataset=Fhdfsgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                //group size is a unsigned int
                ///\todo should generalize type casting
                for (Int_t i=0;i<nsglocal;i++) {numingroup[i]=((unsigned int*)data)[i];Halo[i+noffset].Alloc(numingroup[i]);}
                itemp++;

                dataset=Fhdfsgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsglocal;i++) offset[i]=((unsigned long long*)data)[i];
                itemp++;
                dataset=Fhdfsgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsglocal;i++) uoffset[i]=((unsigned long long*)data)[i];
                itemp++;
                ::operator delete(data);
            }
#endif
                else {
                    for (Int_t i=0;i<nsglocal;i++) Fsgroup>>numingroup[i];
                    for (Int_t i=0;i<nsglocal;i++) Fsgroup>>offset[i];
                    for (Int_t i=0;i<nsglocal;i++) Fsgroup>>uoffset[i];
                }
                //now read bound particle list
                if (ibinary==INBINARY) {
                    Fspart.read((char*)idval,sizeof(Int_t)*nsids);
                    Fsupart.read((char*)&idval[nsids],sizeof(Int_t)*nsuids);
                    if (itypematch!=ALLTYPEMATCH) {
                        Fsparttype.read((char*)typeval,sizeof(UInt_t)*nsids);
                        Fsuparttype.read((char*)&typeval[nsids],sizeof(UInt_t)*nsuids);
                    }
                }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from particle file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nsids+nsuids+1));

                dataset=Fhdfspart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsids;i++) idval[i]=((long long*)data)[i];
                dataset=Fhdfsupart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsuids;i++) idval[i+nids]=((long long*)data)[i];
                if (itypematch!=ALLTYPEMATCH) {
                    dataset=Fhdfsparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nsids;i++) typeval[i]=((unsigned short*)data)[i];
                    dataset=Fhdfsuparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nsuids;i++) typeval[i+nids]=((unsigned short*)data)[i];
                }

                ::operator delete(data);
            }
#endif
                else {
                    for (Int_t i=0;i<nsids;i++) Fspart>>idval[i];
                    for (Int_t i=0;i<nsuids;i++) Fsupart>>idval[i+nsids];
                    if (itypematch!=ALLTYPEMATCH) {
                        for (Int_t i=0;i<nsids;i++) Fsparttype>>typeval[i];
                        for (Int_t i=0;i<nsuids;i++) Fsuparttype>>typeval[i+nsids];
                    }
                }

                for (Int_t i=0;i<nsglocal-1;i++) {
                    numingroupbound[i]=Halo[i+nglocal+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
                }
                numingroupbound[nsglocal-1]=Halo[nglocal+nsglocal-1+noffset].NumberofParticles-(nsuids-uoffset[nsglocal-1]);
                if (itypematch==ALLTYPEMATCH) {
                    for (Int_t i=0;i<nsglocal;i++) {
                        for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+nglocal+noffset].ParticleID[j]=idval[offset[i]+j];
                        nn=numingroup[i]-numingroupbound[i];
                        for (Int_t j=0;j<nn;j++) Halo[i+nglocal+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nsids];
                    }
                }
                else {
                    for (Int_t i=0;i<nsglocal;i++) {
                        counter=0;
                        for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[offset[i]+j];
                        nn=numingroup[i]-numingroupbound[i];
                        for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                        Halo[i+nglocal+noffset].NumberofParticles=counter;
                    }
                }
                delete[] idval;
                delete[] offset;
                delete[] uoffset;
                delete[] numingroupbound;
                delete[] numingroup;
                if (itypematch!=ALLTYPEMATCH) delete[] typeval;
            }
            STFReadHaloProperties(nsglocal, &Halo[noffset+nglocal], Fshaloinfo,
            #ifdef USEHDF
                Fhdfshaloinfo, dataset, dataspace,
            #endif
                ibinary);
        }
        noffset+=nglocal+nsglocal;

        if (ibinary!=INHDF) CloseBinaryorAsciiFiles(Fgroup, Fpart, Fupart,
            Fsgroup, Fspart, Fsupart,
            Fparttype, Fuparttype, Fsparttype, Fsuparttype,
            Fhaloinfo, Fshaloinfo,
            ifieldhalos, itypematch);
#ifdef USEHDF
        else CloseHDFFiles(Fhdfgroup, Fhdfpart, Fhdfupart,
            Fhdfsgroup, Fhdfspart, Fhdfsupart,
            Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
            Fhdfhaloinfo, Fhdfshaloinfo,
            ifieldhalos, itypematch);
#endif

    }
}
//@}

#endif


///Read halo catalog data from an individual snapshot;
HaloData *ReadHaloGroupCatalogData(string &infile, Int_t &numhalos, int mpi_ninput, int ibinary, int ifieldhalos, int itypematch, int iverbose)
{
    HaloData *Halo;
    Int_t itemp;
    unsigned long ltemp;
    Int_t noffset,numfiletypes;
    long unsigned nparts,haloid;
    long unsigned TotalNumberofHalos;
    fstream Fgroup,Fpart,Fupart; //field objects
    fstream Fsgroup,Fspart,Fsupart; //sublevels
    fstream Fparttype,Fuparttype; //field objects
    fstream Fsparttype,Fsuparttype; //sublevels
    fstream Fhaloinfo, Fshaloinfo; //halo properties
#ifdef USEHDF
    H5File Fhdfgroup,Fhdfpart,Fhdfupart; //field objects
    H5File Fhdfsgroup,Fhdfspart,Fhdfsupart; //sublevels
    H5File Fhdfparttype,Fhdfuparttype; //field objects
    H5File Fhdfsparttype,Fhdfsuparttype; //sublevels
    H5File Fhdfhaloinfo, Fhdfshaloinfo; //halo properties
    DataSpace dataspace;
    DataSet dataset;
    HDFCatalogNames hdfnames;
    void * data;
#endif
    unsigned long nids,nuids,nsids,nsuids,nglocal,nsglocal;
    unsigned long nidstot,nuidstot,nsidstot,nsuidstot;
    Int_t *numingroup,*numingroupbound,*offset,*uoffset;
    Int_t counter,nn;
    int nmpicount,itask,nprocs;
    Int_t *idval;
    UInt_t *typeval;
#ifndef USEMPI
    int ThisTask=0,NProcs=1;
#endif
    numfiletypes=3;
    if (ifieldhalos) numfiletypes+=3;
    if (itypematch!=ALLTYPEMATCH) {
        numfiletypes+=2;
        if (ifieldhalos) numfiletypes+=2;
    }
    //halo properties
    numfiletypes+=1;
    if (ifieldhalos) numfiletypes+=1;

    ///adjust files read based on whether mpi output or not
    if (mpi_ninput==0) nmpicount=1;
    else nmpicount=mpi_ninput;
    //to offset the halo structure pointer given that data maybe split among multiple files
    noffset=0;
    for (int k=0;k<nmpicount;k++) {
        nglocal=nsglocal=0;
        //load the group catalogues and the related particle files
        if (ibinary!=INHDF) OpenBinaryorAsciiFiles(infile, ibinary, numfiletypes, k, mpi_ninput, ifieldhalos, itypematch,
            Fgroup, Fpart, Fupart,
            Fsgroup, Fspart, Fsupart,
            Fparttype, Fuparttype, Fsparttype, Fsuparttype,
            Fhaloinfo, Fshaloinfo
        );
#ifdef USEHDF
        else OpenHDFFiles(infile, numfiletypes, k, mpi_ninput, ifieldhalos, itypematch,
            Fhdfgroup, Fhdfpart, Fhdfupart,
            Fhdfsgroup, Fhdfspart, Fhdfsupart,
            Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
            Fhdfhaloinfo, Fhdfshaloinfo
        );
#endif
        //read number of output files
        STFReadNumFileInfoAndCorrectNumFile(itask, nprocs, nmpicount, mpi_ninput, Fgroup, Fsgroup,
#ifdef USEHDF
            Fhdfgroup, dataset, dataspace,
#endif
            ibinary, ifieldhalos);
        STFReadNumGroups(nglocal, TotalNumberofHalos, nsglocal, Fgroup, Fsgroup,
#ifdef USEHDF
            Fhdfgroup, Fhdfsgroup, dataset, dataspace,
#endif
            ibinary, ifieldhalos);

        //allocate memory for halos
        if (k==0) {
            if (iverbose) cout<<"reading "<<infile<<" split into "<<nmpicount<<" files "<<endl<<"Total number of halos in all files is "<<TotalNumberofHalos<<endl;
            Halo=new HaloData[TotalNumberofHalos];
            numhalos=TotalNumberofHalos;
            //if no haloes close file and return
            if (numhalos==0) {
                if (ibinary!=INHDF) CloseBinaryorAsciiFiles(Fgroup, Fpart, Fupart,
                    Fsgroup, Fspart, Fsupart,
                    Fparttype, Fuparttype, Fsparttype, Fsuparttype,
                    Fhaloinfo, Fshaloinfo,
                    ifieldhalos, itypematch);
#ifdef USEHDF
                else CloseHDFFiles(Fhdfgroup, Fhdfpart, Fhdfupart,
                    Fhdfsgroup, Fhdfspart, Fhdfsupart,
                    Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
                    Fhdfhaloinfo, Fhdfshaloinfo,
                    ifieldhalos, itypematch);
#endif
                return Halo;
            }
        }
        //if no local groups no need to read anything else in this file. So only proceed if ngloca>0
        if (nglocal>0) {
            STFReadNumData(nids, nsids, nuids, nsuids,
                nidstot, nsidstot, nuidstot, nsuidstot,
                Fpart, Fupart, Fspart, Fsupart,
                Fparttype, Fuparttype, Fsparttype, Fsuparttype,
#ifdef USEHDF
                Fhdfpart, Fhdfspart, Fhdfupart, Fhdfsupart,
                Fhdfparttype, Fhdfsparttype, Fhdfuparttype, Fhdfsuparttype,
                dataset, dataspace,
#endif
                ibinary, ifieldhalos, itypematch);
            //now files prepped to be read for relavant data

            numingroup=new Int_t[nglocal+1];
            offset=new Int_t[nglocal+1];
            uoffset=new Int_t[nglocal+1];
            numingroupbound=new Int_t[nglocal+1];
            idval=new Int_t[nids+nuids+1];
            if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];

            if (ibinary==INBINARY) {
                Fgroup.read((char*)numingroup,sizeof(Int_t)*nglocal);
                for (Int_t i=0;i<nglocal;i++) {
                    Halo[i+noffset].Alloc(numingroup[i]);
                }
                Fgroup.read((char*)offset,sizeof(Int_t)*nglocal);
                Fgroup.read((char*)uoffset,sizeof(Int_t)*nglocal);
            }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from group file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nglocal+1));
                dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                //group size is a unsigned int
                ///\todo should generalize type casting
                for (Int_t i=0;i<nglocal;i++) {numingroup[i]=((unsigned int*)data)[i];Halo[i+noffset].Alloc(numingroup[i]);}
                itemp++;

                dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nglocal;i++) offset[i]=((unsigned long long*)data)[i];
                itemp++;
                dataset=Fhdfgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nglocal;i++) uoffset[i]=((unsigned long long*)data)[i];
                itemp++;
                ::operator delete(data);
            }
#endif
            else {
                for (Int_t i=0;i<nglocal;i++) {
                    Fgroup>>numingroup[i];
                    Halo[i+noffset].Alloc(numingroup[i]);
                }
                for (Int_t i=0;i<nglocal;i++) Fgroup>>offset[i];
                for (Int_t i=0;i<nglocal;i++) Fgroup>>uoffset[i];
            }

            //now read bound particle list
            if (ibinary==INBINARY) {
                Fpart.read((char*)idval,sizeof(Int_t)*nids);
                Fupart.read((char*)&idval[nids],sizeof(Int_t)*nuids);
                if (itypematch!=ALLTYPEMATCH) {
                    Fparttype.read((char*)typeval,sizeof(UInt_t)*nids);
                    Fuparttype.read((char*)&typeval[nids],sizeof(UInt_t)*nuids);
                }
            }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from particle file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nids+nuids+1));

                dataset=Fhdfpart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nids;i++) idval[i]=((long long*)data)[i];
                dataset=Fhdfupart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nuids;i++) idval[i+nids]=((long long*)data)[i];
                if (itypematch!=ALLTYPEMATCH) {
                    dataset=Fhdfparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nids;i++) typeval[i]=((unsigned short*)data)[i];
                    dataset=Fhdfuparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nuids;i++) typeval[i+nids]=((unsigned short*)data)[i];
                }

                ::operator delete(data);
            }
#endif
            else {
                for (Int_t i=0;i<nids;i++) Fpart>>idval[i];
                for (Int_t i=0;i<nuids;i++) Fupart>>idval[i+nids];
                if (itypematch!=ALLTYPEMATCH) {
                    for (Int_t i=0;i<nids;i++) Fparttype>>typeval[i];
                    for (Int_t i=0;i<nuids;i++) Fuparttype>>typeval[i+nids];
                }
            }
            //now have loaded particle ids and group sizes

            for (Int_t i=0;i<nglocal-1;i++) {
                numingroupbound[i]=Halo[i+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
            }
            numingroupbound[nglocal-1]=Halo[nglocal-1+noffset].NumberofParticles-(nuids-uoffset[nglocal-1]);
            //if don't care about particle type, simply store info
            if (itypematch==ALLTYPEMATCH) {
                for (Int_t i=0;i<nglocal;i++) {
                    for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+noffset].ParticleID[j]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) Halo[i+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nids];
                }
            }
            else {
                for (Int_t i=0;i<nglocal;i++) {
                    counter=0;
                    for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                    Halo[i+noffset].NumberofParticles=counter;
                }
            }
            delete[] idval;
            delete[] offset;
            delete[] uoffset;
            delete[] numingroupbound;
            delete[] numingroup;
            if (itypematch!=ALLTYPEMATCH) delete[] typeval;
            STFReadHaloProperties(nglocal, &Halo[noffset], Fhaloinfo,
            #ifdef USEHDF
                Fhdfhaloinfo, dataset, dataspace,
            #endif
                ibinary);
        }

        for (Int_t i=0;i<nglocal+nsglocal;i++) Halo[i+noffset].haloID=i+1+noffset;
        //if field halos and subhaloes stored in separate files, proceed only if subhaloes present in local file.
        if (ifieldhalos && nsglocal) {
            //now read substructure data
            numingroup=new Int_t[nsglocal];
            offset=new Int_t[nsglocal];
            uoffset=new Int_t[nsglocal];
            numingroupbound=new Int_t[nsglocal];
            idval=new Int_t[nsids+nsuids];
            if (itypematch!=ALLTYPEMATCH) typeval=new UInt_t[nids+nuids+1];
            if (ibinary==INBINARY) {
                Fsgroup.read((char*)numingroup,sizeof(Int_t)*nsglocal);
                for (Int_t i=0;i<nsglocal;i++) {
                    Halo[i+nglocal+noffset].Alloc(numingroup[i]);
                }
                Fsgroup.read((char*)offset,sizeof(Int_t)*nsglocal);
                Fsgroup.read((char*)uoffset,sizeof(Int_t)*nsglocal);
            }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from group file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nsglocal+1));
                dataset=Fhdfsgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                //group size is a unsigned int
                ///\todo should generalize type casting
                for (Int_t i=0;i<nsglocal;i++) {numingroup[i]=((unsigned int*)data)[i];Halo[i+noffset].Alloc(numingroup[i]);}
                itemp++;

                dataset=Fhdfsgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsglocal;i++) offset[i]=((unsigned long long*)data)[i];
                itemp++;
                dataset=Fhdfsgroup.openDataSet(hdfnames.group[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.groupdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsglocal;i++) uoffset[i]=((unsigned long long*)data)[i];
                itemp++;
                ::operator delete(data);
            }
#endif
            else {
                for (Int_t i=0;i<nsglocal;i++) {
                    Fsgroup>>numingroup[i];
                    Halo[i+nglocal+noffset].Alloc(numingroup[i]);
                }
                for (Int_t i=0;i<nsglocal;i++) Fsgroup>>offset[i];
                for (Int_t i=0;i<nsglocal;i++) Fsgroup>>uoffset[i];
            }
            //now read bound particle list
            if (ibinary==INBINARY) {
                Fspart.read((char*)idval,sizeof(Int_t)*nsids);
                Fsupart.read((char*)&idval[nsids],sizeof(Int_t)*nsuids);
                if (itypematch!=ALLTYPEMATCH) {
                    Fsparttype.read((char*)typeval,sizeof(UInt_t)*nsids);
                    Fsuparttype.read((char*)&typeval[nsids],sizeof(UInt_t)*nsuids);
                }
            }
#ifdef USEHDF
            else if (ibinary==INHDF) {
                //read from particle file group size and offsets
                itemp=4;
                //allocate a void pointer large enough to store data
                data=::operator new(sizeof(long long)*(nsids+nsuids+1));

                dataset=Fhdfspart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsids;i++) idval[i]=((long long*)data)[i];
                dataset=Fhdfsupart.openDataSet(hdfnames.part[itemp]);
                dataspace=dataset.getSpace();
                dataset.read(data,hdfnames.partdatatype[itemp],dataspace);
                for (Int_t i=0;i<nsuids;i++) idval[i+nids]=((long long*)data)[i];
                if (itypematch!=ALLTYPEMATCH) {
                    dataset=Fhdfsparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nsids;i++) typeval[i]=((unsigned short*)data)[i];
                    dataset=Fhdfsuparttype.openDataSet(hdfnames.types[itemp]);
                    dataspace=dataset.getSpace();
                    dataset.read(data,hdfnames.typesdatatype[itemp],dataspace);
                    for (Int_t i=0;i<nsuids;i++) typeval[i+nids]=((unsigned short*)data)[i];
                }

                ::operator delete(data);
            }
#endif
            else {
                for (Int_t i=0;i<nsids;i++) Fspart>>idval[i];
                for (Int_t i=0;i<nsuids;i++) Fsupart>>idval[i+nsids];
                if (itypematch!=ALLTYPEMATCH) {
                    for (Int_t i=0;i<nsids;i++) Fsparttype>>typeval[i];
                    for (Int_t i=0;i<nsuids;i++) Fsuparttype>>typeval[i+nsids];
                }
            }

            for (Int_t i=0;i<nsglocal-1;i++) {
                numingroupbound[i]=Halo[i+nglocal+noffset].NumberofParticles-(uoffset[i+1]-uoffset[i]);
            }
            numingroupbound[nsglocal-1]=Halo[nglocal+nsglocal-1+noffset].NumberofParticles-(nsuids-uoffset[nsglocal-1]);
            if (itypematch==ALLTYPEMATCH) {
                for (Int_t i=0;i<nsglocal;i++) {
                    for (Int_t j=0;j<numingroupbound[i];j++) Halo[i+nglocal+noffset].ParticleID[j]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) Halo[i+nglocal+noffset].ParticleID[j+numingroupbound[i]]=idval[uoffset[i]+j+nsids];
                }
            }
            else {
                for (Int_t i=0;i<nsglocal;i++) {
                    counter=0;
                    for (Int_t j=0;j<numingroupbound[i];j++) if (CheckType(typeval[offset[i]+j],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[offset[i]+j];
                    nn=numingroup[i]-numingroupbound[i];
                    for (Int_t j=0;j<nn;j++) if (CheckType(typeval[uoffset[i]+j+nids],itypematch)) Halo[i+nglocal+noffset].ParticleID[counter++]=idval[uoffset[i]+j+nids];
                    Halo[i+nglocal+noffset].NumberofParticles=counter;
                }
            }
            delete[] idval;
            delete[] offset;
            delete[] uoffset;
            delete[] numingroupbound;
            delete[] numingroup;
            if (itypematch!=ALLTYPEMATCH) delete[] typeval;
            STFReadHaloProperties(nsglocal, &Halo[noffset+nglocal], Fshaloinfo,
            #ifdef USEHDF
                Fhdfshaloinfo, dataset, dataspace,
            #endif
                ibinary);
        }
        noffset+=nglocal+nsglocal;

        if (ibinary!=INHDF) CloseBinaryorAsciiFiles(Fgroup, Fpart, Fupart,
            Fsgroup, Fspart, Fsupart,
            Fparttype, Fuparttype, Fsparttype, Fsuparttype,
            Fhaloinfo, Fshaloinfo,
            ifieldhalos, itypematch);
#ifdef USEHDF
        else CloseHDFFiles(Fhdfgroup, Fhdfpart, Fhdfupart,
            Fhdfsgroup, Fhdfspart, Fhdfsupart,
            Fhdfparttype, Fhdfuparttype, Fhdfsparttype, Fhdfsuparttype,
            Fhdfhaloinfo, Fhdfshaloinfo,
            ifieldhalos, itypematch);
#endif

    }

    return Halo;
}
