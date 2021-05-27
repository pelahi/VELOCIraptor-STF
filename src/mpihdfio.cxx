/*! \file mpihdfio.cxx
 *  \brief this file contains routines used with MPI compilation and HDF io and domain construction.
 */

#if defined(USEMPI) && defined(USEHDF)

//-- For MPI

#include "logging.h"
#include "stf.h"
#include "hdfitems.h"

/// \name HDF Domain decomposition
//@{

/*!
    Determine the domain decomposition.\n
    Here the domains are constructured in data units
    only read tasks should call this routine. It is tricky to get appropriate load balancing and correct number of particles per processor.\n

    I could use recursive binary splitting like kd-tree along most spread axis till have appropriate number of volumes corresponding
    to number of processors. Or build a Peno-Hilbert space filling curve.

    BUT for identifying particles on another mpi domain, simple enough to split volume into rectangular cells.
*/

///Determine Domain for HDF input
void MPIDomainExtentHDF(Options &opt){
    char buf[2000];
    HDF_Group_Names hdf_gnames;
    HDF_Header hdf_header_info;
    hid_t Fhdf;
    float floatbuff;
    double doublebuff;
    //FloatType floattype;
    //PredType HDFREALTYPE(PredType::NATIVE_FLOAT);
    int ireaderror = 0;

    if (ThisTask==0) {
        if(opt.num_files>1) sprintf(buf,"%s.0.hdf5",opt.fname);
        else sprintf(buf,"%s.hdf5",opt.fname);

        //Try block to detect exceptions raised by any of the calls inside it
        //try
        {
            //turn off the auto-printing when failure occurs so that we can
            //handle the errors appropriately
            //Exception::dontPrint();

            //Open the specified file and the specified dataset in the file.
            Fhdf = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
            LOG(info) << "Loading HDF header info in header group: " << hdf_gnames.Header_name;

            if (opt.ihdfnameconvention == HDFSWIFTEAGLENAMES || opt.ihdfnameconvention == HDFOLDSWIFTEAGLENAMES)
            {
                /* SWIFT can have non-cubic boxes; but for cosmological runs they will always be cubes.
                * This makes the BoxSize a vector attribute, with it containing three values, but they
                * will always be the same. */
                hdf_header_info.BoxSize = read_attribute_v<double>(Fhdf, hdf_header_info.names[hdf_header_info.IBoxSize])[0];
            }
            else
            {
                hdf_header_info.BoxSize = read_attribute<double>(Fhdf, hdf_header_info.names[hdf_header_info.IBoxSize]);
            }
        }
        /*
        catch(GroupIException error)
        {
            HDF5PrintError(error);
        }
        // catch failure caused by the H5File operations
        catch( FileIException error )
        {
            HDF5PrintError(error);
        }
        // catch failure caused by the DataSet operations
        catch( DataSetIException error )
        {
            HDF5PrintError(error);
            ireaderror=1;
        }
        // catch failure caused by the DataSpace operations
        catch( DataSpaceIException error )
        {
            HDF5PrintError(error);
            ireaderror=1;
        }
        // catch failure caused by the DataSpace operations
        catch( DataTypeIException error )
        {
            HDF5PrintError(error);
            ireaderror=1;
        }
        Fhdf.close();
        */
        H5Fclose(Fhdf);
        for (int i=0;i<3;i++) {mpi_xlim[i][0]=0;mpi_xlim[i][1]=hdf_header_info.BoxSize;}
    }
    //There may be issues with particles exactly on the edge of a domain so before expanded limits by a small amount
    //now only done if a specific compile option passed
#ifdef MPIEXPANDLIM
    for (int i=0;i<3;i++) {
        Double_t dx=0.001*(mpi_xlim[i][1]-mpi_xlim[i][0]);
        mpi_xlim[i][0]-=dx;mpi_xlim[i][1]+=dx;
    }
#endif

    //make sure limits have been found
    MPI_Barrier(MPI_COMM_WORLD);
    if (NProcs==1) {
        for (int i=0;i<3;i++) {
            mpi_domain[ThisTask].bnd[i][0]=mpi_xlim[i][0];
            mpi_domain[ThisTask].bnd[i][1]=mpi_xlim[i][1];
        }
    }
}

void MPIDomainDecompositionHDF(Options &opt){
    Int_t i,j,k,n,m;
    int Nsplit,isplit;

    if (ThisTask==0) {
    }
}

///reads HDF file to determine number of particles in each MPIDomain
void MPINumInDomainHDF(Options &opt)
{
    if (NProcs==1) return;
    MPIDomainExtentHDF(opt);
    MPIInitialDomainDecomposition(opt);
    MPIDomainDecompositionHDF(opt);

    Int_t i,j,k;
    unsigned long long n,nchunk;
    char buf[2000];
    MPI_Status status;

    //structure stores the names of the groups in the hdf input
    HDF_Group_Names hdf_gnames (opt.ihdfnameconvention);
    //structures store names in groups
    HDF_Part_Info hdf_gas_info(HDFGASTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_dm_info(HDFDMTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_extradm_info(HDFDM1TYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_tracer_info(HDFTRACERTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_star_info(HDFSTARTYPE,opt.ihdfnameconvention);
    HDF_Part_Info hdf_bh_info(HDFBHTYPE,opt.ihdfnameconvention);

    HDF_Part_Info *hdf_parts[NHDFTYPE];
    hdf_parts[0]=&hdf_gas_info;
    hdf_parts[1]=&hdf_dm_info;
    #ifdef HIGHRES
    hdf_parts[2]=&hdf_extradm_info;
    hdf_parts[3]=&hdf_extradm_info;
    #else
    hdf_parts[2]=&hdf_extradm_info;
    hdf_parts[3]=&hdf_tracer_info;
    #endif
    hdf_parts[4]=&hdf_star_info;
    hdf_parts[5]=&hdf_bh_info;

    //to store the groups, data sets and their associated data spaces
    vector<HDF_Header> hdf_header_info;
    vector<hid_t>Fhdf;
    vector<hid_t>headerattribs;
    vector<hid_t>headerdataspace;
    vector<hid_t>partsgroup;
    vector<hid_t>partsdataset;
    vector<hid_t>partsdataspace;
    hid_t chunkspace;
    unsigned long long chunksize=opt.inputbufsize;
    //buffers to load data
    double *doublebuff=new double[chunksize*3];
    vector<int> vintbuff;
    vector<long long> vlongbuff;
    //arrays to store number of items to read and offsets when selecting hyperslabs
    //at most one needs a dimensionality of 13 for the tracer particles in Illustris
    hsize_t filespacecount[13],filespaceoffset[13];
    int ifloat,iint;
    int datarank;
    hsize_t datadim[5];
    Int_t Nlocalbuf,ibuf=0,*Nbuf, *Nbaryonbuf;
    int *ireadtask,*readtaskID;
    hid_t plist_id = H5P_DEFAULT;
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    std::vector<int> ireadfile(opt.num_files);
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);
#ifdef USEPARALLELHDF
    MPI_Comm mpi_comm_read;
    MPI_Comm mpi_comm_parallel_read;
    int ThisReadTask, NProcsReadTask, ThisParallelReadTask, NProcsParallelReadTask;
    if (opt.nsnapread > opt.num_files) {
        MPI_Comm_split(MPI_COMM_WORLD, (ireadtask[ThisTask]>=0), ThisTask, &mpi_comm_read);
        int ntaskread = ceil(opt.nsnapread/opt.num_files);
        int ifile = floor(ireadtask[ThisTask]/ntaskread);
        MPI_Comm_rank(mpi_comm_read, &ThisReadTask);
        MPI_Comm_size(mpi_comm_read, &NProcsReadTask);
        MPI_Comm_split(mpi_comm_read, ifile, ThisReadTask, &mpi_comm_parallel_read);
        MPI_Comm_rank(mpi_comm_parallel_read, &ThisParallelReadTask);
        MPI_Comm_size(mpi_comm_parallel_read, &NProcsParallelReadTask);
        plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, mpi_comm_parallel_read, MPI_INFO_NULL);
    }
#endif

    Nbuf=new Int_t[NProcs];
    Nbaryonbuf=new Int_t[NProcs];
    for (j=0;j<NProcs;j++) Nbuf[j]=0;
    for (j=0;j<NProcs;j++) Nbaryonbuf[j]=0;

    ///array listing number of particle types used.
    ///Since Illustris contains an unused type of particles (2) and tracer particles (3) really not useful to iterate over all particle types in loops
    int nusetypes,nbusetypes;
    int usetypes[NHDFTYPE];
    if (ireadtask[ThisTask]>=0) {
        HDFSetUsedParticleTypes(opt,nusetypes,nbusetypes,usetypes);
        hdf_header_info.resize(opt.num_files);
        Fhdf.resize(opt.num_files);
        headerdataspace.resize(opt.num_files);
        headerattribs.resize(opt.num_files);
        partsgroup.resize(opt.num_files*NHDFTYPE,-1);
        partsdataset.resize(opt.num_files*NHDFTYPE,-1);
        partsdataspace.resize(opt.num_files*NHDFTYPE,-1);

        MPISetFilesRead(opt,ireadfile,ireadtask);
        for(i=0; i<opt.num_files; i++) {
    	    if(ireadfile[i] == 0 ) continue;
            if(opt.num_files>1) sprintf(buf,"%s.%d.hdf5",opt.fname,int(i));
            else sprintf(buf,"%s.hdf5",opt.fname);
            //Open the specified file and the specified dataset in the file.
            Fhdf[i]=H5Fopen(buf, H5F_ACC_RDONLY, plist_id);
#ifdef USEPARALLELHDF
            H5Pclose(plist_id);
#endif
            //get number in file
            if (opt.ihdfnameconvention==HDFSWIFTEAGLENAMES || opt.ihdfnameconvention==HDFOLDSWIFTEAGLENAMES) {
                vlongbuff = read_attribute_v<long long>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INuminFile]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=vlongbuff[k];
            }
            else{
                vintbuff = read_attribute_v<int>(Fhdf[i], hdf_header_info[i].names[hdf_header_info[i].INuminFile]);
                for (k=0;k<NHDFTYPE;k++) hdf_header_info[i].npart[k]=vintbuff[k];
            }
            //open particle group structures
            for (j=0;j<nusetypes;j++) {k=usetypes[j]; partsgroup[i*NHDFTYPE+k]=HDF5OpenGroup(Fhdf[i],hdf_gnames.part_names[k]);}
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                 for (j=1;j<=nbusetypes;j++) {k=usetypes[j];partsgroup[i*NHDFTYPE+k]=HDF5OpenGroup(Fhdf[i],hdf_gnames.part_names[k]);}
            }

            //get positions
            for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[0]);
                partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                partsdataset[i*NHDFTYPE+k]=HDF5OpenDataSet(partsgroup[i*NHDFTYPE+k],hdf_parts[k]->names[0]);
                partsdataspace[i*NHDFTYPE+k]=HDF5OpenDataSpace(partsdataset[i*NHDFTYPE+k]);
            }
            for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                unsigned long long nstart = 0, nend = hdf_header_info[i].npart[k];
                unsigned long long nlocalsize;
#ifdef USEPARALLELHDF
                if (opt.num_files<opt.nsnapread) {
                    plist_id = H5Pcreate(H5P_DATASET_XFER);
                    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
                    nlocalsize = nend / NProcsParallelReadTask;
                    nstart = nlocalsize*ThisParallelReadTask;
                    if (ThisParallelReadTask < NProcsParallelReadTask -1)
                        nend = nlocalsize + nstart;
                }
#endif
                if (nend-nstart<chunksize)nchunk=nend-nstart;
                else nchunk=chunksize;
                for(n=nstart;n<nend;n+=nchunk)
                {
                    if (nend - n < chunksize && nend - n > 0) nchunk=nend-n;
                    //setup hyperslab so that it is loaded into the buffer
                    HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n, plist_id);
                    for (auto nn=0;nn<nchunk;nn++) {
                        ibuf=MPIGetParticlesProcessor(opt, doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Nbuf[ibuf]++;
                    }
                }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    unsigned long long nstart = 0, nend = hdf_header_info[i].npart[k];
#ifdef USEPARALLELHDF
                    if (opt.num_files<opt.nsnapread) {
                        unsigned long long nlocalsize = nend / NProcsParallelReadTask;
                        nstart = nlocalsize*ThisParallelReadTask;
                        if (ThisParallelReadTask < NProcsParallelReadTask -1)
                            nend = nlocalsize + nstart;
                    }
#endif
                    if (nend-nstart<chunksize)nchunk=nend-nstart;
                    else nchunk=chunksize;
                    for(n=nstart;n<nend;n+=nchunk)
                    {
                        if (nend - n < chunksize && nend - n > 0) nchunk=nend-n;
                        // setup hyperslab so that it is loaded into the buffer
                        HDF5ReadHyperSlabReal(doublebuff, partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n, plist_id);

                        for (auto nn=0;nn<nchunk;nn++) {
                            ibuf=MPIGetParticlesProcessor(opt, doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                            Nbaryonbuf[ibuf]++;
                        }
                    }
                }
            }
#ifdef USEPARALLELHDF
            H5Pclose(plist_id);
#endif
            //close data spaces
            for (auto &hidval:partsdataspace) HDF5CloseDataSpace(hidval);
            for (auto &hidval:partsdataset) HDF5CloseDataSet(hidval);
            for (auto &hidval:partsgroup) HDF5CloseGroup(hidval);
            H5Fclose(Fhdf[i]);
	   }
    }
    //now having read number of particles, run all gather
    Int_t mpi_nlocal[NProcs];
    MPI_Allreduce(Nbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
    Nlocal=mpi_nlocal[ThisTask];
    if (opt.iBaryonSearch) {
        MPI_Allreduce(Nbaryonbuf,mpi_nlocal,NProcs,MPI_Int_t,MPI_SUM,MPI_COMM_WORLD);
        Nlocalbaryon[0]=mpi_nlocal[ThisTask];
    }
#ifdef USEPARALLELHDF
    if (opt.nsnapread > opt.num_files) {
        MPI_Comm_free(&mpi_comm_parallel_read);
        MPI_Comm_free(&mpi_comm_read);
    }
#endif
    delete[] ireadtask;
    delete[] readtaskID;
    delete[] doublebuff;
    delete[] Nbuf;
    delete[] Nbaryonbuf;

}

//@}

#endif
