/*! \file mpihdfio.cxx
 *  \brief this file contains routines used with MPI compilation and HDF io and domain construction.
 */

#if defined(USEMPI) && defined(USEHDF)

//-- For MPI

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
            cout<<"Loading HDF header info in header group: "<<hdf_gnames.Header_name<<endl;

            hdf_header_info.BoxSize = read_attribute<double>(Fhdf, hdf_header_info.names[hdf_header_info.IBoxSize]);
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
    MPIInitialDomainDecomposition();
    MPIDomainDecompositionHDF(opt);

    Int_t i,j,k,n,nchunk;
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
    HDF_Header *hdf_header_info;
    vector<hid_t>Fhdf;
    vector<hid_t>headerattribs;
    vector<hid_t>headerdataspace;
    vector<hid_t>partsgroup;
    vector<hid_t>partsdataset;
    vector<hid_t>partsdataspace;
    hid_t chunkspace;
    int chunksize=opt.inputbufsize;
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
    int *ireadfile,*ireadtask,*readtaskID;
    ireadtask=new int[NProcs];
    readtaskID=new int[opt.nsnapread];
    ireadfile=new int[opt.num_files];
    MPIDistributeReadTasks(opt,ireadtask,readtaskID);

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
        hdf_header_info=new HDF_Header[opt.num_files];
        Fhdf.resize(opt.num_files);
        headerdataspace.resize(opt.num_files);
        headerattribs.resize(opt.num_files);
        partsgroup.resize(opt.num_files*NHDFTYPE);
        partsdataset.resize(opt.num_files*NHDFTYPE);
        partsdataspace.resize(opt.num_files*NHDFTYPE);

        MPISetFilesRead(opt,ireadfile,ireadtask);
        for(i=0; i<opt.num_files; i++) {
    	    if(ireadfile[i] == 0 ) continue;
            if(opt.num_files>1) sprintf(buf,"%s.%d.hdf5",opt.fname,i);
            else sprintf(buf,"%s.hdf5",opt.fname);
            //Open the specified file and the specified dataset in the file.
            // Fhdf[i].openFile(buf, H5F_ACC_RDONLY);
            Fhdf[i]=H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
            //get number in file
            if (opt.ihdfnameconvention==HDFSWIFTEAGLENAMES) {
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
                //data loaded into memory in chunks
                if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                else nchunk=chunksize;
                for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                {
                    if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                    //setup hyperslab so that it is loaded into the buffer
                    HDF5ReadHyperSlabReal(doublebuff,partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n);
                    for (int nn=0;nn<nchunk;nn++) {
                        ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                        Nbuf[ibuf]++;
                    }
                }
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) {
                for (j=1;j<=nbusetypes;j++) {
                    k=usetypes[j];
                    //data loaded into memory in chunks
                    if (hdf_header_info[i].npart[k]<chunksize)nchunk=hdf_header_info[i].npart[k];
                    else nchunk=chunksize;
                    for(n=0;n<hdf_header_info[i].npart[k];n+=nchunk)
                    {
                        if (hdf_header_info[i].npart[k]-n<chunksize&&hdf_header_info[i].npart[k]-n>0)nchunk=hdf_header_info[i].npart[k]-n;
                        // setup hyperslab so that it is loaded into the buffer
                        HDF5ReadHyperSlabReal(doublebuff, partsdataset[i*NHDFTYPE+k], partsdataspace[i*NHDFTYPE+k], 1, 3, nchunk, n);

                        for (int nn=0;nn<nchunk;nn++) {
                            ibuf=MPIGetParticlesProcessor(doublebuff[nn*3],doublebuff[nn*3+1],doublebuff[nn*3+2]);
                            Nbaryonbuf[ibuf]++;
                        }
                    }
                }
            }
            for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                HDF5CloseDataSpace(partsdataspace[i*NHDFTYPE+k]);
                HDF5CloseDataSet(partsdataset[i*NHDFTYPE+k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                HDF5CloseDataSpace(partsdataspace[i*NHDFTYPE+k]);
                HDF5CloseDataSet(partsdataset[i*NHDFTYPE+k]);
            }
            for (j=0;j<nusetypes;j++) {
                k=usetypes[j];
                HDF5CloseGroup(partsdataspace[i*NHDFTYPE+k]);
            }
            if (opt.partsearchtype==PSTDARK && opt.iBaryonSearch) for (j=1;j<=nbusetypes;j++) {
                k=usetypes[j];
                HDF5CloseGroup(partsdataspace[i*NHDFTYPE+k]);
            }
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
    delete[] doublebuff;

}

//@}

#endif
