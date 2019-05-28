/*! \file hdfitems.h
 *  \brief this file contains definitions and routines for reading HDF5 files
 *
 *   NOTE: the routines are based on reading Illustris HDF outputs
 */

#ifndef HDFITEMS_H
#define HDFITEMS_H


#include "H5Cpp.h"

using namespace H5;

///\name ILLUSTRIS specific constants
//@{
///convert illustris metallicty to ratio to solar
#define ILLUSTRISZMET 1.0/0.0127
//@}

///number of particle types
#define NHDFTYPE 6
///\name define hdf particle types
//@{
#define HDFGASTYPE 0
#define HDFDMTYPE 1
#define HDFWINDTYPE 2
#define HDFTRACERTYPE 3
#define HDFSTARTYPE 4
#define HDFBHTYPE 5
#define HDFEXTRATYPE 6

#define HDFDM1TYPE 2
#define HDFDM2TYPE 3
//@}

///\name number of entries in various data groups
//@{
#define HDFHEADNINFO 12
#define HDFGASNINFO 20
#define HDFDMNINFO 7
#define HDFTRACERNINFO 3
#define HDFSTARNINFO 13
#define HDFBHNINFO 21
#define HDFMAXPINFO 40
//@}

///\name to access where extra baryonic properties are located in the HDF_Part_Info structure that code will use to calcualte object properties
//@{
#define HDFGASIMETAL 0

#define HDFSTARIMETAL 40
#define HDFSTARIAGE 41

#define HDFBHIMDOT 50
//@}

///number of luminosity bands for stars
#define NLUMBANDS 8

///how many particle properties are read from file in one go
#define HDFCHUNKSIZE 100000

///how many data blocks are of  interest. That is what data blocks do I wish to load. pos,vel,mass,pid,U,Zmet,tage?
#define NHDFDATABLOCK 10
///here number shared by all particle types
#define NHDFDATABLOCKALL 4
//Maximum dimensionality of a datablock
///example at most one needs a dimensionality of 13 for the tracer particles in Illustris for fluid related info
#define HDFMAXPROPDIM 13

///\defgroup HDFNAMES labels for HDF naming conventions
//@{
#define HDFNUMNAMETYPES  8
#define HDFILLUSTISNAMES 0
#define HDFGADGETXNAMES  1
#define HDFEAGLENAMES    2
#define HDFGIZMONAMES    3
#define HDFSIMBANAMES    4
#define HDFMUFASANAMES   5
#define HDFSWIFTEAGLENAMES    6
#define HDFEAGLEVERSION2NAMES    7
//@}

///size of chunks in hdf files for Compression
#define HDFOUTPUTCHUNKSIZE 8192

#if H5_VERSION_GE(1,10,1)
#define HDF5_FILE_GROUP_COMMON_BASE H5::Group
#else
#define HDF5_FILE_GROUP_COMMON_BASE H5::CommonFG
#endif

template <typename AttributeHolder>
static inline
H5::Attribute get_attribute(const AttributeHolder &l, const std::string attr_name)
{
	auto exists = H5Aexists(l.getId(), attr_name.c_str());
	if (exists == 0) {
		throw invalid_argument(std::string("attribute not found ") + attr_name);
	}
	else if (exists < 0) {
		throw std::runtime_error("Error on H5Aexists");
	}
	return l.openAttribute(attr_name);
}

static inline
H5::Attribute get_attribute(const HDF5_FILE_GROUP_COMMON_BASE &file_or_group, const std::vector<std::string> &parts)
{
	// This is the attribute name
	if (parts.size() == 1) {
		return get_attribute(static_cast<const H5::Group &>(file_or_group), parts[0]);
	}

	auto n_groups = file_or_group.getNumObjs();

	const auto path = parts.front();
	for(hsize_t i = 0; i < n_groups; i++) {

		auto objname = file_or_group.getObjnameByIdx(i);
		if (objname != path) {
			continue;
		}

		auto objtype = file_or_group.getObjTypeByIdx(i);
		if (objtype == H5G_GROUP) {
			std::vector<std::string> subparts(parts.begin() + 1, parts.end());
			return get_attribute(file_or_group.openGroup(objname), subparts);
		}
		else if (objtype == H5G_DATASET) {
			std::vector<std::string> subparts(parts.begin() + 1, parts.end());
			return get_attribute(file_or_group.openDataSet(objname), parts.back());
		}
	}

	throw invalid_argument("attribute name not found");
}

static inline
vector<string> tokenize(const string &s, const string &delims)
{
	string::size_type lastPos = s.find_first_not_of(delims, 0);
	string::size_type pos     = s.find_first_of(delims, lastPos);

	vector<string> tokens;
	while (string::npos != pos || string::npos != lastPos) {
		tokens.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delims, pos);
		pos = s.find_first_of(delims, lastPos);
	}
	return tokens;
}

static inline
H5::Attribute get_attribute(const H5::H5File &file, const string &name)
{
	std::vector<std::string> parts = tokenize(name, "/");
	return get_attribute(file, parts);
}

template<typename T>
static inline
void _do_read(const H5::Attribute &attr, const H5::DataType type, T &val)
{
	attr.read(type, &val);
}

template<>
void _do_read<std::string>(const H5::Attribute &attr, const H5::DataType type, std::string &val)
{
	attr.read(type, val);
}

template<typename T>
const T read_attribute(const H5::H5File &filename, const std::string &name) {
	std::string attr_name;
	H5::Attribute attr = get_attribute(filename, name);
	H5::DataType type = attr.getDataType();
	T val;
	_do_read(attr, type, val);
	attr.close();
	return val;
}

template<typename T>
const T read_attribute(const std::string &filename, const std::string &name) {
	H5::H5File file(filename, H5F_ACC_RDONLY);
	return read_attribute<T>(file, name);
}

static inline void HDF5PrintError(const H5::Exception &error) {
#if H5_VERSION_GE(1,10,1)
	error.printErrorStack();
#else
	error.printError();
#endif
}

inline
H5::DataType _datatype_string(const std::string &val)
{
    return H5::StrType(H5::PredType::C_S1, val.size());
}

///This structures stores the strings defining the groups of data in the hdf input. NOTE: HERE I show the strings for Illustris format
struct HDF_Group_Names {
    //define the strings associated with the types of structures contained in the hdf file.
    H5std_string Header_name;
    H5std_string GASpart_name;
    H5std_string DMpart_name;
    H5std_string EXTRApart_name;
    H5std_string TRACERpart_name;
    H5std_string STARpart_name;
    H5std_string BHpart_name;
    H5std_string part_names[NHDFTYPE];
    H5std_string names[NHDFTYPE+1];

    ///constructor
    HDF_Group_Names(int hdfnametype=HDFEAGLENAMES){
        switch (hdfnametype) {
          case HDFSWIFTEAGLENAMES:
            Header_name=H5std_string("Header");
            GASpart_name=H5std_string("PartType0");
            DMpart_name=H5std_string("PartType1");
            EXTRApart_name=H5std_string("PartType2");
            TRACERpart_name=H5std_string("PartType3");
            STARpart_name=H5std_string("PartType4");
            BHpart_name=H5std_string("PartType5");
          break;

          default:
            Header_name=H5std_string("Header");
            GASpart_name=H5std_string("PartType0");
            DMpart_name=H5std_string("PartType1");
            EXTRApart_name=H5std_string("PartType2");
            TRACERpart_name=H5std_string("PartType3");
            STARpart_name=H5std_string("PartType4");
            BHpart_name=H5std_string("PartType5");
          break;
        }

        part_names[0]=GASpart_name;
        part_names[1]=DMpart_name;
        part_names[2]=EXTRApart_name;
        part_names[3]=TRACERpart_name;
        part_names[4]=STARpart_name;
        part_names[5]=BHpart_name;

        names[0]=Header_name;
        names[1]=GASpart_name;
        names[2]=DMpart_name;
        names[3]=EXTRApart_name;
        names[4]=TRACERpart_name;
        names[5]=STARpart_name;
        names[6]=BHpart_name;
    }
};

///data stored in the header group structure in the HDF format
struct HDF_Header {

    double      BoxSize;
    unsigned long int npart[NHDFTYPE];
    unsigned int npartTotal[NHDFTYPE];
    unsigned int npartTotalHW[NHDFTYPE];
    double      mass[NHDFTYPE];
    double      Omega0, OmegaLambda, HubbleParam;
    double      redshift, time;
    int         iscosmological;
    int         num_files;

    H5std_string names[HDFHEADNINFO];
    const static int IBoxSize  =0;
    const static int IMass     =1;
    const static int INuminFile=2;
    const static int INumTot   =3;
    const static int INumTotHW =4;
    const static int IOmega0   =5;
    const static int IOmegaL   =6;
    const static int IRedshift =7;
    const static int ITime     =8;
    const static int INumFiles =9;
    const static int IHubbleParam =10;
    const static int IIsCosmological =11;

    ///constructor
    HDF_Header(int hdfnametype=HDFEAGLENAMES) {
        int itemp=0;
        switch (hdfnametype) {
          case HDFSWIFTEAGLENAMES:
            names[itemp++]=H5std_string("Header/BoxSize");
            names[itemp++]=H5std_string("Header/MassTable");
            names[itemp++]=H5std_string("Header/NumPart_ThisFile");
            names[itemp++]=H5std_string("Header/NumPart_Total");
            names[itemp++]=H5std_string("Header/NumPart_Total_HighWord");
            names[itemp++]=H5std_string("Cosmology/Omega_m");
            names[itemp++]=H5std_string("Cosmology/Omega_lambda");
            names[itemp++]=H5std_string("Header/Redshift");
            names[itemp++]=H5std_string("Header/Time");
            names[itemp++]=H5std_string("Header/NumFilesPerSnapshot");
            names[itemp++]=H5std_string("Cosmology/h");
            names[itemp++]=H5std_string("Cosmology/Cosmological run");
            break;

          default:
            names[itemp++]=H5std_string("Header/BoxSize");
            names[itemp++]=H5std_string("Header/MassTable");
            names[itemp++]=H5std_string("Header/NumPart_ThisFile");
            names[itemp++]=H5std_string("Header/NumPart_Total");
            names[itemp++]=H5std_string("Header/NumPart_Total_HighWord");
            names[itemp++]=H5std_string("Header/Omega0");
            names[itemp++]=H5std_string("Header/OmegaLambda");
            names[itemp++]=H5std_string("Header/Redshift");
            names[itemp++]=H5std_string("Header/Time");
            names[itemp++]=H5std_string("Header/NumFilesPerSnapshot");
            names[itemp++]=H5std_string("Header/HubbleParam");
            break;
        }
    }
};

struct HDF_Part_Info {
    H5std_string names[HDFMAXPINFO];
    int ptype;
    int nentries;
    //store where properties are located
    int propindex[100];

    //the HDF naming convenction for the data blocks. By default assumes ILLUSTRIS nameing convention
    //for simplicity, all particles have basic properties listed first, x,v,ids,mass in this order
    HDF_Part_Info(int PTYPE, int hdfnametype=HDFEAGLENAMES) {
        ptype=PTYPE;
        int itemp=0;
        //gas
        if (ptype==HDFGASTYPE) {
            names[itemp++]=H5std_string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Velocity");
            else names[itemp++]=H5std_string("Velocities");
            names[itemp++]=H5std_string("ParticleIDs");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Mass");
            else names[itemp++]=H5std_string("Masses");
            names[itemp++]=H5std_string("Density");
            names[itemp++]=H5std_string("InternalEnergy");
            names[itemp++]=H5std_string("StarFormationRate");
            //always place the metacallity at position 7 in naming array
            if (hdfnametype==HDFILLUSTISNAMES) {
                propindex[HDFGASIMETAL]=itemp;
                names[itemp++]=H5std_string("GFM_Metallicity");
                names[itemp++]=H5std_string("ElectronAbundance");
                names[itemp++]=H5std_string("NeutralHydrogenAbundance");
                names[itemp++]=H5std_string("Volume");
                names[itemp++]=H5std_string("SmoothingLength");
                names[itemp++]=H5std_string("Potential");
                names[itemp++]=H5std_string("SubfindDensity");
                names[itemp++]=H5std_string("SubfindHsml");
                names[itemp++]=H5std_string("SubfindVelDisp");
                names[itemp++]=H5std_string("GFM_AGNRadiation");
                names[itemp++]=H5std_string("GFM_CoolingRate");
                names[itemp++]=H5std_string("GFM_WindDMVelDisp");
                names[itemp++]=H5std_string("NumTracers");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                propindex[HDFGASIMETAL]=itemp;
                //names[itemp++]=H5std_string("Metallicity");//11 metals stored in this data set
                names[itemp++]=H5std_string("Metallicity_00");//only grab the first of the 11, which is total
                names[itemp++]=H5std_string("ElectronAbundance");
                names[itemp++]=H5std_string("FractionH2");
                names[itemp++]=H5std_string("GrackleHI");
                names[itemp++]=H5std_string("GrackleHII");
                names[itemp++]=H5std_string("GrackleHM");
                names[itemp++]=H5std_string("GrackleHeI");
                names[itemp++]=H5std_string("GrackleHeII");
                names[itemp++]=H5std_string("GrackleHeIII");
                names[itemp++]=H5std_string("NWindLaunches");
                names[itemp++]=H5std_string("NeutralHydrogenAbundance");
                names[itemp++]=H5std_string("ParticleChildIDsNumber");
                names[itemp++]=H5std_string("ParticleIDGenerationNumber");
                names[itemp++]=H5std_string("Potential");
                names[itemp++]=H5std_string("Sigma");
                names[itemp++]=H5std_string("SmoothingLength");
            }
            else if (hdfnametype==HDFSIMBANAMES || hdfnametype==HDFMUFASANAMES) {
                propindex[HDFGASIMETAL]=itemp;
                names[itemp++]=H5std_string("Metallicity");//11 metals stored in this data set
                names[itemp++]=H5std_string("ElectronAbundance");
                names[itemp++]=H5std_string("FractionH2");
                names[itemp++]=H5std_string("GrackleHI");
                names[itemp++]=H5std_string("GrackleHII");
                names[itemp++]=H5std_string("GrackleHM");
                names[itemp++]=H5std_string("GrackleHeI");
                names[itemp++]=H5std_string("GrackleHeII");
                names[itemp++]=H5std_string("GrackleHeIII");
                names[itemp++]=H5std_string("NWindLaunches");
                names[itemp++]=H5std_string("NeutralHydrogenAbundance");
                names[itemp++]=H5std_string("Potential");
                names[itemp++]=H5std_string("Sigma");
                names[itemp++]=H5std_string("SmoothingLength");
                names[itemp++]=H5std_string("DelayTime");
                names[itemp++]=H5std_string("Dust_Masses");
                names[itemp++]=H5std_string("Dust_Metallicity");//11 metals stored in this data set
            }
            else if (hdfnametype==HDFEAGLENAMES) {
                propindex[HDFGASIMETAL]=itemp;
                names[itemp++]=H5std_string("Metallicity");
            }
        }
        //dark matter
        if (ptype==HDFDMTYPE) {
            names[itemp++]=H5std_string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Velocity");
            else names[itemp++]=H5std_string("Velocities");
            names[itemp++]=H5std_string("ParticleIDs");
            if (hdfnametype==HDFSWIFTEAGLENAMES) {
                names[itemp++]=H5std_string("Masses");
            }
            if (hdfnametype==HDFSIMBANAMES||hdfnametype==HDFMUFASANAMES) {
                names[itemp++]=H5std_string("Masses");
                names[itemp++]=H5std_string("Potential");
            }
            if (hdfnametype==HDFILLUSTISNAMES) {
                names[itemp++]=H5std_string("Potential");
                names[itemp++]=H5std_string("SubfindDensity");
                names[itemp++]=H5std_string("SubfindHsml");
                names[itemp++]=H5std_string("SubfindVelDisp");
            }
        }
        //also dark matter particles
        if (ptype==HDFDM1TYPE ||ptype==HDFDM2TYPE) {
            names[itemp++]=H5std_string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Velocity");
            else names[itemp++]=H5std_string("Velocities");
            names[itemp++]=H5std_string("ParticleIDs");
            names[itemp++]=H5std_string("Masses");
            if (hdfnametype==HDFSIMBANAMES||hdfnametype==HDFMUFASANAMES) {
                names[itemp++]=H5std_string("Potential");
            }
        }
        if (ptype==HDFTRACERTYPE) {
            names[itemp++]=H5std_string("FluidQuantities");
            names[itemp++]=H5std_string("ParentID");
            names[itemp++]=H5std_string("TracerID");
        }
        if (ptype==HDFSTARTYPE) {
            names[itemp++]=H5std_string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Velocity");
            else names[itemp++]=H5std_string("Velocities");
            names[itemp++]=H5std_string("ParticleIDs");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Mass");
            else names[itemp++]=H5std_string("Masses");
            //for stars assume star formation and metallicy are position 4, 5 in name array
            if (hdfnametype==HDFILLUSTISNAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=H5std_string("GFM_StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                names[itemp++]=H5std_string("GFM_Metallicity");
                names[itemp++]=H5std_string("Potential");
                names[itemp++]=H5std_string("SubfindDensity");
                names[itemp++]=H5std_string("SubfindHsml");
                names[itemp++]=H5std_string("SubfindVelDisp");
                names[itemp++]=H5std_string("GFM_InitialMass");
                names[itemp++]=H5std_string("GFM_StellarPhotometrics");
                names[itemp++]=H5std_string("NumTracers");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=H5std_string("StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                //names[itemp++]=H5std_string("Metallicity");//11 metals stored in this data set
                names[itemp++]=H5std_string("Metallicity_00");//only grab the first of the 11, which is total
                names[itemp++]=H5std_string("AGS-Softening Dataset");
                names[itemp++]=H5std_string("ParticleChildIDsNumber");
                names[itemp++]=H5std_string("ParticleIDGenerationNumber");
                names[itemp++]=H5std_string("Potential");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=H5std_string("StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                names[itemp++]=H5std_string("Metallicity");
                names[itemp++]=H5std_string("Potential");
                names[itemp++]=H5std_string("Dust_Masses");
                names[itemp++]=H5std_string("Dust_Metallicity");//11 metals stored in this data set
            }
            else if (hdfnametype==HDFEAGLENAMES) {
                propindex[HDFSTARIAGE]=itemp;
                names[itemp++]=H5std_string("StellarFormationTime");
                propindex[HDFSTARIMETAL]=itemp;
                names[itemp++]=H5std_string("Metallicity");
            }
        }
        if (ptype==HDFBHTYPE) {
            names[itemp++]=H5std_string("Coordinates");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Velocity");
            else names[itemp++]=H5std_string("Velocities");
            names[itemp++]=H5std_string("ParticleIDs");
            if(hdfnametype==HDFEAGLENAMES) names[itemp++]=H5std_string("Mass");
            else names[itemp++]=H5std_string("Masses");
            if (hdfnametype==HDFILLUSTISNAMES) {
                names[itemp++]=H5std_string("HostHaloMass");
                names[itemp++]=H5std_string("Potential");
                names[itemp++]=H5std_string("SubfindDensity");
                names[itemp++]=H5std_string("SubfindHsml");
                names[itemp++]=H5std_string("SubfindVelDisp");
                names[itemp++]=H5std_string("BH_CumEgyInjection_QM");
                names[itemp++]=H5std_string("BH_CumMassGrowth_QM");
                names[itemp++]=H5std_string("BH_Density");
                names[itemp++]=H5std_string("BH_Hsml");
                names[itemp++]=H5std_string("BH_Mass");
                names[itemp++]=H5std_string("BH_Mass_bubbles");
                names[itemp++]=H5std_string("BH_Mass_ini");
                names[itemp++]=H5std_string("BH_Mdot");
                names[itemp++]=H5std_string("BH_Pressure");
                names[itemp++]=H5std_string("BH_Progs");
                names[itemp++]=H5std_string("BH_U");
                names[itemp++]=H5std_string("NumTracers");
            }
            else if (hdfnametype==HDFGIZMONAMES) {
                names[itemp++]=H5std_string("AGS-Softening");
                names[itemp++]=H5std_string("Potential");
            }
            else if (hdfnametype==HDFMUFASANAMES) {
                names[itemp++]=H5std_string("StellarFormationTime");
                names[itemp++]=H5std_string("BH_AccretionLength");
                names[itemp++]=H5std_string("BH_Mass");
                names[itemp++]=H5std_string("BH_Mass_AlphaDisk");
                names[itemp++]=H5std_string("BH_Mass_Mdot");
                names[itemp++]=H5std_string("BH_NProgs");
                names[itemp++]=H5std_string("Potential");
            }
            else if (hdfnametype==HDFEAGLENAMES) {
                //names[itemp++]=H5std_string("StellarFormationTime");
                //names[itemp++]=H5std_string("Metallicity");
            }
        }
        nentries=itemp;
    }
};
//@}

/// \name Get the number of particles in the hdf files
//@{
inline Int_t HDF_get_nbodies(char *fname, int ptype, Options &opt)
{
    char buf[2000],buf1[2000],buf2[2000];
    sprintf(buf1,"%s.0.hdf5",fname);
    sprintf(buf2,"%s.hdf5",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }

    H5File Fhdf;
    HDF_Group_Names hdf_gnames;
    //to store the groups, data sets and their associated data spaces
    Attribute headerattribs;
    HDF_Header hdf_header_info = HDF_Header(opt.ihdfnameconvention);
    //buffers to load data
    string stringbuff;
    string swift_str = "SWIFT";
    int intbuff[NHDFTYPE];
    long long longbuff[NHDFTYPE];
    unsigned int uintbuff[NHDFTYPE];
    int j,k,ireaderror=0;
    Int_t nbodies=0;
    DataSpace headerdataspace;

    //to determine types
    IntType inttype;
    StrType stringtype;
    int nusetypes,usetypes[NHDFTYPE];

    if (ptype==PSTALL) {
        //lets assume there are dm/gas.
        nusetypes=0;
        usetypes[nusetypes++]=HDFGASTYPE;usetypes[nusetypes++]=HDFDMTYPE;
        if (opt.iuseextradarkparticles) {usetypes[nusetypes++]=HDFDM1TYPE;usetypes[nusetypes++]=HDFDM2TYPE;}
        if (opt.iusestarparticles) usetypes[nusetypes++]=HDFSTARTYPE;
        if (opt.iusesinkparticles) usetypes[nusetypes++]=HDFBHTYPE;
        if (opt.iusewindparticles) usetypes[nusetypes++]=HDFWINDTYPE;
        if (opt.iusetracerparticles) usetypes[nusetypes++]=HDFTRACERTYPE;
    }
    else if (ptype==PSTDARK) {
        nusetypes=1;usetypes[0]=HDFDMTYPE;
        if (opt.iuseextradarkparticles) usetypes[nusetypes++]=HDFDM1TYPE;usetypes[nusetypes++]=HDFDM2TYPE;
    }
    else if (ptype==PSTGAS) {nusetypes=1;usetypes[0]=HDFGASTYPE;}
    else if (ptype==PSTSTAR) {nusetypes=1;usetypes[0]=HDFSTARTYPE;}
    else if (ptype==PSTBH) {nusetypes=1;usetypes[0]=HDFBHTYPE;}
    //else if (ptype==PSTNOBH) {nusetypes=3;usetypes[0]=0;usetypes[1]=1;usetypes[2]=4;}

    //Try block to detect exceptions raised by any of the calls inside it
    try
    {
        //turn off the auto-printing when failure occurs so that we can
        //handle the errors appropriately
        Exception::dontPrint();

        //Open the specified file and the specified dataset in the file.
        Fhdf.openFile(buf, H5F_ACC_RDONLY);
        cout<<"Loading HDF header info in header group: "<<hdf_gnames.Header_name<<endl;

        if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) {

          // Check if it is a SWIFT snapshot.
          headerattribs=get_attribute(Fhdf, "Header/Code");
          stringtype = headerattribs.getStrType();
          headerattribs.read(stringtype, stringbuff);

          // Read SWIFT parameters
          if(!swift_str.compare(stringbuff)) {
            // Is it a cosmological simulation?
            headerattribs=get_attribute(Fhdf, hdf_header_info.names[hdf_header_info.IIsCosmological]);
            headerdataspace=headerattribs.getSpace();

            if (headerdataspace.getSimpleExtentNdims()!=1) ireaderror=1;
            inttype=headerattribs.getIntType();
            if (inttype.getSize()==sizeof(int)) {
              headerattribs.read(PredType::NATIVE_INT,&intbuff[0]);
              hdf_header_info.iscosmological = intbuff[0];
            }
            if (inttype.getSize()==sizeof(long long)) {
              headerattribs.read(PredType::NATIVE_LONG,&longbuff[0]);
              hdf_header_info.iscosmological = longbuff[0];
            }

            if (!hdf_header_info.iscosmological && opt.icosmologicalin) {
              cout<<"Error: cosmology is turned on in the config file but the snaphot provided is a non-cosmological run."<<endl;
#ifdef USEMPI
              MPI_Abort(MPI_COMM_WORLD, 8);
#else
              exit(0);
#endif

            }
            else if (hdf_header_info.iscosmological && !opt.icosmologicalin) {
              cout<<"Error: cosmology is turned off in the config file but the snaphot provided is a cosmological run."<<endl;
#ifdef USEMPI
              MPI_Abort(MPI_COMM_WORLD, 8);
#else
              exit(0);
#endif

            }
          }
          // If the code is not SWIFT
          else {
            cout<<"SWIFT EAGLE HDF5 naming convention chosen in config file but the snapshot was not produced by SWIFT. The string read was: "<<stringbuff<<endl;
#ifdef USEMPI
            MPI_Abort(MPI_COMM_WORLD, 8);
#else
            exit(0);
#endif
          }
        }

        headerattribs=get_attribute(Fhdf, hdf_header_info.names[hdf_header_info.INumTot]);

        headerattribs.read(PredType::NATIVE_UINT,&uintbuff);
        for (j=0;j<NHDFTYPE;j++) hdf_header_info.npartTotal[j]=uintbuff[j];

        headerattribs=get_attribute(Fhdf, hdf_header_info.names[hdf_header_info.INumTotHW]);
        headerattribs.read(PredType::NATIVE_UINT,&uintbuff);
        for (j=0;j<NHDFTYPE;j++) hdf_header_info.npartTotalHW[j]=uintbuff[j];
    }
    catch(GroupIException &error)
    {
        HDF5PrintError(error);
        cerr<<"Error in group might suggest config file has the incorrect HDF naming convention. ";
        cerr<<"Check HDF_name_convetion or add new naming convention updating hdfitems.h in the source code. "<<endl;
        Fhdf.close();
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD,8);
#else
        exit(8);
#endif
    }
    // catch failure caused by the H5File operations
    catch( FileIException &error )
    {
      HDF5PrintError(error);
      cerr<<"Error reading file. Exiting "<<endl;
      Fhdf.close();
#ifdef USEMPI
      MPI_Abort(MPI_COMM_WORLD,8);
#else
      exit(8);
#endif
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException &error )
    {
      HDF5PrintError(error);
      cerr<<"Error in data set might suggest config file has the incorrect HDF naming convention. ";
      cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
      Fhdf.close();
#ifdef USEMPI
      MPI_Abort(MPI_COMM_WORLD,8);
#else
      exit(8);
#endif
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException &error )
    {
      HDF5PrintError(error);
      cerr<<"Error in data space might suggest config file has the incorrect HDF naming convention. ";
      cerr<<"Check HDF_name_convetion or update hdfio.cxx in the source code to read correct format"<<endl;
      Fhdf.close();
#ifdef USEMPI
      MPI_Abort(MPI_COMM_WORLD,8);
#else
      exit(8);
#endif
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException &error )
    {
      HDF5PrintError(error);
      cerr<<"Error in data type might suggest need to update hdfio.cxx in the source code to read correct format"<<endl;
      Fhdf.close();
#ifdef USEMPI
      MPI_Abort(MPI_COMM_WORLD,8);
#else
      exit(8);
#endif
    }
    // catch failure caused by missing attribute
    catch( invalid_argument error )
    {
      if(opt.ihdfnameconvention == HDFSWIFTEAGLENAMES) {
        cerr<<"Reading SWIFT EAGLE HDF5 file: "<<error.what()<<endl;
#ifdef USEMPI
        MPI_Abort(MPI_COMM_WORLD, 8);
#else
        exit(8);
#endif
      }
    }
    Fhdf.close();

    for(j=0, nbodies=0; j<nusetypes; j++) {
        k=usetypes[j];
        nbodies+=hdf_header_info.npartTotal[k];
        nbodies+=((long long)(hdf_header_info.npartTotalHW[k]) << 32);
    }
    return nbodies;

}//@}



/// \name Get the number of hdf files per snapshot
//@{
inline Int_t HDF_get_nfiles(char *fname, int ptype)
{
    char buf[2000],buf1[2000],buf2[2000];
    sprintf(buf1,"%s.0.hdf5",fname);
    sprintf(buf2,"%s.hdf5",fname);
    if (FileExists(buf1)) sprintf(buf,"%s",buf1);
    else if (FileExists(buf2)) sprintf(buf,"%s",buf2);
    else {
        printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n", buf1, buf2);
        exit(9);
    }

    H5File Fhdf;
    HDF_Group_Names hdf_gnames;
    //to store the groups, data sets and their associated data spaces
    Attribute headerattribs;
    HDF_Header hdf_header_info;
    //buffers to load data
    int intbuff;
    long long longbuff;
    int ireaderror=0;
    Int_t nfiles = 0;
    IntType inttype;

    //Try block to detect exceptions raised by any of the calls inside it
    try
    {
        //turn off the auto-printing when failure occurs so that we can
        //handle the errors appropriately
        Exception::dontPrint();

        //Open the specified file and the specified dataset in the file.
        Fhdf.openFile(buf, H5F_ACC_RDONLY);
        //get header group

        headerattribs = get_attribute(Fhdf, hdf_header_info.names[hdf_header_info.INumFiles]);
        inttype = headerattribs.getIntType();
        if (inttype.getSize() == sizeof(int))
        {
          headerattribs.read(PredType::NATIVE_INT,&intbuff);
          hdf_header_info.num_files = intbuff;
        }
        if (inttype.getSize() == sizeof(long long))
        {
          headerattribs.read(PredType::NATIVE_LONG,&longbuff);
          hdf_header_info.num_files = longbuff;
        }
    }
    catch(GroupIException &error)
    {
        HDF5PrintError(error);
    }
    // catch failure caused by the H5File operations
    catch( FileIException &error )
    {
        HDF5PrintError(error);
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException &error )
    {
        HDF5PrintError(error);
        ireaderror=1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException &error )
    {
        HDF5PrintError(error);
        ireaderror=1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException &error )
    {
        HDF5PrintError(error);
        ireaderror=1;
    }
    Fhdf.close();

    return nfiles = hdf_header_info.num_files;

}
//@}

#endif
