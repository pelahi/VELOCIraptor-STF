/*! \file nchiladaitems.h
 *  \brief this file contains definitions and routines for reading Nchilada used by Changa nbody code
 *
 * Nchilada format is structure in such a fashion that each particle type has its own
 * directory which then contains a file for every property of that particle type
 * ex: outdir/dark/pos is the file that stores the positions of dark matter particles
 *
 * Some of this code is based off of NBodyShop utiliy package, (like xdr_template.h)
 * and ChanGa N-body code, which uses this data format
 */

#ifndef NCHILADA_STRUCTS_H
#define NCHILADA_STRUCTS_H

#ifdef USEXDR

#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "endianutils.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "tipsy_structs.h"

#define NCHILADAMAXDIM 3

///number of particle types
#define NNCHILADATYPE 4
///\name define Nchilada particle types
//@{
#define NCHILADAGASTYPE 0
#define NCHILADADMTYPE 1
#define NCHILADASTARTYPE 2
#define NCHILADABHTYPE 3
//@}

///\name Nchilada data structures
//@{
///header structure that is in every nchilada file
struct nchilada_dump {
    int    magic;
    double time;
    int    iHighWord;
    int    nbodies;
    int    ndim;
    int    code;
};
//an enumeration that associates a code value with a data type
enum NCDataTypeCode {
    int8 = 1,
    uint8,
    int16,
    uint16,
    int32,
    uint32,
    int64,
    uint64,
    float32,
    float64
};


//@}

///\name XDR read routines, assist in reading tipsy-nchilada data
//@{
///read the tipsy header, and if it fails returns 0
int xdr_header(XDR *pxdrs,struct tipsy_dump *ph);

///read the gas particle, and if it fails returns 0
int xdr_gas(XDR *pxdrs,struct tipsy_gas_particle *ph);

///read the dark particle, and if it fails returns 0
int xdr_dark(XDR *pxdrs,struct tipsy_dark_particle *ph);

///read the star particle, and if it fails returns 0
int xdr_star(XDR *pxdrs,struct tipsy_star_particle *ph);

///reads nchilada header, and if fails returns 0
int xdr_NCHeader(XDR *pxdrs,struct nchilada_dump *ph);

///opens all the associated files for nchilada data format
void xdr_NC_Open(XDR *xdrs, enum xdr_op op, FILE **fp, char *achInFile, long lstart);

///returns the data format pointed to by the XDR which can be then
///used by the \ref NCDataTypeCode to set the data type
int xdr_type(XDR *pxdrs, void *data, int type, int num);

//@}

/*!\name XDR_templates
 Provides inline functions xdr_template() which perform
 XDR conversion of a value.  Numerous overloads
 of this function are provided for many common types, including
*/
//@{
inline bool_t xdr_template(XDR* xdrs, unsigned char* val) {
    return xdr_u_char(xdrs, val);
}
inline bool_t xdr_template(XDR* xdrs, char* val) {
    return xdr_char(xdrs, val);
}
inline bool_t xdr_template(XDR* xdrs, unsigned short* val) {
    return xdr_u_short(xdrs, val);
}
inline bool_t xdr_template(XDR* xdrs, short* val) {
    return xdr_short(xdrs, val);
}
inline bool_t xdr_template(XDR* xdrs, unsigned int* val) {
    return xdr_u_int(xdrs, val);
}
inline bool_t xdr_template(XDR* xdrs, int* val) {
    return xdr_int(xdrs, val);
}
inline bool_t xdr_template(XDR* xdrs, unsigned long* val) {
    return xdr_u_long(xdrs, (unsigned long *)val);
}
inline bool_t xdr_template(XDR* xdrs, long * val) {
    return xdr_long(xdrs, (long *)val);
}
inline bool_t xdr_template(XDR* xdrs, float* val) {
    return xdr_float(xdrs, val);
}
inline bool_t xdr_template(XDR* xdrs, double* val) {
    return xdr_double(xdrs, val);
}

inline bool_t xdr_template(XDR* xdrs, struct nchilada_dump* h) {
    return (xdr_template(xdrs, &(h->magic))
            && xdr_template(xdrs, &(h->time))
            && xdr_template(xdrs, &(h->nbodies))
            && xdr_template(xdrs, &(h->ndim))
            && xdr_template(xdrs, reinterpret_cast<enum_t *>(&(h->code))));
}
//@}


/*! Allocate for and read in a field from an XDR stream.  You need to have
 read the header already.  The min/max pair are put at the end of the array.
 */
template <typename T> inline T* readField(XDR* xdrs, const u_int64_t N, const u_int64_t startParticle = 0) {
    T* data = new T[N + 2];
    //put min/max at the end
    if(!xdr_template(xdrs, data + N) || !xdr_template(xdrs, data + N + 1)) {
        delete[] data;
        return 0;
    }
    if(data[N] == data[N + 1]) {
        //if all elements are the same, just copy the value into the array
        for(u_int64_t i = 0; i < N; ++i)
            data[i] = data[N];
    }
    else {
        ///\todo correct the offseting need for very large arrays
        /*
        ///add the offset from the nchilada header
        off_t offset = sizeof(struct nchilada_dump)
            + (startParticle + 2) * TypeHandling::Type2Dimensions<T>::dimensions * mySizeof(TypeHandling::Type2Code<T>::code);
// XXX NASTY kludge to get around 4 byte limit of xdr functions
#define __BILLION 1000000000
        fseek((FILE *)xdrs->x_private, 0, 0);
        while(offset > __BILLION) {
            offset -= __BILLION;
            fseek((FILE *)xdrs->x_private, __BILLION,SEEK_CUR);
        }
        if(fseek((FILE *)xdrs->x_private, offset, SEEK_CUR) != 0) {
            cerr << "readField seek failed: " << offset << endl;
            delete[] data;
            return 0;
        }
#if 0
        if(!xdr_setpos(xdrs, sizeof(struct nchilada_dump) + (startParticle + 2) * TypeHandling::Type2Dimensions<T>::dimensions * mySizeof(TypeHandling::Type2Code<T>::code))) {
            delete[] data;
            return 0;
        }
#endif
        */
        for(u_int64_t i = 0; i < N; ++i) {
            if(!xdr_template(xdrs, data + i)) {
                delete[] data;
                return 0;
            }
        }
    }
    return data;
}
template <typename T> inline T* readField3D(XDR* xdrs, const u_int64_t N, const u_int64_t startParticle = 0) {
    T* data = new T[3*(N + 2)];
    //put min/max at the end
    Int_t ioffset;
    for (int idim=0;idim<3;idim++) {
        ioffset=idim*(N+2);
        if(!xdr_template(xdrs, data + N+ioffset) || !xdr_template(xdrs, data + N + 1+ioffset)) {
            delete[] data;
            return 0;
        }
        if(data[N+ioffset] == data[N + 1+ioffset]) {
            //if all elements are the same, just copy the value into the array
            for(u_int64_t i = 0; i < N; ++i)
                data[i+ioffset] = data[N+ioffset];
        }
        else {
            /*
            off_t offset = FieldHeader::sizeBytes
                + (startParticle + 2) * TypeHandling::Type2Dimensions<T>::dimensions * mySizeof(TypeHandling::Type2Code<T>::code);
// XXX NASTY kludge to get around 4 byte limit of xdr functions
#define __BILLION 1000000000
            fseek((FILE *)xdrs->x_private, 0, 0);
            while(offset > __BILLION) {
                offset -= __BILLION;
                fseek((FILE *)xdrs->x_private, __BILLION,SEEK_CUR);
            }
            if(fseek((FILE *)xdrs->x_private, offset, SEEK_CUR) != 0) {
                cerr << "readField seek failed: " << offset << endl;
                delete[] data;
                return 0;
            }
#if 0
            if(!xdr_setpos(xdrs, FieldHeader::sizeBytes + (startParticle + 2) * TypeHandling::Type2Dimensions<T>::dimensions * mySizeof(TypeHandling::Type2Code<T>::code))) {
                delete[] data;
                return 0;
            }
#endif
            */
            for(u_int64_t i = 0; i < N; ++i) {
                if(!xdr_template(xdrs, data + i+ioffset)) {
                    delete[] data;
                    return 0;
                }
            }
        }
    }
    return data;
}

template <typename T>
inline void* readField(XDR* xdrs, const unsigned int dimensions, const u_int64_t N, const u_int64_t startParticle = 0) {
    if(dimensions == 3) {
        return readField3D<T>(xdrs, N, startParticle+N);
    }
    else return readField<T>(xdrs, N, startParticle);
}

/** Given the type code in the header, reads in the correct type of data.
 */
inline void* readField(const struct nchilada_dump& fh, XDR* xdrs, u_int64_t numParticles = 0, const u_int64_t startParticle = 0) {
    if(fh.ndim != 1 && fh.ndim != 3) return 0;
    if(numParticles == 0) numParticles = fh.nbodies;
    switch(fh.code) {
        case int8:
            return readField<char>(xdrs, fh.ndim, numParticles, startParticle);
        case uint8:
            return readField<unsigned char>(xdrs, fh.ndim, numParticles, startParticle);
        case int16:
            return readField<short>(xdrs, fh.ndim, numParticles, startParticle);
        case uint16:
            return readField<unsigned short>(xdrs, fh.ndim, numParticles, startParticle);
        case int32:
            return readField<int>(xdrs, fh.ndim, numParticles, startParticle);
        case uint32:
            return readField<unsigned int>(xdrs, fh.ndim, numParticles, startParticle);
        case int64:
            return readField<int64_t>(xdrs, fh.ndim, numParticles, startParticle);
        case uint64:
            return readField<u_int64_t>(xdrs, fh.ndim, numParticles, startParticle);
        case float32:
            return readField<float>(xdrs, fh.ndim, numParticles, startParticle);
        case float64:
            return readField<double>(xdrs, fh.ndim, numParticles, startParticle);
        default:
            return 0;
    }
}

///read a particular data field from a nchilada file
void *readFieldData(FILE *&infile, nchilada_dump &fh, unsigned int dim, u_int64_t numParticles, u_int64_t startParticle);
void *readFieldData(const string filename, nchilada_dump &fh, unsigned int dim, u_int64_t numParticles, u_int64_t startParticle);

/*!\name XDRException handler class
 *
 */
//@{
/*
class XDRException : public std::exception {
  public:
    ~XDRException() throw() { };
    XDRException();
    XDRException(const std::string & desc);
    virtual std::string getText() const throw();
    virtual const char* what() const throw();
  private:
    std::string d;
  };

class XDRReadError : public XDRException {
  public:
    ~XDRReadError() throw() { };
    XDRReadError(std::string s, int64_t w);
    std::string getText() const throw();
 private:
    std::string oper;
    int64_t where;
  };

inline std::ostream & operator <<(std::ostream &os, XDRException &e) {
   os << e.getText();
   return os;
}

inline XDRException::XDRException() : d("") {
}

inline XDRException::XDRException(const std::string & desc)  : d(desc) {
}

inline std::string XDRException::getText() const throw() {
  if(d=="")
    return "Unknown XDR exception";
  else
    return d;
}

inline const char* XDRException::what() const throw() {
  return getText().c_str();
}

inline XDRReadError::XDRReadError(std::string s, int64_t w) : oper(s),
     where(w) {
}

inline std::string XDRReadError::getText() const throw() {
    char sWhere[128];
    sprintf(sWhere, "%ld", (long) where) ;
  return oper + " Error at " + sWhere;
}
*/
//@}

/*! \name Nchilada FieldHeader class
 * The header present in every file of attribute data (a field).
 */
//@{
/*
class FieldHeader {
public:

    static const unsigned int sizeBytes = 28;

    static const int MagicNumber = 1062053;

    int magic;
    double time;
    u_int64_t numParticles;
    unsigned int dimensions; //1 for scalar, 3 for vector
    TypeHandling::DataTypeCode code;

    FieldHeader() : magic(MagicNumber) { }

    FieldHeader(const TypeHandling::TypedArray& arr) : magic(MagicNumber), numParticles(arr.length), dimensions(arr.ndim), code(arr.code) { }

    /// Output operator, used for formatted display
    friend std::ostream& operator<< (std::ostream& os, const FieldHeader& h) {
        return os << "Time: " << h.time
                << "\nTotal number of particles: " << h.numParticles
                << "\nDimensions: " << h.ndim
                << "\nData Type: " << h.code;
    }
};

inline bool_t xdr_template(XDR* xdrs, FieldHeader* h) {
    return (xdr_template(xdrs, &(h->magic))
            && xdr_template(xdrs, &(h->time))
            && xdr_template(xdrs, &(h->numParticles))
            && xdr_template(xdrs, &(h->dimensions))
            && xdr_template(xdrs, reinterpret_cast<enum_t *>(&(h->code))));
}
*/
//@}


///This structures stores the strings defining the names of the data
struct Nchilada_Part_Names {
    //define the strings associated with the types of structures contained in the nchilada file.
    string Header_name;
    string GASpart_name;
    string DMpart_name;
    string STARpart_name;
    string BHpart_name;
    string part_names[NNCHILADATYPE];
    string names[NNCHILADATYPE+1];

    ///constructor
    Nchilada_Part_Names(){
        Header_name=string("descrption.xml");
        GASpart_name=string("gas/");
        DMpart_name=string("dark/");
        STARpart_name=string("star/");
        BHpart_name=string("star/");
        part_names[0]=GASpart_name;
        part_names[1]=DMpart_name;
        part_names[2]=STARpart_name;
        part_names[3]=BHpart_name;

        names[0]=Header_name;
        names[1]=GASpart_name;
        names[2]=DMpart_name;
        names[3]=STARpart_name;
        names[4]=BHpart_name;
    }
};

/// \name Get the number of particles in the nchilada files
//@{
Int_t ncGetCount(string filename);
Int_t Nchilada_get_nbodies(char *fname, int ptype, Options &opt);
//@}
#endif

#endif
