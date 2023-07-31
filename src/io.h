#ifndef SRC_IO_H
#define SRC_IO_H

///Math code
#include <NBodyMath.h>

/// Return true if if given file exists
bool FileExists(const char *fname);

/// @brief Wrap input for peroid 
/// @param p Period
/// @param x x coord
/// @param y y coord
/// @param z z coord
template<class T> void PeriodWrapInput(Double_t p, T &x, T &y, T &z){
    Double_t pfac;
    pfac = std::floor(x/p);
    x -=  pfac*p;
    pfac = std::floor(y/p);
    y -=  pfac*p;
    pfac = std::floor(z/p);
    z -=  pfac*p;
};
/// @brief Wrap input for peroid
/// @param p 
/// @param pos 
template<class T> void PeriodWrapInput(Double_t p, T *pos) {
    Double_t pfac;
    for (auto i=0;i<3;i++) {
        pfac = std::floor(pos[i]/p);
        pos[i] -=  pfac*p;
    }
}

#endif /* SRC_IO_H */