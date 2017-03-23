/*

compconv.f

C++ wrapper for the comptonise XSPEC convoluation model.
Solves the Kompaneets equation using the seed photons from the
input spectrum.

Based on nthcomp; the thermally comptonized continuum model of 
Zdziarski, Johnson & Magdziarz 1996, MNRAS, 283, 193
as extended by Zycki, Done & Smith 1999, MNRAS 309, 561

Calculation subroutines (Fortran) are in docompconf.f

v2.0, July 2016, D.R. Wilkins

*/

#include <xsTypes.h>
#include <functionMap.h>
#include <stlToCArrays.h>
#include <XSUtil/Utils/XSutility.h>
#include <cfortran.h>

// Functions called FROM here to Fortran:

PROTOCCALLSFSUB6(DOCOMPCONV,docompconv,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)
#define DOCOMPCONV(ear,ne,param,ifl,photar,prim) \
           CCALLSFSUB6(DOCOMPCONV,docompconv,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV,\
             ear,ne,param,ifl,photar,prim)


extern "C" void compconv(const RealArray& energyArray, const RealArray& params,
        int spectrumNumber, RealArray& flux, RealArray& fluxErr, 
        const string& initString)
{
   // Actual work is done in donthcomp.f

   using namespace XSutility;

   int ne = static_cast<int>(energyArray.size()) - 1;

   float *ear=0, *pars=0, *photar=0, *photer=0;
   XSFunctions::stlToFloatArrays<float>(energyArray, params, flux, fluxErr,
           ear, pars, photar, photer);
   auto_array_ptr<float> apEar(ear);
   auto_array_ptr<float> apPars(pars);
   auto_array_ptr<float> apPhotar(photar);
   auto_array_ptr<float> apPhoter(photer);

   DOCOMPCONV(ear, ne, pars, spectrumNumber, photar, photer);

   XSFunctions::floatFluxToStl<float>(photar, photer, ne, false, flux, fluxErr);       

}

