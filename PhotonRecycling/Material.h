#ifndef MATERIAL_HEADER
#define MATERIAL_HEADER

#include <TROOT.h>
#include <TSpline.h>
#include <iostream>

class Material
{
public:
    Material();
    virtual ~Material();
    virtual void Setn(double *x, double *y, int n);
    virtual void Setkapa(double *x, double *y, int n);
    virtual void Setalpha(double *x, double *y, int n);
    virtual double n(double lambda);    // real part of refractive index
    virtual double kapa(double lambda); // imaginary part of refractive index
    virtual double alpha(double lambda);

private:
    bool doncurve;
    bool dokapacurve;
    bool doalphacurve;
    TSpline3 *ncurve;
    TSpline3 *kapacurve;
    TSpline3 *alphacurve;
};

#endif