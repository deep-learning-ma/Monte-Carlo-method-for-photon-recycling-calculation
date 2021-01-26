#include "Material.h"

using namespace std;

Material::Material()
{
    doncurve = false;
    dokapacurve = false;
    doalphacurve = false;
}

Material::~Material()
{
    if (doncurve)
        delete ncurve;
    if (dokapacurve)
        delete kapacurve;
    if (doalphacurve)
        delete alphacurve;
}

void Material::Setn(double *x, double *y, int n)
{
    ncurve = new TSpline3("ncurve", x, y, n);
    doncurve = true;
}

void Material::Setkapa(double *x, double *y, int n)
{
    kapacurve = new TSpline3("kapacurve", x, y, n);
    dokapacurve = true;
}

void Material::Setalpha(double *x, double *y, int n)
{
    alphacurve = new TSpline3("alphacurve", x, y, n);
    doalphacurve = true;
}

double Material::n(double lambda)
{
    return ncurve->Eval(lambda);
}

double Material::kapa(double lambda)
{
    return kapacurve->Eval(lambda);
}

double Material::alpha(double lambda)
{
    return alphacurve->Eval(lambda);
}