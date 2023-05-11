
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDEigen.h"
#include "TMatrixTUtils.h"
#include "TVectorF.h"

#include "boost/multi_array.hpp"

#include <array>
#include <iostream>
#include <math.h>


//Request doubles
#define TMatrixF TMatrixD
#define TMatrixFSym TMatrixDSym
#define TVectorF TVectorD
#define TMatrixFColumn TMatrixDColumn

typedef std::array<std::array<double,4>,4> Mat4x4;
typedef std::array<double,3> Vec3;

class HydroTmunu
{
private:
public:
    double eps;
    double Tr;
    Vec3 u;

    double pitautau;
    double pitaux;
    double pitauy;
    double pitaueta;
    double pixx;
    double pixy;
    double pixeta;
    double piyy;
    double piyeta;
    double pietaeta;
    HydroTmunu();
    HydroTmunu(double epsilon, double trace, Vec3 u, std::array<double, 10> pi);
    ~HydroTmunu();
};

HydroTmunu::HydroTmunu() :
eps(.00f),Tr(.0f),u({.0f,.0f,.0f}),
pitautau(.0f), pitaux(.0f), pitauy(.0f),
pitaueta(.0f), pixx(.0f), pixy(.0f),
pixeta(.0f), piyy(.0f), piyeta(.0f),
pietaeta(.0f) {}

 HydroTmunu::HydroTmunu(double epsilon, double trace, Vec3 u, std::array<double, 10> pi)
{
    this->eps = epsilon;
    this->Tr = trace;
    this->u = u;
    this->pitautau = pi[0];
    this->pitaux = pi[1];
    this->pitauy = pi[2];
    this->pitaueta = pi[3];
    this->pixx = pi[4];
    this->pixy = pi[5];
    this->pixeta = pi[6];
    this->piyy = pi[7];
    this->piyeta = pi[8];
    this->pietaeta = pi[9];
}

HydroTmunu::~HydroTmunu()
{
}


class Hydrodynamizer
{
private:
    boost::multi_array<Mat4x4,3> Tmunu;
    int nx, ny, neta;
    TMatrixFSym gmunu;


public:
    boost::multi_array<HydroTmunu,3> TmunuOut;
    Hydrodynamizer(boost::multi_array<Mat4x4,3> Tmunu, double tau,
                   int nx, int ny, int neta);
    ~Hydrodynamizer();

    void diagonalize();
};

Hydrodynamizer::Hydrodynamizer(boost::multi_array<Mat4x4,3> Tmunu,
                               double tau,
                               int nx, int ny, int neta) :
Tmunu(Tmunu),
nx(nx),
ny(ny),
neta(neta),
gmunu(TMatrixFSym(4)),
TmunuOut(boost::extents[1][1][1])
{
    //Define metric
    std::array<double,16> gmunu_data;
    for(int mu=0; mu<4; ++mu)
    for(int nu=0; nu<4; ++nu)
        gmunu_data[mu*4+nu] = (mu == nu ? (mu > 0 ? -1 : 1) : 0);

    gmunu_data[15] = (double) -pow(tau,2);

    gmunu.SetMatrixArray(gmunu_data.data());

    //Allocate array of final results
    TmunuOut.resize(boost::extents[nx][ny][neta]);
}

///\brief diagonalize an energy momentum tensor, to get energy density
/// and velocity
///\param T A 4x4 matrix
void Hydrodynamizer::diagonalize(){

    TMatrixFSym CellEMT(4);

    int nx = Tmunu.shape()[0];
    int ny = Tmunu.shape()[1];
    int neta = Tmunu.shape()[2];

    TVectorD zero_vector(4);
    for(int mu=0;mu<4;++mu) zero_vector[mu]=.0;

    int failed_counter = 0;

    std::array<double,16> data;
    TMatrixF eigenvectors(4,4);
    for(int ix=0; ix<nx;++ix)
    for(int iy=0; iy<ny;++iy)
    for(int ieta=0; ieta<neta;++ieta){
        bool is_zero = true;
        for (int mu=0;mu<4;++mu)
        for (int nu=0;nu<4;++nu){
            data[nu*4+mu] = Tmunu[ix][iy][ieta][mu][nu];
            if (fabs(data[nu*4+mu]) > 1.E-10) is_zero = false;
        }

        if( is_zero ){
            TmunuOut[ix][iy][ieta].eps = 0;
            TmunuOut[ix][iy][ieta].u[0] = 0;
            TmunuOut[ix][iy][ieta].u[1] = 0;
            TmunuOut[ix][iy][ieta].u[2] = 0;

            TmunuOut[ix][iy][ieta].Tr = 0;

            TmunuOut[ix][iy][ieta].pitautau = 0;
            TmunuOut[ix][iy][ieta].pitaux = 0;
            TmunuOut[ix][iy][ieta].pitauy = 0;
            TmunuOut[ix][iy][ieta].pitaueta = 0;
            TmunuOut[ix][iy][ieta].pixx = 0;
            TmunuOut[ix][iy][ieta].pixy = 0;
            TmunuOut[ix][iy][ieta].pixeta = 0;
            TmunuOut[ix][iy][ieta].piyy = 0;
            TmunuOut[ix][iy][ieta].piyeta = 0;
            TmunuOut[ix][iy][ieta].pietaeta = 0;

            continue;
        }



        //Allocate data
        CellEMT.SetMatrixArray(data.data());
        //Multiply by the metric
        TMatrixF Tmixed(4,4);
        Tmixed.Mult(CellEMT, gmunu);

        TMatrixDEigen EigenSolver((TMatrixD) Tmixed);
        double eigenval = EigenSolver.GetEigenValuesRe()[0];
        double eigenval_im = EigenSolver.GetEigenValuesIm()[0];

        //Check if matrix eigenvals are real
        if ( fabs(eigenval_im) > 1.E-3){
            failed_counter++;
            std::cout << "[WARNING]: Complex energy density found. Discarding complex part";
        }

        eigenvectors = EigenSolver.GetEigenVectors();
        //TMatrixDRow u = TMatrixDRow(eigenvectors,0);
        TMatrixFColumn u = TMatrixFColumn(eigenvectors,0);


        if (u[0] < 0) u *= -1.;

        if (eigenval < 0) std::cerr<< "Error: Negative energy density detected" << std::endl;

        //Normalizes the vector
        double umu_umu = pow(u[0],2) - pow(u[1],2) - pow(u[2],2) + gmunu(3,3)*pow(u[3],2);
        if (umu_umu > 0){
            double alpha = 1./sqrt(umu_umu);
            u *= alpha;
        } else {
            failed_counter++;
            std::cout << "u^2 < 0 found." << std::endl;
            u[0] = 1;
            u[1] = 0;
            u[2] = 0;
            u[3] = 0;
        }

        if (u[0] < 0) std::cout << "Normalization failed. Negative u0 found" << std::endl;

        if( fabs(1.-sqrt(pow(u[0],2) - pow(u[1],2) - pow(u[2],2) + gmunu(3,3)*pow(u[3],2))) > 1.E-3){
            std::cout << "Normalization failed. Not normalized vector found" << std::endl;
            std::cout << u[0] << ", " << u[1] << ", " << u[2] << ", " << u[3] << std::endl;
        }

        TmunuOut[ix][iy][ieta].eps = eigenval;
        TmunuOut[ix][iy][ieta].u[0] = u[1];
        TmunuOut[ix][iy][ieta].u[1] = u[2];
        TmunuOut[ix][iy][ieta].u[2] = u[3];

        double Tr = Tmixed(0,0)+Tmixed(1,1)+Tmixed(2,2)+Tmixed(3,3);
        TmunuOut[ix][iy][ieta].Tr = Tr;

        TmunuOut[ix][iy][ieta].pitautau
            = CellEMT(0,0) + (u[0]*u[0]*(Tr-4*eigenval) + gmunu(0,0)*(eigenval-Tr))/3.;
        TmunuOut[ix][iy][ieta].pitaux
            = CellEMT(0,1) + (u[0]*u[1]*(Tr-4*eigenval))/3.;
        TmunuOut[ix][iy][ieta].pitauy
            = CellEMT(0,2) + (u[0]*u[2]*(Tr-4*eigenval))/3. ;
        TmunuOut[ix][iy][ieta].pitaueta
            = CellEMT(0,3) + (u[0]*u[3]*(Tr-4*eigenval))/3.;
        TmunuOut[ix][iy][ieta].pixx
            = CellEMT(1,1) + (u[1]*u[1]*(Tr-4*eigenval) + gmunu(1,1)*(eigenval-Tr))/3.;
        TmunuOut[ix][iy][ieta].pixy
            = CellEMT(1,2) + (u[1]*u[2]*(Tr-4*eigenval))/3.;
        TmunuOut[ix][iy][ieta].pixeta
            = CellEMT(1,3) + (u[1]*u[3]*(Tr-4*eigenval))/3.;
        TmunuOut[ix][iy][ieta].piyy
            = CellEMT(2,2) + (u[2]*u[2]*(Tr-4*eigenval) + gmunu(2,2)*(eigenval-Tr))/3.;
        TmunuOut[ix][iy][ieta].piyeta
            = CellEMT(2,3) + (u[2]*u[3]*(Tr-4*eigenval))/3.;
        TmunuOut[ix][iy][ieta].pietaeta
            = CellEMT(3,3) + (u[3]*u[3]*(Tr-4*eigenval) + (eigenval-Tr)/gmunu(3,3))/3.; //We are computing contravariant components


        //Check tensor:
        double shear_Tr =   TmunuOut[ix][iy][ieta].pitautau - TmunuOut[ix][iy][ieta].pixx
                          - TmunuOut[ix][iy][ieta].piyy + TmunuOut[ix][iy][ieta].pietaeta*gmunu(3,3);
        double umu_pimutau = u[0]*TmunuOut[ix][iy][ieta].pitautau - u[1]*TmunuOut[ix][iy][ieta].pitaux
                          - u[2]*TmunuOut[ix][iy][ieta].pitauy + TmunuOut[ix][iy][ieta].pitaueta*u[3]*gmunu(3,3);
        double umu_pimux = u[0]*TmunuOut[ix][iy][ieta].pitaux - u[1]*TmunuOut[ix][iy][ieta].pixx
                          - u[2]*TmunuOut[ix][iy][ieta].pixy + TmunuOut[ix][iy][ieta].pixeta*u[3]*gmunu(3,3);
        double umu_pimuy = u[0]*TmunuOut[ix][iy][ieta].pitauy - u[1]*TmunuOut[ix][iy][ieta].pixy
                          - u[2]*TmunuOut[ix][iy][ieta].piyy + TmunuOut[ix][iy][ieta].piyeta*u[3]*gmunu(3,3);
        double umu_pimueta = u[0]*TmunuOut[ix][iy][ieta].pitaueta - u[1]*TmunuOut[ix][iy][ieta].pixeta
                          - u[2]*TmunuOut[ix][iy][ieta].piyeta + TmunuOut[ix][iy][ieta].pietaeta*u[3]*gmunu(3,3);
        if(abs(shear_Tr) > 1.E-7) std::cout << "Trace of shear tensor: "<<shear_Tr<<std::endl;
        if(abs(umu_pimutau) > 1.E-7) std::cout << "u not orthogonal to pi^{mu tau}. Value is: "<<umu_pimutau<<std::endl;
        if(abs(umu_pimux) > 1.E-7) std::cout << "u not orthogonal to pi^{mu x}: "<<umu_pimux<<std::endl;
        if(abs(umu_pimuy) > 1.E-7) std::cout << "u not orthogonal to pi^{mu y}: "<<umu_pimuy<<std::endl;
        if(abs(umu_pimueta) > 1.E-7) std::cout << "u not orthogonal to pi^{mu eta}: "<<umu_pimueta<<std::endl;


        //if(abs)


    }
    if (failed_counter > 0)
        std::cout << "[WARNING]: Decomposition failed in "<< failed_counter << " cells out of "<< nx*ny*neta << std::endl;
}

Hydrodynamizer::~Hydrodynamizer()
{
}
