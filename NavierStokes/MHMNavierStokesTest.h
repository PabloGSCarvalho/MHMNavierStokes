/*
 *  MHMNavierStokesTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__MHMNavierStokesTest__
#define __PZ__MHMNavierStokesTest__

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZMHMNavierStokesMaterial.h"
#include "TPZMultiphysicsCompMesh.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzinterpolationspace.h"
#include "pztrnsform.h"
#include "TPZSimulationData.h"
#include "TPZMHMNavierStokesMeshControl.h"
#include "TPZAnalyticSolution.h"
#include "TPZNSAnalysis.h"

using namespace std;
using namespace pzshape;

class MHMNavierStokesTest{
private:
    
    int fdim; //Dimensão do problema
    int fmatID; //Materia do elemento volumétrico
        
    //Materiais das condições de contorno
    int fmatBCbott;
    int fmatBCtop;
    int fmatBCleft;
    int fmatBCright;

    int fmatBCleft_2;
    int fmatBCright_2;

    int fmatBChole;

    int fmatBCtop_z;
    int fmatBCbott_z;//normal negativa
    
    //Material do elemento de interface
    int fmatLambda; //Multiplier material
    int fmatLambdaBC;
    
    int fmatLambdaBC_bott;
    int fmatLambdaBC_top;
    int fmatLambdaBC_left;
    int fmatLambdaBC_right;

    int fmatLambdaBC_hole;
    
    int fmatLambdaBC_top_z;
    int fmatLambdaBC_bott_z;

    int fmatLambdaBC_left_2;
    int fmatLambdaBC_right_2;
    
    int fmatInterfaceLeft;
    int fmatInterfaceRight;
    int fmatWrap;
    
    //Materiais das condições de contorno (elementos de interface)

    int fmatIntBCtop_z;
    int fmatIntBCbott_z;

    std::map<int,TPZManVector<REAL,3>> f_HoleCoord; //Dado o indice do elemento 2D, devolve a coord do hole associado

    std::map<int,int> f_ArcCentralNode; //Dado o indice do elemento 2D, devolve o indice do nó central do arco
    
    //Materia de um ponto
    int fmatPoint;
    
    //Condições de contorno do problema
    int fdirichlet_v;
    int fneumann_v;

    int fdirichlet_sigma;
    int fneumann_sigma;

    int f_BJS_condition;

    int fpenetration;
    int fpointtype;
    int fdirichletvar;
    
    int fquadmat1; //Parte inferior do quadrado
    int fquadmat2; //Parte superior do quadrado
    int fquadmat3; //Material de interface
    
    STATE fviscosity;
    STATE fpermeability;
    STATE ftheta;
    
    REAL fphi_r;
    
    bool f_is_hdivFull;
    
    bool f_hdivPlus;
    
    MElementType feltype;

    TStokesAnalytic::MProblemType f_problemtype;

    TStokesAnalytic::EExactSol f_domaintype;

    TStokesAnalytic f_ExactSol;

    TStokesAnalytic f_ExactSol_2;

    bool f_Holemesh = false;
    
    TPZVec<TPZCompMesh * > f_mesh_vector;
    
    static TPZTransform<STATE> f_T;
    
    static TPZTransform<STATE> f_InvT;
    
    bool f_allrefine = false;
    
    TPZGeoMesh *f_mesh0;

    TPZStack<TPZGeoElSide> f_skellNeighs;
    
    bool f_3Dmesh = false;

    TPZSimulationData *f_sim_data;
    
public:

    MHMNavierStokesTest();
    
    ~MHMNavierStokesTest();
  
    MHMNavierStokesTest(const MHMNavierStokesTest &cp)
    {
        DebugStop();
    }

    MHMNavierStokesTest &operator =(const MHMNavierStokesTest &cp)
    {
        DebugStop();
        return *this;
    }


    void Run();
    
    void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix );

    void SolveNonLinearProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix );

    std::ostream &ConfigPrint(std::ostream &out);
    
    /*  Malhas geometricas */
    TPZGeoMesh *CreateGMesh(TPZVec<int> &n_s, TPZVec<REAL> &h_s);

    TPZGeoMesh *CreateGMesh3D(TPZVec<int> &n_s, TPZVec<REAL> &h_s);

    TPZGeoMesh *CreateGMeshCoupling(TPZVec<int> &n_s, TPZVec<REAL> &h_s);

    TPZGeoMesh *CreateGMeshCurve();

    TPZGeoMesh *CreateGMeshRefPattern(TPZVec<int> &n_div, TPZVec<REAL> &h_s);

    TPZAutoPointer<TPZRefPattern> CreateGMeshObstacle(int nrefs, TPZManVector<REAL,6> &FirstCoord, TPZManVector<REAL,6> &h_el);

    TPZGeoMesh *CreateGMeshVugsRefPattern(TPZVec<int> &n_div, TPZVec<REAL> &h_s);

    TPZAutoPointer<TPZRefPattern> CreateGMeshVugs(int nrefs, TPZManVector<REAL,6> &FirstCoord, TPZManVector<REAL,6> &h_el);

    TPZAutoPointer<TPZRefPattern> CreateGMeshQuadRef(int nrefs, TPZManVector<REAL,6> &h_el);

    TPZManVector<REAL,3> ParametricCircle(REAL radius,REAL theta);
    
    static void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);
    
    void TetrahedralMeshCubo(TPZVec<int> &n_s);
    
    void SubdomainRefine(int nDiv, TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices);
    
    void UniformRefine4(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction); //Elemento escolhido pela coordenada -> grande estrutura
    
    void UniformRefine3(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div); //Refinamento padrão para triangulo
    
    void UniformRefine2(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div); //Refinamento padrão para quadrilatero
    
    void UniformRefine(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction); //Elemento escolhido pela coordenada
    //   TPZGeoMesh *GMeshDeformed(int dim, bool ftriang, int ndiv);
    
    void InsertLowerDimMaterial(TPZGeoMesh *gmesh);
    
    void ChangeExternalOrderConnects(TPZCompMesh *mesh, int addToOrder);
    /* Malhas computacionais */
    
    TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh, int64_t &index);
    
    void InsertMaterialObjects(TPZMHMeshControl *control);
        
    
    //TPZCompMesh *CMesh_St(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZMultiphysicsCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);

    void SetOriginalMesh(TPZGeoMesh *gmesh){
        f_mesh0 = gmesh;
        ComputeSkelNeighbours(); //papapa ?????
    };

    void SetSimulationData(TPZSimulationData *simdata){
        f_sim_data = simdata;
    };
    
    
    void SetAllRefine(){
        f_allrefine = true;
    };

    void Set3Dmesh(){
        f_3Dmesh = true;
        fdim = 3;
        f_ExactSol.fDimension = 3;
    };
    
    void SetHdivPlus(){
        f_hdivPlus = true;
    };
    
    
    void SetFullHdiv(){
        f_is_hdivFull = true;
    };
    
    void SetElType(MElementType eltype){
        feltype = eltype;
    };

    void SetProblemType(TStokesAnalytic::MProblemType type){
        f_problemtype = type;
        f_ExactSol.fProblemType = type;
    };

    void SetDomainType(TStokesAnalytic::EExactSol type){
        f_domaintype = type;
        f_ExactSol.fExactSol = type;
    };

    // Set transform object and its transformation
    void SetTransform(TPZTransform<STATE> Transf, TPZTransform<STATE> InvTransf){
        f_T = Transf;
        f_InvT = InvTransf;
    }
    
    void GetElIndexCoarseMesh(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindex);
   
    void SetRotationMatrix(REAL theta){
        TPZFMatrix<STATE> TfMatrix(3,3,0.);
        TfMatrix(0,0)= cos(theta);
        TfMatrix(0,1)= -sin(theta);
        TfMatrix(1,0)= sin(theta);
        TfMatrix(1,1)= cos(theta);
        TfMatrix(2,2)= 1.;

        TPZFMatrix<REAL> sum(3,1,0.); //only rotation
        TPZFMatrix<STATE> Inv(0,0,0.);
        f_T.SetMatrix(TfMatrix, sum);
        
        TfMatrix.Inverse(Inv, ENoDecompose);
        f_InvT.SetMatrix(Inv, sum);
        
    }
    
    void SetRotation3DMatrix(REAL rot_x, REAL rot_y, REAL rot_z){
        TPZFMatrix<STATE> RotX(3,3,0.),RotY(3,3,0.),RotZ(3,3,0.),R(3,3,0.),RZY(3,3,0.);
        RotX(0,0)= 1.;
        RotX(1,1)= cos(rot_x);
        RotX(1,2)= -sin(rot_x);
        RotX(2,1)= sin(rot_x);
        RotX(2,2)= cos(rot_x);
        
        RotY(0,0)= cos(rot_y);
        RotY(0,2)= sin(rot_y);
        RotY(1,1)= 1.;
        RotY(2,0)= -sin(rot_y);
        RotY(2,2)= cos(rot_y);
        
        RotZ(0,0)= cos(rot_z);
        RotZ(0,1)= -sin(rot_z);
        RotZ(1,0)= sin(rot_z);
        RotZ(1,1)= cos(rot_z);
        RotZ(2,2)= 1.;
        
        RotZ.Multiply(RotY,RZY);
        RZY.Multiply(RotX,R);
        
        TPZFMatrix<REAL> sum(3,1,0.); //only rotation
        TPZFMatrix<STATE> Inv(0,0,0.);
        f_T.SetMatrix(R, sum);
        
        R.Inverse(Inv, ENoDecompose);
        f_InvT.SetMatrix(Inv, sum);
    }
    
    //solucao exata
    static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //lado direito da equacao
    static void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
    // static void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);
    void AddMultiphysicsInterfaces(TPZMultiphysicsCompMesh &cmesh);
    
    void AddMultiphysicsInterfaces(TPZMultiphysicsCompMesh &cmesh, int matfrom, int mattarget);
    
    void AddMultiphysicsInterfacesLeftNRight(TPZMultiphysicsCompMesh &cmesh, int matfrom);
    
    // Rotate function
    void Rotate(TPZVec<REAL> &co, TPZVec<REAL> &co_r, bool rotate);
    
    // Insere interfaces na malha multifísica
    void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh);
    
    // Agrupa elementos e realiza condensação estática
    void GroupAndCondense(TPZMultiphysicsCompMesh *cmesh_m);
    
    void ComputeSkelNeighbours();

    void InsertArcInterface(TPZAutoPointer<TPZGeoMesh> gmesh); //Insert arcs between Stokes and Darcy materials
    
    bool IsSkellNeighbour(TPZGeoElSide neigh);

};


#endif 

