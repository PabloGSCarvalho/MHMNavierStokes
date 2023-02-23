/*
 *  NavierStokesTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__NavierStokesTest__
#define __PZ__NavierStokesTest__

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "TPZLinearAnalysis.h"
#include "TPZBndCondT.h"
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
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZLinearAnalysis.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "pzinterpolationspace.h"
#include "pztrnsform.h"
#include "TPZSimulationData.h"
#include "TPZNSAnalysis.h"
#include "TPZAnalyticSolution.h"

using namespace std;
using namespace pzshape;



//enum NSDomainType {EObstacle,EOneCurve,ESinCos};

class NavierStokesTest{
    
private:
    
    int fdim; //Dimensão do problema
    std::set<int> fmatIDs; //Materia do elemento volumétrico
        
    //Materiais das condições de contorno
    int fmatBCbott;
    int fmatBCtop;
    int fmatBCleft;
    int fmatBCright;
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
    
    int fmatWrapBC_bott;
    int fmatWrapBC_top;
    int fmatWrapBC_left;
    int fmatWrapBC_right;
    
    int fmatInterfaceLeft;
    int fmatInterfaceRight;
    int fmatWrap;
    
    //Materiais das condições de contorno (elementos de interface)
    int fmatIntBCbott;
    int fmatIntBCtop;
    int fmatIntBCleft;
    int fmatIntBCright;

    int fmatIntBCtop_z;
    int fmatIntBCbott_z;

    //Materia de um ponto
    int fmatPoint;
    
    //Condições de contorno do problema
    int fdirichlet_v;
    int fneumann_v;
    
    int fdirichlet_sigma;
    int fneumann_sigma;
    
    int fpenetration;
    int fpointtype;
    int fdirichletvar;
    
    int fquadmat1; //Parte inferior do quadrado
    int fquadmat2; //Parte superior do quadrado
    int fquadmat3; //Material de interface
    
    STATE fviscosity;
    STATE fpermeability;
    STATE ftheta;
    
    int fSpaceV;
    
    REAL fphi_r;
    
    bool f_is_hdivFull;
    
    int f_fluxOrder = 1;
    
    int f_tractionOrder = 0;
    
    int f_hdivPlus = 0;
    
    MElementType feltype;
    
    TStokesAnalytic::MProblemType f_problemtype;
    
    TStokesAnalytic::EExactSol f_domaintype;
    
    TStokesAnalytic f_ExactSol;
    
    TPZManVector<TPZCompMesh *, 6> f_mesh_vector;
    
    std::map<int,TPZManVector<REAL,3>> f_HoleCoord; //Dado o indice do elemento 2D, devolve a coord do hole associado

    std::map<int,int> f_ArcCentralNode; //Dado o indice do elemento 2D, devolve o indice do nó central do arco
    
    static TPZTransform<STATE> f_T;
    
    static TPZTransform<STATE> f_InvT;
    
    bool f_allrefine = false;
    
    int f_inter_ref = 0.;
    
    TPZGeoMesh *f_mesh0;
    
    TPZStack<TPZGeoElSide> f_skellNeighs;
    
    bool f_3Dmesh = false;

    bool f_StokesTest = false;
    
    bool f_OseenTest = false;
    
    bool f_CurveTest = false;

    TPZSimulationData *f_sim_data;
    
public:

    NavierStokesTest();
    
    ~NavierStokesTest();
    
    NavierStokesTest(const NavierStokesTest &cp)
    {
        DebugStop();
    }

    NavierStokesTest &operator =(const NavierStokesTest &cp)
    {
        DebugStop();
        return *this;
    }

    void Run(int pOrder, TPZVec<int> &n_s, TPZVec<REAL> &h_s);
    
    /*  Malhas geometricas */
    
    TPZGeoMesh *CreateGMeshCurve();
    
    //TPZGeoMesh *CreateGMeshObstacle();
    
    TPZAutoPointer<TPZRefPattern> CreateGMeshObstacle(TPZManVector<REAL,6> &FirstCoord, TPZManVector<REAL,6> &h_el);
    
    TPZGeoMesh *CreateGMeshCurveBlend();
    
    TPZGeoMesh *CreateGMeshCurveBlendSimple();
    
    TPZGeoMesh *CreateGMeshRefPattern(TPZVec<int> &n_s, TPZVec<REAL> &h_s);
    
    TPZGeoMesh *CreateGMesh(TPZVec<int> &n_s, TPZVec<REAL> &h_s);

    TPZGeoMesh *CreateGMesh3D(TPZVec<int> &n_s, TPZVec<REAL> &h_s);
    
    TPZManVector<REAL,3>  ParametricCircle(REAL radius,REAL theta);
    
    static void GenerateNodes(TPZGeoMesh *gmesh, int64_t nelem);
    
    void TetrahedralMeshCubo(TPZVec<int> &n_s);
    
    void UniformRefine4(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction); //Elemento escolhido pela coordenada -> grande estrutura
    
    void UniformRefine3(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div); //Refinamento padrão para triangulo
    
    void UniformRefine2(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div); //Refinamento padrão para quadrilatero
    
    void UniformRefine(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction); //Elemento escolhido pela coordenada
    //   TPZGeoMesh *GMeshDeformed(int dim, bool ftriang, int ndiv);
    
    void InsertLowerDimMaterial(TPZGeoMesh *gmesh);
    
    void ChangeExternalOrderConnects(TPZCompMesh *mesh, int addToOrder);
    /* Malhas computacionais */
    
    TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh);
    
    TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh);
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh);
    
    TPZCompMesh *CMesh_pM(TPZGeoMesh *gmesh, int pOrder);
    TPZCompMesh *CMesh_gM(TPZGeoMesh *gmesh, int pOrder);

    TPZCompMesh *CMesh_pM_0(TPZGeoMesh *gmesh, int pOrder);
    TPZCompMesh *CMesh_gM_0(TPZGeoMesh *gmesh, int pOrder);
    
    
    //TPZCompMesh *CMesh_St(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZMultiphysicsCompMesh *CMesh_m(TPZGeoMesh *gmesh, int pOrder);

    TPZMultiphysicsCompMesh *CMesh_m_cavity(TPZGeoMesh *gmesh, int pOrder);

    void SetOriginalMesh(TPZGeoMesh *gmesh){
        f_mesh0 = gmesh;
        ComputeSkelNeighbours();
    };

    void SetSimulationData(TPZSimulationData *sim_data){
        f_sim_data = sim_data;
    }
    
    void SetAllRefine(){
        f_allrefine = true;
    }


    void Set3Dmesh(){
        f_3Dmesh = true;
        fdim = 3;
    }
    
    void SetHdivPlus(int order){
        f_hdivPlus = order;
    }
    
    void SetTractionOrder(int order)
    {
        f_tractionOrder = order;
    }
    
    void SetFluxOrder(int order)
    {
        f_fluxOrder = order;
    }
    
    
    void SetFullHdiv(){
        f_is_hdivFull = true;
    };

    void SetProblemType(TStokesAnalytic::MProblemType type){
        f_problemtype = type;
        f_ExactSol.fProblemType = type;
    };

    void SetDomainType(TStokesAnalytic::EExactSol type){

        f_ExactSol.fExactSol = type;
    };
    
    void SetElType(MElementType eltype){
        feltype = eltype;
    };

    // Set transform object and its transformation
    void SetTransform(TPZTransform<STATE> Transf, TPZTransform<STATE> InvTransf){
        f_T = Transf;
        f_InvT = InvTransf;
    }
    
    void SetInternRef(int inter_ref){
        f_inter_ref = inter_ref;
    }
    
//    void SetTransfMatrix(TPZFMatrix<STATE> TfMatrix){
//        f_Tmatrix = TfMatrix;
//        TPZFMatrix<STATE> Inv(0,0,0.);
//        f_Tmatrix.Inverse(Inv, ENoDecompose);
//        f_InvTmatrix = Inv;
//        TPZFMatrix<REAL> sum(3,3,0.); //only rotation
//        if(f_T.Mult().Rows()>0 && f_InvT.Mult().Rows()>0){
//            f_T.SetMatrix(f_Tmatrix, sum);
//            f_InvT.SetMatrix(f_InvTmatrix, sum);
//        }
//    }

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

    //solucao exata
    static void Sol_exact_Curve(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //solucao exata
    static void Sol_exact_Oseen(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //lado direito da equacao
    static void F_source_Oseen(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
    
    //solucao exata
    static void Sol_exact_Stokes(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //lado direito da equacao
    static void F_source_Stokes(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
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
    
    bool IsSkellNeighbour(TPZGeoElSide neigh);

    void VerifyEquilibrium(TPZCompMesh *cmesh);
    
};


#endif 
