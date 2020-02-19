

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
#include "StokesTest.h"
#include "BrinkmanTest.h"
#include "HybridBrinkmanTest.h"
#include "NavierStokesTest.h"
#include "MHMStokesTest.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "pzgengrid.h"

#include "TPZStokesMaterial.h"
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
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "pztrnsform.h"

#include "fad.h"
#include "fadType.h"

//HDivPiola = 1;
const int SpaceHDiv = 1; //Velocidade em subespaço de H(div)
const int SpaceContinuous = 2; //Velocidade em subespaço de [H1]ˆ2
const int SpaceDiscontinuous = 3; //Velociadade em subespaço de H(Ph) - Ph: partição
const REAL Pi=M_PI;

//Verificação dos modelos:

const REAL visco=1., permeability=1., theta=-1.; //Coeficientes: viscosidade, permeabilidade, fator simetria

bool StokesDomain = false , BrinkmanDomain = false;

bool HybridBrinkmanDomain = false, MHMStokesDomain = false, NavierStokesDomain = true;

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e16;
//    gRefDBase.InitializeAllUniformRefPatterns();
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    REAL hx=2.,hy=2.; //Dimensões em x e y do domínio
    //double hx=Pi,hy=2.;
    int h_level = 0;
    int nx=h_level+1 ,ny=h_level+1; //Número de nos em x  y
    int pOrder = 0; //Ordem polinomial de aproximação
    
    TPZVec<REAL> h_s(3,0);
    h_s[0]=2.,h_s[1]=2.,h_s[2]=2.; //Dimensões em x e y do domínio

    if (NavierStokesDomain){
        
        int pOrder = 2;
        
        for (int it=0; it<=0; it++) {
            //h_level = pow(2., 2+it);
            h_level = 4;
            
            TPZVec<int> n_s(3,0.);
            n_s[0]=h_level ,n_s[1]=h_level;
            
            n_s[2]=h_level; //Obs!!
            
            REAL visc = 1.0; //->Darcy
            
            NavierStokesTest  * Test2 = new NavierStokesTest();
            //Test2->Set3Dmesh();
            //Test2->SetProblemType(NSObstacle);
            
            
            //Test2->SetElType(ETriangle);
            Test2->SetInternRef(0);
            Test2->SetStokesTest();
            //Test2->SetHdivPlus();

            TPZTransform<STATE> Transf(3,3), InvTransf(3,3);
            Test2->SetTransform(Transf, InvTransf);

            REAL rot_x = 5.;
            REAL rot_z = 44.;
            REAL rot_y = -85.;
            rot_z = rot_z*Pi/180.;
            rot_y = rot_y*Pi/180.;
            rot_z = rot_z*Pi/180.;
            
            //Test2->SetRotation3DMatrix(rot_x,rot_y,rot_z);
            //Test2->SetAllRefine();
            Test2->Run(SpaceHDiv, pOrder, n_s, h_s,visc);
            
        }
        
    }
    else if (MHMStokesDomain)
    {
        
        for (int it=0; it<=0; it++) {
            //h_level = pow(2., 1+it);
            h_level = 4;
            
            TPZVec<int> n_s(3,0.);
            n_s[0]=h_level,n_s[1]=h_level;
            n_s[2]=h_level; //Obs!!
            
            MHMStokesTest  * Test2 = new MHMStokesTest();
            //Test2->Set3Dmesh();
            //Test2->SetElType(ECube);
            //Test2->SetHdivPlus();
            
            TPZTransform<STATE> Transf(3,3), InvTransf(3,3);
            Test2->SetTransform(Transf, InvTransf);
            
            REAL rot_x = 5.;
            REAL rot_z = 44.;
            REAL rot_y = -85.;
            rot_z = rot_z*Pi/180.;
            rot_y = rot_y*Pi/180.;
            rot_z = rot_z*Pi/180.;
            
            //Test2->SetRotation3DMatrix(rot_x,rot_y,rot_z);
            TPZSimulationData simdata;
            simdata.SetInternalOrder(2);
            simdata.SetSkeletonOrder(1);
            simdata.SetCoarseDivisions(n_s);
            simdata.SetDomainSize(h_s);
            simdata.SetNInterRefs(0);
            simdata.SetViscosity(1.);
            simdata.SetNthreads(0);
            //simdata.SetShapeTest(); // Test for shape functions
            
            Test2->SetSimulationData(simdata);
            Test2->Run();
            
        }
        
    }
    else if (HybridBrinkmanDomain){
        
        int pOrder = 1;
        
        for (int it=0; it<=0; it++) {
            //h_level = pow(2., 2+it);
            h_level = 2;
            
            TPZVec<int> n_s(3,0.);
            n_s[0]=h_level ,n_s[1]=h_level;
            
            n_s[2]=h_level; //Obs!!
            
            REAL visc = 1.0; //->Darcy
            
            HybridBrinkmanTest  * Test2 = new HybridBrinkmanTest();
            //Test2->Set3Dmesh();
            Test2->SetProblemType(ObstacleP);
            
            
            //Test2->SetElType(ETriangle);
            Test2->SetInternRef(1);
            //Test2->SetHdivPlus();
            
            TPZTransform<STATE> Transf(3,3), InvTransf(3,3);
            Test2->SetTransform(Transf, InvTransf);
            
            REAL rot_x = 5.;
            REAL rot_z = 44.;
            REAL rot_y = -85.;
            rot_z = rot_z*Pi/180.;
            rot_y = rot_y*Pi/180.;
            rot_z = rot_z*Pi/180.;
            
            //Test2->SetRotation3DMatrix(rot_x,rot_y,rot_z);
            //Test2->SetAllRefine();
            Test2->Run(SpaceHDiv, pOrder, n_s, h_s,visc);
            
        }
        
    }
    
    return 0;
}

