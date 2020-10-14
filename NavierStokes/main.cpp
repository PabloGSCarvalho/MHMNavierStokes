

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
#include "NavierStokesTest.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "TPZGenGrid2D.h"

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
#include "TPZAnalyticSolution.h"
#include "MHMNavierStokesTest.h"

#include "fad.h"
#include "fadType.h"

//HDivPiola = 1;
const int SpaceHDiv = 1; //Velocidade em subespaço de H(div)
const int SpaceContinuous = 2; //Velocidade em subespaço de [H1]ˆ2
const int SpaceDiscontinuous = 3; //Velociadade em subespaço de H(Ph) - Ph: partição
const REAL Pi=M_PI;

//Verificação dos modelos:

const REAL visco=1., permeability=1., theta=-1.; //Coeficientes: viscosidade, permeabilidade, fator simetria
//bool MHMProblem = true; //True for MHM problem, False for hybrid formulation problem
//bool HybridProblem = false;
enum Simulation_case {MHMProblem, HybridProblem, Coupling, Vugs2D};

int main(int argc, char *argv[])
{
    
    TPZMaterial::gBigNumber = 1.e10;
//    gRefDBase.InitializeAllUniformRefPatterns();
    Simulation_case sim_case = Vugs2D;
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    //Dados do problema:
    
    REAL hx=2.,hy=2.; //Dimensões em x e y do domínio
    //double hx=Pi,hy=2.;
    int h_level = 0;
    int nx=h_level+1 ,ny=h_level+1; //Número de nos em x  y
    
    TPZVec<REAL> h_s(3,0);
    h_s[0]=2.,h_s[1]=2.,h_s[2]=2.; //Dimensões em x e y do domínio

    switch(sim_case) {

        case MHMProblem: //Pressure
        {
            int pOrder = 1;
            for (pOrder=2; pOrder<=2; pOrder++){
                for (int it=4; it<=4; it++) {
                    h_level = pow(2., it);

                    std::cout<< " ---- Runnig level = " << h_level << " ------ "<<std::endl;
                    TPZVec<int> n_s(3,0.);
                    n_s[0]=h_level,n_s[1]=h_level;
                    n_s[2]=h_level; //Obs!!
                    //h_s[0]=4*h_s[0];

                    MHMNavierStokesTest  *Test2 = new MHMNavierStokesTest();
                    //Test2->Set3Dmesh();
                    //Test2->SetHdivPlus();
                    //Test2->SetElType(ECube);

                    TPZTransform<STATE> Transf(3,3), InvTransf(3,3);
                    Test2->SetTransform(Transf, InvTransf);

                    REAL rot_x = 5.;
                    REAL rot_z = 44.;
                    REAL rot_y = -85.;
                    rot_z = rot_z*Pi/180.;
                    rot_y = rot_y*Pi/180.;
                    rot_z = rot_z*Pi/180.;

                    //Test2->SetRotation3DMatrix(rot_x,rot_y,rot_z);
                    TPZSimulationData *sim_data= new TPZSimulationData;
                    sim_data->SetInternalOrder(pOrder);
                    sim_data->SetSkeletonOrder(pOrder);
                    sim_data->SetCoarseDivisions(n_s);
                    sim_data->SetDomainSize(h_s);
                    sim_data->SetNInterRefs(0);
                    sim_data->SetViscosity(1.);
                    sim_data->SetBrinkmanCoef(0.); //For Brinkman
                    sim_data->SetNthreads(8);
                    //simdata.SetShapeTest(); // Test for shape functions

                    sim_data->SetOptimizeBandwidthQ(true);
                    //sim_data->SetStaticCondensation(false);
                    sim_data->Set_n_iterations(40);
                    sim_data->Set_epsilon_cor(0.0000001);
                    sim_data->Set_epsilon_res(0.0000001);
                    sim_data->SetPardisoSolver();
                    if(h_level==4){
                        sim_data->ActivatePostProcessing();
                    }
                    if(h_level>=32){
                        //    sim_data->SetPardisoSolver();
                    }

                    sim_data->SetProblemType(TStokesAnalytic::ENavierStokes);
                    sim_data->SetDomainType(TStokesAnalytic::EKovasznay);

                    //Transient parameters:
                    //sim_data->SetTimeTotal(11.);
                    //sim_data->SetTimeStep(0.1);

                    Test2->SetSimulationData(sim_data);
                    Test2->Run();

                }
            }
        }
            break;

        case Vugs2D:
        {

            int pOrder = 2;
            h_level = 1;

            TPZVec<int> n_s(3,0.);
            n_s[0]=h_level,n_s[1]=h_level;
            MHMNavierStokesTest  *Test2 = new MHMNavierStokesTest();
            //Test2->SetElType(ETriangle);

            TPZTransform<STATE> Transf(3,3), InvTransf(3,3);
            Test2->SetTransform(Transf, InvTransf);
            TPZSimulationData *sim_data= new TPZSimulationData;
            sim_data->SetInternalOrder(pOrder);
            sim_data->SetSkeletonOrder(pOrder);
            sim_data->SetCoarseDivisions(n_s);
            sim_data->SetDomainSize(h_s);
            sim_data->SetNInterRefs(0);
            sim_data->SetViscosity(1.);  //
            sim_data->SetBrinkmanCoef(0.); //Material 1 => Stokes
            sim_data->SetPermeability(1.); //Material 2 => Darcy
            sim_data->SetNthreads(8);

            sim_data->SetOptimizeBandwidthQ(true);
            //sim_data->SetStaticCondensation(false);
            sim_data->Set_n_iterations(40);
            sim_data->Set_epsilon_cor(0.0000001);
            sim_data->Set_epsilon_res(0.00001);
            sim_data->SetPardisoSolver();
            sim_data->ActivatePostProcessing();

            sim_data->SetProblemType(TStokesAnalytic::EStokes);
            sim_data->SetDomainType(TStokesAnalytic::EVugs2D);

            //Transient parameters:
            //sim_data->SetTimeTotal(11.);
            //sim_data->SetTimeStep(0.1);

            Test2->SetSimulationData(sim_data);
            Test2->Run();

        }
            break;

        case Coupling:
        {

            int pOrder = 1;
            for (pOrder=2; pOrder<=2; pOrder++){
                for (int it=6; it<=6; it++) {
                    h_level = pow(2., it);

                    std::cout<< " ---- Runnig level = " << h_level << " ------ "<<std::endl;
                    TPZVec<int> n_s(3,0.);
                    n_s[0]=h_level,n_s[1]=h_level;
                    MHMNavierStokesTest  *Test2 = new MHMNavierStokesTest();
                    Test2->SetElType(ETriangle);

                    TPZTransform<STATE> Transf(3,3), InvTransf(3,3);
                    Test2->SetTransform(Transf, InvTransf);
                    TPZSimulationData *sim_data= new TPZSimulationData;
                    sim_data->SetInternalOrder(pOrder);
                    sim_data->SetSkeletonOrder(pOrder);
                    sim_data->SetCoarseDivisions(n_s);
                    sim_data->SetDomainSize(h_s);
                    sim_data->SetNInterRefs(0);
                    sim_data->SetViscosity(1.);  //
                    sim_data->SetBrinkmanCoef(0.); //Material 1 => Stokes
                    sim_data->SetPermeability(1.); //Material 2 => Darcy
                    sim_data->SetNthreads(8);

                    sim_data->SetOptimizeBandwidthQ(true);
                    //sim_data->SetStaticCondensation(false);
                    sim_data->Set_n_iterations(40);
                    sim_data->Set_epsilon_cor(0.0000001);
                    sim_data->Set_epsilon_res(0.00001);
                    sim_data->SetPardisoSolver();
                    if(h_level==16){
                        sim_data->ActivatePostProcessing();
                    }

                    sim_data->SetProblemType(TStokesAnalytic::ENavierStokes);
                    sim_data->SetDomainType(TStokesAnalytic::ECouplingSD);

                    //Transient parameters:
                    //sim_data->SetTimeTotal(11.);
                    //sim_data->SetTimeStep(0.1);

                    Test2->SetSimulationData(sim_data);
                    Test2->Run();

                }
            }

        }
            break;

        case HybridProblem:
        {
            int pOrder = 2;

            for (int it=1; it<=1; it++) {
                h_level = pow(2,it+3);

//                2 << (it+3);
                h_level = 1;
                std::cout<< " ---- Runnig level = " << h_level << " ------ "<<std::endl;

                TPZVec<int> n_s(3,0.);
                n_s[0]=h_level ,n_s[1]=h_level;

                n_s[2]=h_level; //Obs!!

                NavierStokesTest  * Test2 = new NavierStokesTest();
                //Test2->Set3Dmesh();
                //Test2->SetElType(ETriangle);
                //Test2->SetHdivPlus();
                Test2->SetFluxOrder(pOrder);
                Test2->SetHdivPlus(0);
                Test2->SetTractionOrder(pOrder-1);
                Test2->SetInternRef(1);

                //Simulation Data
                TPZSimulationData *sim_data= new TPZSimulationData;
                sim_data->SetViscosity(.1);

                sim_data->SetNthreads(0);
                sim_data->SetOptimizeBandwidthQ(true);
                sim_data->Set_n_iterations(20);
                sim_data->Set_epsilon_cor(0.0002);
                sim_data->Set_epsilon_res(0.0002);
                sim_data->SetPardisoSolver();
                sim_data->ActivatePostProcessing();

                if(h_level==64&&pOrder==3){
                    sim_data->ActivatePostProcessing();
                }
                if(h_level>=128){
                    sim_data->SetPardisoSolver();
                }

                //Select problem type (ENavierStokes,ENavierStokesCDG, EOseen,EStokes,EBrinkman)
                sim_data->SetProblemType(TStokesAnalytic::ENavierStokes);
                //Select domain type (EObstacles,EOneCurve,ERetangular,EPconst,EKovasznay,EKovasznayCDG)
                sim_data->SetDomainType(TStokesAnalytic::EObstacles);

                //Transient parameters:
                sim_data->SetTimeTotal(1);
                sim_data->SetTimeStep(0.1);


                //sim_data->ActivatePostProcessing();
                Test2->SetSimulationData(sim_data);

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
                Test2->Run(pOrder, n_s, h_s);

            }
        }
            break;

        default: {
            DebugStop();
        }
    }
}

