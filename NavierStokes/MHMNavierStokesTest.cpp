/*
 *  MHMStokesTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#include "MHMNavierStokesTest.h"
#include "pzcheckgeom.h"
#include "pzstack.h"
#include "TPZParSkylineStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include "TPZInterfaceInsertion.h"
#include "pzinterpolationspace.h"
#include "pzcompel.h"
#include "TPZVecL2.h"
#include "pzintel.h"
#include "TPZNullMaterial.h"
#include "TPZGenGrid2D.h"
#include "TPZLagrangeMultiplier.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "TPZExtendGridDimension.h"
#include "tpzgeoelrefpattern.h"
#include "TPZMHMNavierStokesMeshControl.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"

using namespace std;

const REAL Pi=M_PI;

const REAL phi_r = 0.;

TPZTransform<REAL> MHMNavierStokesTest::f_T(3,3);

TPZTransform<REAL> MHMNavierStokesTest::f_InvT(3,3);

MHMNavierStokesTest::MHMNavierStokesTest()
{
    
    fdim=2; //Dimensão do problema
    fmatID=1; //Materia do elemento volumétrico
    
    //Materiais das condições de contorno
    fmatBCbott=-1;
    fmatBCtop=-2;
    fmatBCleft=-3;
    fmatBCright=-4;
    fmatBCtop_z=-5; //3D
    fmatBCbott_z=-6; //3D normal negativa

    fmatBChole=-7;

    fmatBCleft_2=-23;
    fmatBCright_2=-24;

    //Material do elemento de interface
    fmatLambda=3; // Multiplier material
//    fmatLambdaBC=3;
    
    fmatLambdaBC_bott=-11;
    fmatLambdaBC_top=-12;
    fmatLambdaBC_left=-13;
    fmatLambdaBC_right=-14;
    fmatLambdaBC_top_z=-15;
    fmatLambdaBC_bott_z=-16;

    fmatLambdaBC_hole=-17;

    fmatLambdaBC_left_2=-33;
    fmatLambdaBC_right_2=-34;
    
    fmatInterfaceLeft=5;
    fmatInterfaceRight=6;
    fmatWrap = 7;
    
    //Materia de um ponto
    fmatPoint=-15;
    
    //Condições de contorno do problema

    fdirichlet_v=0;
    fneumann_v=1;

    fdirichlet_sigma=1;
    fneumann_sigma=0;

    f_BJS_condition =7;

    fpenetration=2;
    fpointtype=5;
    fdirichletvar=4;
    
    
    fquadmat1=1; //Parte inferior do quadrado
    fquadmat2=2; //Parte superior do quadrado
    fquadmat3=3; //Material de interface
    
    fviscosity=1.;
    fpermeability=1.;
    ftheta=-1.;
    
    fphi_r=0;
    
    f_is_hdivFull = false;
    
    f_hdivPlus = false;
    
    feltype = EQuadrilateral;

    f_HoleCoord.clear();
    f_ArcCentralNode.clear();

    f_problemtype = TStokesAnalytic::EStokes;

    f_domaintype = TStokesAnalytic::ENone;

    f_ExactSol.fProblemType = f_problemtype;

    f_ExactSol_2.fProblemType = TStokesAnalytic::EBrinkman;

    f_mesh_vector.resize(4);
    
    f_T = TPZTransform<>(3,3);
    f_InvT = TPZTransform<>(3,3);
    
}

MHMNavierStokesTest::~MHMNavierStokesTest()
{
    
}

void MHMNavierStokesTest::Run()
{
    int int_order = f_sim_data->GetInternalOrder();
    int skeleton_order = f_sim_data->GetSkeletonOrder();
    TPZVec<int> n_s = f_sim_data->GetCoarseDivisions();
    TPZVec<REAL> h_s = f_sim_data->GetDomainSize();
    int nrefs = f_sim_data->GetNInterRefs();
    SetProblemType(f_sim_data->GetProblemType());
    SetDomainType(f_sim_data->GetDomainType());
    f_ExactSol.fvisco = f_sim_data->GetViscosity();
    f_ExactSol.fcBrinkman = f_sim_data->GetBrinkmanCoef();

    if(feltype==ECube||feltype==EPrisma||feltype==ETetraedro){
        Set3Dmesh();
    }
    //Gerando malha geométrica com elementos coarse, os quais serão subdomínios MHM
    TPZGeoMesh *gmesh;

    switch(f_domaintype) {

        case TStokesAnalytic::EObstacles: //Pressure
        {
            //    gmesh = CreateGMeshCurve();
            gmesh = CreateGMeshRefPattern(n_s,h_s);
        }
            break;

        case TStokesAnalytic::EVugs2D: //Pressure
        {
            //mesh = CreateGMeshVugsRefPattern(n_s,h_s);
            gmesh = CreateGmshMesh();
            //gmesh = CreateGMesh3D(n_s, h_s);
        }
            break;

        case TStokesAnalytic::ECouplingSD: //Pressure
        case TStokesAnalytic::ECouplingNSD:
        {
            gmesh = CreateGMeshCoupling(n_s, h_s);
        }
            break;

        default: {
            if (f_3Dmesh) {
                gmesh = CreateGMesh3D(n_s, h_s);
            } else {
                gmesh = CreateGMesh(n_s, h_s);
            }
        }
    }

    //Vetor com os indices dos elementos coarse
    TPZVec<int64_t> coarseindex;
    GetElIndexCoarseMesh(gmesh, coarseindex);
    
    //Refinamento de subelemntos
    SubdomainRefine(nrefs,gmesh,coarseindex);
    //InsertLowerDimMaterial(gmesh);

    
    //Criando objeto para gerenciar a malha MHM
    TPZAutoPointer<TPZGeoMesh> gmeshpointer(gmesh);
    TPZMHMNavierStokesMeshControl *StokesControl;
    StokesControl = new TPZMHMNavierStokesMeshControl(gmeshpointer);

    StokesControl->SetStaticCondensation(f_sim_data->IsStaticCondensedQ());
    StokesControl->DefinePartitionbyCoarseIndices(coarseindex); //Define the MHM partition by the coarse element indices

    if(f_domaintype==TStokesAnalytic::EVugs2D){
        StokesControl->SetBJSInterface(true);
    //    InsertArcInterface(StokesControl->GMesh());
    }

    if(0){
        std::ofstream fileg("MalhaGeo_0.txt"); //Impressão da malha geométrica (formato txt)
        std::ofstream filegvtk("MalhaGeo_0.vtk"); //Impressão da malha geométrica (formato vtk)
        StokesControl->GMesh()->Print(fileg);
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk,true);
    }

    std::set<int> matids;
    matids.insert(fmatID);
    if(f_domaintype==TStokesAnalytic::ECouplingSD||f_domaintype==TStokesAnalytic::ECouplingNSD||f_domaintype==TStokesAnalytic::EVugs2D){
        matids.insert(2);
    }
    StokesControl->fMaterialIds = matids;
    matids.clear();
    matids.insert(fmatBCtop);
    matids.insert(fmatBCbott);
    matids.insert(fmatBCleft);
    matids.insert(fmatBCright);
    if (f_3Dmesh) {
        matids.insert(fmatBCtop_z);
        matids.insert(fmatBCbott_z);
    }
    if(f_domaintype==TStokesAnalytic::EObstacles){
        matids.insert(fmatBChole);
    }
    if(f_domaintype==TStokesAnalytic::ECouplingSD||f_domaintype==TStokesAnalytic::ECouplingNSD){
        matids.insert(fmatBCleft_2);
        matids.insert(fmatBCright_2);
    }
    
    StokesControl->fMaterialBCIds = matids;
    StokesControl->fBCTractionMatIds[fmatBCtop]=fmatLambdaBC_top;
    StokesControl->fBCTractionMatIds[fmatBCbott]=fmatLambdaBC_bott;
    StokesControl->fBCTractionMatIds[fmatBCleft]=fmatLambdaBC_left;
    StokesControl->fBCTractionMatIds[fmatBCright]=fmatLambdaBC_right;
    if (f_3Dmesh) {
        StokesControl->fBCTractionMatIds[fmatBCtop_z]=fmatLambdaBC_top_z;
        StokesControl->fBCTractionMatIds[fmatBCbott_z]=fmatLambdaBC_bott_z;
    }
    if(f_domaintype==TStokesAnalytic::EObstacles){
        StokesControl->fBCTractionMatIds[fmatBChole]=fmatLambdaBC_hole;
    }
    if(f_domaintype==TStokesAnalytic::ECouplingSD||f_domaintype==TStokesAnalytic::ECouplingNSD){
        StokesControl->fBCTractionMatIds[fmatBCleft_2]=fmatLambdaBC_left_2;
        StokesControl->fBCTractionMatIds[fmatBCright_2]=fmatLambdaBC_right_2;
    }
    
    InsertMaterialObjects(StokesControl);
    
    StokesControl->SetInternalPOrder(int_order);
    StokesControl->SetSkeletonPOrder(skeleton_order);
    //StokesControl->DivideSkeletonElements(0); //Insere material id do skeleton wrap

    //if (fsimData.GetNInterRefs()>0) {
    StokesControl->SetCoarseAverageMultipliers(true);
    //}

    //Malha computacional
    StokesControl->BuildComputationalMesh(0);

    if(1){
        std::ofstream fileg1("MalhaGeo.txt"); //Impressão da malha geométrica (formato txt)
        std::ofstream filegvtk1("MalhaGeo.vtk"); //Impressão da malha geométrica (formato vtk)
        StokesControl->GMesh()->Print(fileg1);
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, filegvtk1,true);
    }


    if(0){
#ifdef PZDEBUG
    std::ofstream filecm("MalhaC_MHM.txt");
    StokesControl->CMesh()->Print(filecm);

    TPZCompMesh *cmeshP = StokesControl->GetMeshes()[1].operator->();
    std::ofstream outp("Malha_P_MHM.vtk");
    cmeshP->LoadReferences();
    TPZVTKGeoMesh::PrintCMeshVTK(cmeshP, outp, false);
        
    TPZCompMesh *cmeshV = StokesControl->GetMeshes()[0].operator->();
    std::ofstream outv("MalhaC_V_MHM.vtk");
    cmeshV->LoadReferences();
    TPZVTKGeoMesh::PrintCMeshVTK(cmeshV, outv, false);
        
    TPZCompMesh *cmeshM = StokesControl->CMesh().operator->();
    std::ofstream out("MalhaC_MHM.vtk");
    cmeshM->LoadReferences();
    TPZVTKGeoMesh::PrintCMeshVTK(cmeshM, out, false);
#endif
    }


    std::cout << "MHM Hdiv Computational meshes created\n";
    std::cout << "Number of equations MHMStokes " << StokesControl->CMesh()->NEquations() << std::endl;
    std::string configuration;
    
    std::stringstream MHMStokesPref;
    MHMStokesPref << "MHMStokes";

    SolveProblem(StokesControl->CMesh(), StokesControl->GetMeshes(), MHMStokesPref.str());

    //SolveNonLinearProblem(StokesControl->CMesh(), StokesControl->GetMeshes(), MHMStokesPref.str());
    std::cout << "FINISHED!" << std::endl;
    
}




void MHMNavierStokesTest::SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix){

    bool shapetest = f_sim_data->GetShapeTest();
    //calculo solution
    bool shouldrenumber = f_sim_data->GetOptimizeBandwidthQ();
    TPZAnalysis an(cmesh,shouldrenumber);

    if(f_sim_data->IsPardisoSolverQ()){
        TPZSymetricSpStructMatrix strmat(cmesh.operator->());
        strmat.SetNumThreads(f_sim_data->GetNthreads());
        an.SetStructuralMatrix(strmat);
    }else{
        //TPZSkylineStructMatrix strmat(cmesh.operator->());
        //strmat.SetNumThreads(f_sim_data->GetNthreads());
        TPZFStructMatrix strmat(cmesh.operator->());
        strmat.SetNumThreads(f_sim_data->GetNthreads());
        an.SetStructuralMatrix(strmat);
    }

    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();

    // Shape functions plot :
    if(shapetest){
        
        TPZVec<int64_t> equ_indexes(1);
        equ_indexes[0] = 38;
        //        for (int i=37; i<38; i++) {
        //            equ_indexes[i-37] = i;
        //        }
        std::string name_phi = "MHMStokes_shape.vtk";
        TPZVec<std::string> var_name(2);
        var_name[0]="V";
        var_name[1]="P";
        
        //TPZBuildMultiphysicsMesh::ShowShape(f_mesh_vector,cmesh_m, an, name_phi, equ_indexes);
        an.ShowShape(name_phi, equ_indexes, 1, var_name);
        return;
        
    }
    
    std::cout << "Solving\n";
    an.Solve();
    

#ifdef PZDEBUG
    if(0)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Kg = ",global,EMathematicaInput);
        an.Rhs().Print("Fg = ",global,EMathematicaInput);
    }
#endif

    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);

    if(f_sim_data->IsPostProcessingActivatedQ()){

        {
            std::ofstream out(prefix+"_MeshWithSol.txt");
            cmesh->Print(out);
        }


        std::cout << "Post Processing " << std::endl;
        std::string plotfile;
        std::stringstream sout_geo;
        std::stringstream sout;
        {
            sout << prefix << "Approx_";
            ConfigPrint(sout);
            plotfile = sout.str() + "_dim2.vtk";
        }
        plotfile = "StokesMHMPlot.vtk";

        {
            sout_geo << prefix << "Geo_";
            ConfigPrint(sout_geo) << "_dim2.vtk";
        }

        std::ofstream plotfile3(sout_geo.str());
        TPZVTKGeoMesh::PrintGMeshVTK(cmesh.operator->()->Reference(), plotfile3, true);


        std::cout << "plotfiles " << " " << plotfile.c_str() << std::endl;
        TPZStack<std::string> scalnames,vecnames;
        TPZMaterial *mat = cmesh->FindMaterial(1);
        if (!mat) {
            DebugStop();
        }
        else if(mat->NStateVariables() == 1)
        {
            scalnames.Push("P");
            vecnames.Push("V");
            vecnames.Push("f");
            vecnames.Push("V_exact");
            scalnames.Push("P_exact");
            scalnames.Push("Div");
        }

        an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile);
        int resolution = 1;
        an.PostProcess(resolution,cmesh->Dimension());

    }

    if(f_domaintype==TStokesAnalytic::ECavity||f_domaintype==TStokesAnalytic::EObstacles||f_domaintype==TStokesAnalytic::EVugs2D){
        return;
    }

    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,6> Errors;
    ofstream ErroOut("Error_Results_Linear.txt", std::ofstream::app);
    an.SetExact(f_ExactSol.ExactSolution());
    int n_threads_sim =f_sim_data->GetNthreads();
    an.SetThreadsForError(4);
//    an.PostProcessError(Errors,false);

    auto old_buffer = std::cout.rdbuf(nullptr);
    an.PostProcessError(Errors,false);
    std::cout.rdbuf(old_buffer);

    ConfigPrint(ErroOut);
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"Norma L2 - V = "<< Errors[0] << std::endl;
    ErroOut <<"Semi-norma H1/Hdiv - V = "<< Errors[1] << std::endl;
    ErroOut <<"Norma L2 - P = "<< Errors[2] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut <<"Mat2:" << std::endl;
    ErroOut <<"Darcy - Norma L2 - V = "<< Errors[3] << std::endl;
    ErroOut <<"Darcy - Semi-norma H1/Hdiv - V = "<< Errors[4] << std::endl;
    ErroOut <<"Darcy - Norma L2 - P = "<< Errors[5] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();
    
    
}

void MHMNavierStokesTest::SolveNonLinearProblem(TPZAutoPointer<TPZCompMesh> cmesh_m, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix){

    TPZNSAnalysis *NS_analysis = new TPZNSAnalysis;

    TPZVec<std::string> var_name(2);
    var_name[0]="V";
    var_name[1]="P";

    DecomposeType decomposeType = ELU;

    if (f_problemtype==TStokesAnalytic::EStokes) {
        decomposeType = ELDLt;
    }

    int compmesh_size = compmeshes.size();
    TPZManVector<TPZCompMesh *,6> meshVecPtr(compmesh_size);
    for (int i=0; i<compmeshes.size(); i++) {
        meshVecPtr[i] = compmeshes[i].operator->();
    }

    NS_analysis->ConfigureAnalysis(decomposeType, f_sim_data, cmesh_m.operator->(), meshVecPtr, var_name);

    if(f_sim_data->IsTransientQ()){
        NS_analysis->SolveSystemTransient();
    }else{
        NS_analysis->SolveSystem();
    }

    if(f_domaintype==TStokesAnalytic::ECavity||f_domaintype==TStokesAnalytic::EObstacles||f_domaintype==TStokesAnalytic::EVugs2D){
        return;
    }

    //TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh_m);

    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,6> Errors;
    ofstream ErroOut("Error_NavierStokes.txt", std::ofstream::app);
    //papapapapap
    //std::cout << cmesh_m.operator->()->Solution() << std::endl;

    cmesh_m->SolutionN();

    NS_analysis->SetExact(f_ExactSol.ExactSolution());
    int n_threads_sim =f_sim_data->GetNthreads();
    NS_analysis->SetThreadsForError(8);

    auto old_buffer = std::cout.rdbuf(nullptr);
    NS_analysis->PostProcessError(Errors,false);
    std::cout.rdbuf(old_buffer);

#ifdef PZDEBUG
    std::cout <<"-------------" << std::endl;
    std::cout <<"Order = "<< f_sim_data->GetInternalOrder() << "  //  N internal refs = " << f_sim_data->GetNInterRefs() << "  //  Coarse divisions = " << f_sim_data->GetCoarseDivisions()[0] << " x " << f_sim_data->GetCoarseDivisions()[1]  << std::endl;
    std::cout <<"L2-norm - V = "<< Errors[0] << std::endl;
    std::cout  <<"H1/Hdiv semi-norm - V = "<< Errors[1] << std::endl;
    std::cout <<"L2-norm - P = "<< Errors[2] << std::endl;
    if(f_problemtype==TStokesAnalytic::ENavierStokesCDG||f_problemtype==TStokesAnalytic::EOseenCDG){
        std::cout <<"L2-norm - P - CDG formulation = "<< Errors[5] << std::endl;
    }
    std::cout <<"-------------" << std::endl;
#endif

    ConfigPrint(ErroOut);
    ErroOut <<"  //  Order = "<< f_sim_data->GetInternalOrder() << "  //  N internal refs = " << f_sim_data->GetNInterRefs() << "  //  Coarse divisions = " << f_sim_data->GetCoarseDivisions()[0] << " x " << f_sim_data->GetCoarseDivisions()[1]  << std::endl;
    ErroOut <<" " << std::endl;
    //ErroOut <<"Norma H1/HDiv - V = "<< Errors[0] << std::endl;
    ErroOut <<"L2-norm - V = "<< Errors[0] << std::endl;
    ErroOut <<"H1/Hdiv semi-norm - V = "<< Errors[1] << std::endl;
    ErroOut <<"L2-norm - P = "<< Errors[2] << std::endl;
    if(f_problemtype==TStokesAnalytic::ENavierStokesCDG||f_problemtype==TStokesAnalytic::EOseenCDG){
        ErroOut <<"L2-norm - P - CDG formulation = "<< Errors[5] << std::endl;
    }
    ErroOut <<"-------------" << std::endl;
    ErroOut <<"Mat2:" << std::endl;
    ErroOut <<"Darcy - Norma L2 - V = "<< Errors[3] << std::endl;
    ErroOut <<"Darcy - Semi-norma H1/Hdiv - V = "<< Errors[4] << std::endl;
    ErroOut <<"Darcy - Norma L2 - P = "<< Errors[5] << std::endl;
    ErroOut <<"-------------" << std::endl;
    ErroOut.flush();

    std::cout << "FINISHED!" << std::endl;



}


std::ostream &MHMNavierStokesTest::ConfigPrint(std::ostream &out)
{
    int int_order = f_sim_data->GetInternalOrder();
    int skeleton_order = f_sim_data->GetSkeletonOrder();
    TPZVec<int> n_s = f_sim_data->GetCoarseDivisions();
    TPZVec<REAL> h_s = f_sim_data->GetDomainSize();
    int nrefs = f_sim_data->GetNInterRefs();
    
    std::string elemName;
    
    if (feltype==EQuadrilateral) {
        elemName = " Quadrilateral elements : ";
        n_s[2] = 0;
    }else if(feltype==ETriangle){
        elemName = " Triangular elements : ";
        n_s[2] = 0;
    }else if(feltype==ETetraedro){
        elemName = " Tetrahedral elements : ";
    }else if(feltype==ECube){
        elemName = " Cubic elements : ";
    }
    
    out << elemName << n_s[0] <<" x "<< n_s[1] << " x " << n_s[2] << " - N refs : " << nrefs << " - Order Skel : " << skeleton_order << " - Order Intern : " << int_order<<"\n";
    return out;
}


void MHMNavierStokesTest::Rotate(TPZVec<REAL> &co, TPZVec<REAL> &co_r, bool rotate){
    
    if (rotate==true) {
        //rotação +
        co_r[0] = co[0]*cos(phi_r) - co[1]*sin(phi_r);
        co_r[1] = co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }else{
        
        co_r[0] = co[0]*cos(phi_r) + co[1]*sin(phi_r);
        co_r[1] = - co[0]*sin(phi_r) + co[1]*cos(phi_r);
        
    }
    
}

void MHMNavierStokesTest::InsertLowerDimMaterial(TPZGeoMesh *gmesh){
    
    // Inserir elmentos fmatLambda and fmatLambdaBCs

            int64_t nel = gmesh->NElements();
            for (int64_t el = 0; el<nel; el++) {
                TPZGeoEl *gel = gmesh->Element(el);
                if(gel->HasSubElement()&&f_allrefine)
                {
                    continue;
                }
                if (gel->Dimension() != gmesh->Dimension()) {
                    continue;
                }
                int nsides = gel->NSides();
                for (int is = 0; is<nsides; is++) {
                    if (gel->SideDimension(is) != gmesh->Dimension() - 1) {
                        continue;
                    }
                    
                    TPZGeoElSide gelside(gel,is);
                    TPZGeoElSide neighbour = gelside.Neighbour();
                    
                    if (neighbour == gelside && f_allrefine == false) {
                        continue;
                    }
                    

                    
                    while (neighbour != gelside) {
                        if (neighbour.Element()->Dimension() == gmesh->Dimension() - 1) {
                            int neigh_matID = neighbour.Element()->MaterialId();
        
                            if(neigh_matID==fmatBCbott){
                                    TPZGeoElBC(gelside, fmatLambdaBC_bott);
                            }else if(neigh_matID==fmatBCtop){
                                    TPZGeoElBC(gelside, fmatLambdaBC_top);
                            }else if(neigh_matID==fmatBCleft){
                                    TPZGeoElBC(gelside, fmatLambdaBC_left);
                            }else if(neigh_matID==fmatBCright){
                                    TPZGeoElBC(gelside, fmatLambdaBC_right);
                            }else if(f_3Dmesh && neigh_matID==fmatBCbott_z){
                                    TPZGeoElBC(gelside, fmatLambdaBC_bott_z);
                            }else if(f_3Dmesh && neigh_matID==fmatBCtop_z){
                                    TPZGeoElBC(gelside, fmatLambdaBC_top_z);
                            }
        
                            break;
        
                        }
                        if(neighbour.Element()->HasSubElement()){
                            break;
                        }
                        
                        if (gel->LowestFather()->Index()!=neighbour.Element()->LowestFather()->Index()) {
                            break;
                        }
                        
                        neighbour = neighbour.Neighbour();
        
                    }
        
        
                    if (neighbour == gelside) {
                            TPZGeoElBC(gelside, fmatLambda);
                    }
                }
            }

}

bool MHMNavierStokesTest::IsSkellNeighbour(TPZGeoElSide neighbour){

    if (neighbour.Element()->Dimension() == f_mesh0->Dimension()) {
        int nskellneighs = f_skellNeighs.NElements();
    
        for (int iskell = 0; iskell < nskellneighs; iskell++) {
            TPZStack<TPZGeoElSide> sonSides;
            f_skellNeighs[iskell].GetAllSiblings(sonSides);
            for (int ison=0; ison<sonSides.NElements(); ison++) {
                if (neighbour == sonSides [ison]) {
                    return true;
                }
            }
        }
    
    }
    
    return false;
}

TPZGeoMesh *MHMNavierStokesTest::CreateGMesh(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{
    
    int dimmodel = 2;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    x0[0] = -0.5, x0[1] = 0.;
    x1[0] = 1.5, x1[1] = 2.;

    if(f_domaintype==TStokesAnalytic::ECavity){
        x0[0] = 0., x0[1] = 0.;
        x1[0] = 1., x1[1] = 1.;
    }

    if(f_domaintype==TStokesAnalytic::ESinCos||f_domaintype==TStokesAnalytic::ESinCosBDS){
        x0[0] = 0., x0[1] = -1.;
        x1[0] = 2., x1[1] = 1.;
    }

    TPZGenGrid2D grid(n_div,x0,x1);
    
    //grid.SetDistortion(0.5);
    grid.SetRefpatternElements(true);
    if (feltype==ETriangle) {
        grid.SetElementType(MMeshType::ETriangular);
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    grid.Read(gmesh);
    grid.SetBC(gmesh, 4, fmatBCbott);
    grid.SetBC(gmesh, 5, fmatBCright);
    grid.SetBC(gmesh, 6, fmatBCtop);
    grid.SetBC(gmesh, 7, fmatBCleft);
    
    //Save the original mesh
    
    //SetAllRefine();
    
    TPZVec<REAL> centerCo(2,0.);
    centerCo[0]=1.;
    centerCo[1]=0.;
   // UniformRefine(1, gmesh, centerCo, true);

    //UniformRefine2(1, gmesh, n_div);
//    InsertLowerDimMaterial(gmesh);
    SetOriginalMesh(gmesh);
//    UniformRefine2(1, gmesh, n_div);
//    InsertLowerDimMaterial(gmesh);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();

        if(0){
            std::ofstream Dummyfile("GeometricMesh2d.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
        }
    
    return gmesh;
    
}

TPZGeoMesh *MHMNavierStokesTest::CreateGmshMesh()
{
    int64_t id, index;

    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:

    TPZGeoMesh *gmesh = new TPZGeoMesh();

    //Aqui é implementado um método para malhas criadas no GMSH

    //std::string dirname = PZSOURCEDIR;
    std::string meshSource, gmshFolder;
    meshSource = PZSOURCEDIR;
    gmshFolder = "MHMNavierStokes/GmshRefs/Example3D.msh";
    meshSource.replace( meshSource.end()-5, meshSource.end(), gmshFolder);


    int dim = 2;
    if(f_3Dmesh) {
        dim = 3;
    }
//    grid = "/Users/pablocarvalho/Documents/GitHub/geomec_bench/Fase_1/Benchmark1a/gmsh/GeometryBenchP21Original00.msh";

    TPZGmshReader Geometry;
    REAL s = 1.0;
    Geometry.SetCharacteristiclength(s);
    Geometry.GetDimNamePhysical()[dim-1]["bottom"] = fmatBCbott;
    Geometry.GetDimNamePhysical()[dim-1]["right"] = fmatBCright;
    Geometry.GetDimNamePhysical()[dim-1]["top"] = fmatBCtop;
    Geometry.GetDimNamePhysical()[dim-1]["left"] = fmatBCleft;

    if(f_3Dmesh){
        Geometry.GetDimNamePhysical()[dim-1]["top_z"] = fmatBCtop_z;
        Geometry.GetDimNamePhysical()[dim-1]["bottom_z"] = fmatBCbott_z;
    }

    Geometry.GetDimNamePhysical()[dim]["Omega"] = fmatID; // Stokes Vug
    Geometry.GetDimNamePhysical()[dim]["Omega2"] = 2; //Porous Media
    Geometry.SetFormatVersion("4.1");
    gmesh = Geometry.GeometricGmshMesh(meshSource);

    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();

    gmesh->BuildConnectivity();

    //    int n_div = 0;
    //    UniformRefine(gmesh,n_div);

        gmesh->ResetConnectivities();

        //Seting Blend quadrilateral elements:

        int64_t elementid = 0;
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (!gel) {
                delete gmesh->ElementVec()[iel];
                continue;
            }
            TPZManVector<int64_t> nodeindices;
            MElementType elType = gel->Type();
            int matID = gel->MaterialId();
            if (elType == ETriangle) {
                if (gel->HasSubElement()) {
                    DebugStop();
                }
                gel->GetNodeIndices(nodeindices);
                elementid = gel->Id();
                gmesh->ElementVec()[elementid]->SetFatherIndex(-1);
                delete gmesh->ElementVec()[elementid];
                gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoTriangle >>(
                        elementid, nodeindices, matID, *gmesh);

            } else if (elType == EQuadrilateral) {
                if (gel->HasSubElement()) {
                    DebugStop();
                }
                gel->GetNodeIndices(nodeindices);
                elementid = gel->Id();
                gmesh->ElementVec()[elementid]->SetFatherIndex(-1);
                delete gmesh->ElementVec()[elementid];
                gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad >>(
                        elementid, nodeindices, matID, *gmesh);

            } else if (elType == ECube) {
                if (gel->HasSubElement()) {
                    DebugStop();
                }
                gel->GetNodeIndices(nodeindices);
                elementid = gel->Id();
                gmesh->ElementVec()[elementid]->SetFatherIndex(-1);
                delete gmesh->ElementVec()[elementid];
                gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern<pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube >>(
                        elementid, nodeindices, matID, *gmesh);
            }

        }

        gmesh->BuildConnectivity();


    ofstream bf("before.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bf);
    return gmesh;

}


TPZGeoMesh *MHMNavierStokesTest::CreateGMeshCoupling(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{

    int dimmodel = 2;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    x0[0] = 0., x0[1] = -1.;
    x1[0] = Pi, x1[1] = 1.;

    int y_interface =0;
    if(f_domaintype==TStokesAnalytic::ECouplingNSD){
        x0[0] = 0., x0[1] = 0.;
        x1[0] = 1., x1[1] = 2.;
        y_interface = 1.;
    }

    TPZGenGrid2D grid(n_div,x0,x1);

    //grid.SetDistortion(0.5);
    grid.SetRefpatternElements(true);
    if (feltype==ETriangle) {
        grid.SetElementType(MMeshType::ETriangular);
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    grid.Read(gmesh,1);
    TPZGeoMesh *gmesh2 = new TPZGeoMesh;
    grid.ReadAndMergeGeoMesh(gmesh,gmesh2,2);

    grid.SetBC(gmesh, 4, fmatBCbott);
    grid.SetBC(gmesh, 5, fmatBCright);
    grid.SetBC(gmesh, 6, fmatBCtop);
    grid.SetBC(gmesh, 7, fmatBCleft);

    //Save the original mesh

    //SetAllRefine();

    TPZVec<REAL> centerCo(2,0.);
    centerCo[0]=1.;
    centerCo[1]=0.;
    // UniformRefine(1, gmesh, centerCo, true);

    //UniformRefine2(1, gmesh, n_div);
//    InsertLowerDimMaterial(gmesh);
    SetOriginalMesh(gmesh);
//    UniformRefine2(1, gmesh, n_div);
//    InsertLowerDimMaterial(gmesh);

    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    int nel = gmesh->NElements();
    TPZVec<REAL> center_coord(3,0.);
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        int nsides = gel->NSides();
        if(gel->MaterialId()==fmatBCright||gel->MaterialId()==fmatBCleft){
            TPZGeoElSide gelside_bc(gel,nsides-1);
            gelside_bc.CenterX(center_coord);
            if(center_coord[1]<y_interface){
                gel->SetMaterialId(gel->MaterialId()-20);
            }
        }
        if(gel->Dimension()!=2){
            continue;
        }
        TPZGeoElSide gelside(gel,nsides-1);
        gelside.CenterX(center_coord);
        if(center_coord[1]<y_interface){
            gel->SetMaterialId(2);
        }
    }

    gmesh->BuildConnectivity();

    if(0){
        std::ofstream Dummyfile("GeometricMesh2d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }

    return gmesh;

}

TPZGeoMesh *MHMNavierStokesTest::CreateGMesh3D(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{
    
    int dimmodel = 2;
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    x0[0] = 0., x0[1] = -1.;
    x1[0] = 2., x1[1] = 1.;
    
    x0[2] = 0.;
    x1[2] = 2.;
//    x0[2] = -1.;
//    x1[2] = 1.;

    TPZGenGrid2D grid(n_div,x0,x1);
    
    
    if (feltype == ETriangle|| feltype == EPrisma ) {
        grid.SetElementType(MMeshType::ETriangular);
    }
    TPZGeoMesh *gmesh = new TPZGeoMesh;

    //grid.SetDistortion(0.2);
    
    if (feltype==EQuadrilateral||feltype==ECube||feltype==EPrisma||feltype==ETriangle) {
    
        grid.Read(gmesh);
        grid.SetBC(gmesh, 4, fmatBCbott);
        grid.SetBC(gmesh, 5, fmatBCright);
        grid.SetBC(gmesh, 6, fmatBCtop);
        grid.SetBC(gmesh, 7, fmatBCleft);
        
        
        REAL thickness = h_s[2]/n_div[2];
        TPZExtendGridDimension extend(gmesh,thickness);
        int numlayers = n_div[2];
        extend.SetElType(1);///???
        gmesh = extend.ExtendedMesh(numlayers,fmatBCbott_z,fmatBCtop_z);
        gmesh->SetDimension(3);
    
    } else if(feltype==ETetraedro){
        
        gmesh->SetDimension(3);
        int64_t id = 0;
        int nx = n_div[0]+1, ny = n_div[1]+1, nz = n_div[2]+1;
        TPZVec<REAL> coord(3,0.);
        int nnodes=nx*ny*nz;
        gmesh->NodeVec().Resize(nx*ny*nz);
        for (int k=0; k<nz; k++) {
            for(int i = 0; i < ny; i++){
                for(int j = 0; j < nx; j++){
                    id = i*nx + j+ k*nx*ny;
                    coord[0] = (j)*h_s[0]/(nx - 1);
                    coord[1] = (i)*h_s[1]/(ny - 1)-1;
                    coord[2] = (k)*h_s[2]/(nz - 1);
                    gmesh->NodeVec()[id].Initialize(coord, *gmesh);
                }
            }
        }
        
        
        TPZVec<int64_t> nodeindD1(4,0), nodeindD2(4,0), nodeindU1(4,0), nodeindU2(4,0), nodeindL1(4,0), nodeindL2(4,0);
        int64_t index=0;
        
        for(int kq=0; kq<n_div[2]; kq++){
            for(int iq = 0; iq < n_div[1]; iq++){
                for(int jq = 0; jq < n_div[0]; jq++){
                    
        
                    // Plano xy
                    nodeindD1[0] = (iq)*ny + (jq) + kq*nx*ny;
                    nodeindD1[1] = nodeindD1[0]+1;
                    nodeindD1[2] = nodeindD1[0]+nx;
                    nodeindD1[3] = nodeindD1[1] + (1)*nx*ny;
                    gmesh->CreateGeoElement(ETetraedro, nodeindD1, fmatID, index,1);
                    
                    index++;
                    
                    nodeindD2[0] = nodeindD1[1];
                    nodeindD2[1] = nodeindD1[2];
                    nodeindD2[2] = nodeindD1[1]+nx;
                    nodeindD2[3] = nodeindD1[1] + (1)*nx*ny;
                    gmesh->CreateGeoElement(ETetraedro, nodeindD2, fmatID, index,1);
                    
                    index++;
                    
                    nodeindU1[0] = nodeindD1[0] + (1)*nx*ny;
                    nodeindU1[1] = nodeindU1[0]+1;
                    nodeindU1[2] = nodeindU1[0]+nx;
                    nodeindU1[3] = nodeindD1[2];
                    gmesh->CreateGeoElement(ETetraedro, nodeindU1, fmatID, index,1);
                    
                    index++;
                    
                    nodeindU2[0] = nodeindU1[1];
                    nodeindU2[1] = nodeindU1[2];
                    nodeindU2[2] = nodeindU1[1]+nx;
                    nodeindU2[3] = nodeindD1[2];
                    gmesh->CreateGeoElement(ETetraedro, nodeindU2, fmatID, index,1);
                    
                    index++;
                    
                    // Plano xz
                    
                    nodeindL1[0] = nodeindD1[0];
                    nodeindL1[1] = nodeindD1[2];
                    nodeindL1[2] = nodeindL1[0]+nx*ny;
                    nodeindL1[3] = nodeindU1[1];
                    gmesh->CreateGeoElement(ETetraedro, nodeindL1, fmatID, index,1);
                    
                    index++;
                    
                    
                    nodeindL2[0] = nodeindD2[2];
                    nodeindL2[1] = nodeindD1[1]+nx*ny;
                    nodeindL2[2] = nodeindU2[2];
                    nodeindL2[3] = nodeindU2[3];
                    gmesh->CreateGeoElement(ETetraedro, nodeindL2, fmatID, index,1);
                    
                    index++;
                    
                }
            }
        }

        gmesh->BuildConnectivity();
        
        // Boundary Conditions
        const int numelements = gmesh->NElements();
        for(int el=0; el<numelements; el++)
        {
            TPZManVector <TPZGeoNode,4> Nodefinder(4);
            TPZManVector <REAL,3> nodecoord(3);
            TPZGeoEl *tetra = gmesh->ElementVec()[el];
            
            // na face x = 0
            TPZVec<int64_t> ncoordVec(0); int64_t sizeOfVec = 0;
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[0]-x0[0])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCleft);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face x = 1
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[0]-x1[0])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCright);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face y = 0
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                
                
                if (fabs(nodecoord[1]-x0[1])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCbott);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face y = 1
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[1]-x1[1])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCtop);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face z = 0
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[2]-x0[2])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCbott_z);
            }
            
            ncoordVec.clear();
            sizeOfVec = 0;
            // na face z = 1
            for (int i = 0; i < 4; i++)
            {
                int64_t pos = tetra->NodeIndex(i);
                Nodefinder[i] = gmesh->NodeVec()[pos];
                Nodefinder[i].GetCoordinates(nodecoord);
                if (fabs(nodecoord[2]-x1[2])<1.e-5)
                {
                    sizeOfVec++;
                    ncoordVec.Resize(sizeOfVec);
                    ncoordVec[sizeOfVec-1] = pos;
                }
            }
            if(ncoordVec.NElements() == 3)
            {
                int lado = tetra->WhichSide(ncoordVec);
                TPZGeoElSide tetraSide(tetra, lado);
                TPZGeoElBC(tetraSide,fmatBCtop_z);
            }
            
        }
        
    
    }
    //Save the original mesh
    grid.SetRefpatternElements(true);
    
    //SetAllRefine();
    
    TPZVec<REAL> centerCo(2,0.);
    centerCo[0]=1.;
    centerCo[1]=0.;
    // UniformRefine(1, gmesh, centerCo, true);
    
    //UniformRefine2(1, gmesh, n_div);
    
    //InsertLowerDimMaterial(gmesh);
    SetOriginalMesh(gmesh);
    
    //UniformRefine2(1, gmesh, n_div);
    //InsertLowerDimMaterial(gmesh);
    
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();
    
    gmesh->BuildConnectivity();
    
    
    if(0){
        std::ofstream Dummyfile("GeometricMesh3D.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }
    
    return gmesh;
    
}

TPZGeoMesh *MHMNavierStokesTest::CreateGMeshRefPattern(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{

    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGeoMesh *gmeshRef = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmeshRef->SetDimension(2);

    int nquadnods_x=n_div[0]+1;
    int nquadnods_y=n_div[1]+1;
    int nodes = nquadnods_x*nquadnods_y;
    gmesh->SetMaxNodeId(nodes-1);
    gmesh->NodeVec().Resize(nodes);

    gmeshRef->SetMaxNodeId(nodes-1);
    gmeshRef->NodeVec().Resize(nodes);

    TPZManVector<TPZGeoNode,7> Node(nodes);

    TPZManVector<int64_t,6> TopolQTriangle(3);
    TPZManVector<int64_t,6> TopolQQuadrilateral(4);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<REAL,3> coord(3,0.),coordrot(3,0.);
    TPZVec<REAL> xc(3,0.);

    int64_t nodeindex = 0;


    //Exterior coordinates
    for(int i = 0; i < nquadnods_y; i++){
        for(int j = 0; j < nquadnods_x; j++){
            coord[0] = (j)*h_s[0]/(nquadnods_x-1);
            coord[1] = (i)*h_s[1]/(nquadnods_y-1);
            gmesh->NodeVec()[nodeindex].SetCoord(coord);
            gmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }

    //Ponto 1
    int64_t elementid = 0;

    //Tringular elements with GeoBlend:

//    TopolQQuadrilateral[0] =0;
//    TopolQQuadrilateral[1] =n_div[0];
//    TopolQQuadrilateral[2] =n_div[1]*(n_div[0]+1);
//    TopolQQuadrilateral[3] =(n_div[0]+1)*(n_div[1]+1)-1;
//    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQQuadrilateral, fmatID,*gmesh);
//    elementid++;



    //Conectividade dos elementos:

    for(int i = 0; i < n_div[1]; i++){
        for(int j = 0; j < n_div[0]; j++){
            TopolQQuadrilateral[0] = (i)*(n_div[1]+1) + (j);
            TopolQQuadrilateral[1] = TopolQQuadrilateral[0]+1;
            TopolQQuadrilateral[2] = TopolQQuadrilateral[1]+n_div[0]+1;
            TopolQQuadrilateral[3] = TopolQQuadrilateral[0]+n_div[0]+1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQQuadrilateral, fmatID,*gmesh);
            elementid++;
        }
    }

    gmesh->BuildConnectivity();


    //if(f_Holemesh){
        TPZManVector<REAL,6> FirstCoord(3,0.), h_el(3,0.);
        h_el[0]=h_s[0]/n_div[0];
        h_el[1]=h_s[1]/n_div[1];

        int nreft = 3;
        TPZVec<TPZGeoEl *> sons1;
        //for (int iref = 0; iref < nreft; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }

            MElementType elType = gel->Type();
            if(gel->MaterialId()==1&&elType==EQuadrilateral&&iel==0){

                TPZAutoPointer<TPZRefPattern> RefpObst = CreateGMeshObstacle(nreft,FirstCoord,h_el);
                gel->SetRefPattern(RefpObst);
                gel->Divide(sons1);

            }else{
                TPZAutoPointer<TPZRefPattern> uniform = CreateGMeshQuadRef(nreft,h_el);
                gel->SetRefPattern(uniform);
                gel->Divide(sons1);
            }
        }
        //}

  //  }

    int nref1 = 0;
    TPZVec<TPZGeoEl *> sons2;
    for (int iref = 0; iref < nref1; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (!gel) {
                continue;
            }
            if (gel->HasSubElement()) {
                continue;
            }
            if (gel->MaterialId()==1){
                gel->Divide(sons2);
            }
        }
    }

    // BC and Hole elements
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if (gel->HasSubElement()) {
            continue;
        }
        if (gel->MaterialId()==1){

            //Create Hole elements:
            TPZGeoElSide gelside(gel,7);
            TPZGeoElSide neighbour = gelside.Neighbour();

            bool near_hole = true;
            while (neighbour != gelside) {
                if(neighbour.Element()->Dimension()==2){
                    near_hole = false;
                }
                neighbour = neighbour.Neighbour();
            }

            if (near_hole) {
                TPZGeoElBC(gelside, fmatBChole);
            }


            TPZManVector<int64_t> TopolPlate(4);

            int64_t totalnodes = gel->NNodes();
            for (int i=0; i<totalnodes; i++){
                TopolPlate[i] = gel->NodeIndex(i);
            }

            //Set BC elements:
            TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
            TPZManVector <REAL,3> nodecoord(3);

            TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
            TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;
            TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
            TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;

            for (int64_t i = 0; i < totalnodes; i++)
            {
                Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
                Nodefinder[i].GetCoordinates(nodecoord);

                if (nodecoord[1] == 0.)
                {
                    sizeOfbottVec++;
                    ncoordzbottVec.Resize(sizeOfbottVec);
                    ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
                }
                if (nodecoord[1] == h_s[1])
                {
                    sizeOftopVec++;
                    ncoordztopVec.Resize(sizeOftopVec);
                    ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == 0.)
                {
                    sizeOfleftVec++;
                    ncoordzleftVec.Resize(sizeOfleftVec);
                    ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == h_s[0])
                {
                    sizeOfrightVec++;
                    ncoordzrightVec.Resize(sizeOfrightVec);
                    ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
                }
            }

            if (sizeOfbottVec == 2) {
                int sidesbott = gel->WhichSide(ncoordzbottVec);
                TPZGeoElSide platesidebott(gel, sidesbott);
                TPZGeoElBC(platesidebott,fmatBCbott);
            }
            if (sizeOftopVec == 2) {
                int sidestop = gel->WhichSide(ncoordztopVec);
                TPZGeoElSide platesidetop(gel, sidestop);
                TPZGeoElBC(platesidetop,fmatBCtop);
            }
            if (sizeOfleftVec == 2) {
                int sidesleft = gel->WhichSide(ncoordzleftVec);
                TPZGeoElSide platesideleft(gel, sidesleft);
                TPZGeoElBC(platesideleft,fmatBCleft);
            }
            if (sizeOfrightVec == 2) {
                int sidesright = gel->WhichSide(ncoordzrightVec);
                TPZGeoElSide platesideright(gel, sidesright);
                TPZGeoElBC(platesideright,fmatBCright);
            }
        }
    }


    TPZManVector<REAL,3> coordFather(3,0.);
    //Remove filiation:
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(gel->HasSubElement()) {
            int nsubel = gel->NSubElements();
            for (int isub=0; isub<nsubel; isub++) {
                //set center coord of father in subelements
                int sub_index = gel->SubElement(isub)->Index();
                TPZManVector<REAL> xcenter(3,0.);
                TPZManVector<REAL,3> xicenter(gel->Dimension(),0.);
                gel->CenterPoint(gel->NSides()-1, xicenter);
                gel->X(xicenter,xcenter);
                f_HoleCoord[sub_index] = xcenter;

                int lowest_index = gel->LowestFather()->Index();

                //Verifiy if father of father index was created
                if (f_HoleCoord.find(lowest_index)!=f_HoleCoord.end()) {
                    f_HoleCoord[sub_index] = f_HoleCoord[lowest_index];
                }


                gel->SetSubElement(isub, 0);
                gel->SetFatherIndex(-1);
            }
            gmesh->DeleteElement(gel);
        }
    }

    //Verify filiation removal
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            continue;
        }
        if (gel->HasSubElement()) {
            DebugStop();
        }
    }

    //Build Arcs in BC holes
    int nnodes = 0;
    TPZManVector<REAL,3> coord0(3,0.),coord1(3,0.),coordHole(3,0),coordNode(3,0);
    for(int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            continue;
        }
        if(gel->MaterialId()==fmatBChole){

            TPZManVector<int64_t> nodeindices;
            TPZManVector<int64_t,3> TopolArc(3);
            gel->GetNodeIndices(nodeindices);
            TopolArc[0]=nodeindices[0];
            TopolArc[1]=nodeindices[1];


            //Add new node:
            nnodes = gmesh->NNodes();
            gmesh->NodeVec().Resize(nnodes+1);
            nodeindex = nnodes;
            gmesh->SetNodeIdUsed(nodeindex);
            REAL radius = h_s[0]/(4.*n_div[0]);

            //Find hole central coord
            int nsides = gel->NSides();

            TPZGeoElSide gelside(gel,nsides-1);
            TPZGeoElSide neighbour = gelside.Neighbour();

            while (neighbour != gelside) {
                if(neighbour.Element()->Dimension()==2){
                    int index_neigh = neighbour.Element()->Index();
                    coordHole = f_HoleCoord[index_neigh];
                    f_ArcCentralNode[index_neigh] = nodeindex;
                    break;
                }
                neighbour = neighbour.Neighbour();
            }

            gel->Node(0).GetCoordinates(coord0);
            gel->Node(1).GetCoordinates(coord1);

            REAL rel = (coord0[1]-coordHole[1])/radius;
            REAL theta0H = asin(rel);
            rel = (coord1[1]-coordHole[1])/radius;
            REAL theta1H = asin(rel);
            REAL theta = theta0H - (theta0H-theta1H)/2;

            if (fabs(coord1[0]-coordHole[0])<1.e-6) {
                coord1[0] = coordHole[0];
            }
            if (fabs(coord0[1]-coordHole[1])<1.e-6) {
                coord0[1] = coordHole[1];
            }

            if(coord0[0]<coordHole[0]||coord1[0]<coordHole[0]){
                theta = Pi - theta;
            }

            coordNode = ParametricCircle(radius, theta);
            coordNode[0] +=coordHole[0];
            coordNode[1] +=coordHole[1];
            gmesh->NodeVec()[nodeindex].SetCoord(coordNode);
            gmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);

            TopolArc[2]=nodeindex;

            elementid = gel->Id();

            delete gmesh->ElementVec()[elementid];
            gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBChole,*gmesh);
        }
    }

    gmesh->ResetConnectivities();

    //Seting Blend quadrilateral elements:


    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            delete gmesh->ElementVec()[iel];
            continue;
        }
        TPZManVector<int64_t> nodeindices;
        MElementType elType = gel->Type();
        if(gel->MaterialId()==1&&elType==EQuadrilateral){
            if (gel->HasSubElement()) {
                DebugStop();
            }
            gel->GetNodeIndices(nodeindices);
            elementid = gel->Id();
            gmesh->ElementVec()[elementid]->SetFatherIndex(-1);
            delete gmesh->ElementVec()[elementid];
            gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad >> (elementid, nodeindices, fmatID,*gmesh);

        }
    }


    gmesh->BuildConnectivity();


    int nref2 = 0; //refinamentos internos
    TPZVec<TPZGeoEl *> sons3;
    for (int iref = 0; iref < nref2; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (!gel) {
                continue;
            }
            if (gel->HasSubElement()) {
                continue;
            }
            //if(gel->MaterialId()==1){
            gel->Divide(sons3);
            //}
        }
    }

    gmesh->BuildConnectivity();
    //InsertLowerDimMaterial(gmesh);
    //SetOriginalMesh(gmesh);
    //gmesh->BuildConnectivity();
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();



    {
        std::ofstream Dummyfile("GeometricMesh2d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }

    return gmesh;


}


TPZGeoMesh *MHMNavierStokesTest::CreateGMeshVugsRefPattern(TPZVec<int> &n_div, TPZVec<REAL> &h_s)
{

    TPZGeoMesh *gmesh = new TPZGeoMesh;
    TPZGeoMesh *gmeshRef = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gmeshRef->SetDimension(2);

    int nquadnods_x=n_div[0]+1;
    int nquadnods_y=n_div[1]+1;
    int nodes = nquadnods_x*nquadnods_y;
    gmesh->SetMaxNodeId(nodes-1);
    gmesh->NodeVec().Resize(nodes);

    gmeshRef->SetMaxNodeId(nodes-1);
    gmeshRef->NodeVec().Resize(nodes);

    TPZManVector<TPZGeoNode,7> Node(nodes);

    TPZManVector<int64_t,6> TopolQTriangle(3);
    TPZManVector<int64_t,6> TopolQQuadrilateral(4);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<REAL,3> coord(3,0.),coordrot(3,0.);
    TPZVec<REAL> xc(3,0.);

    int64_t nodeindex = 0;


    //Exterior coordinates
    for(int i = 0; i < nquadnods_y; i++){
        for(int j = 0; j < nquadnods_x; j++){
            coord[0] = (j)*h_s[0]/(nquadnods_x-1);
            coord[1] = (i)*h_s[1]/(nquadnods_y-1);
            gmesh->NodeVec()[nodeindex].SetCoord(coord);
            gmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }

    //Ponto 1
    int64_t elementid = 0;

    //Conectividade dos elementos:

    for(int i = 0; i < n_div[1]; i++){
        for(int j = 0; j < n_div[0]; j++){
            TopolQQuadrilateral[0] = (i)*(n_div[1]+1) + (j);
            TopolQQuadrilateral[1] = TopolQQuadrilateral[0]+1;
            TopolQQuadrilateral[2] = TopolQQuadrilateral[1]+n_div[0]+1;
            TopolQQuadrilateral[3] = TopolQQuadrilateral[0]+n_div[0]+1;
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQQuadrilateral, fmatID,*gmesh);
            elementid++;
        }
    }

    gmesh->BuildConnectivity();


    //if(f_Holemesh){
    TPZManVector<REAL,6> FirstCoord(3,0.), h_el(3,0.);
    h_el[0]=h_s[0]/n_div[0];
    h_el[1]=h_s[1]/n_div[1];

    int nreft = 4; //number of refinements (perimetral direction)
    TPZVec<TPZGeoEl *> sons1;
    //for (int iref = 0; iref < nreft; iref++) {
    int nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if (gel->HasSubElement()) {
            continue;
        }

        MElementType elType = gel->Type();
        if(gel->MaterialId()==1&&elType==EQuadrilateral&&iel==0){

            TPZAutoPointer<TPZRefPattern> RefpObst = CreateGMeshVugs(nreft,FirstCoord,h_el);
            gel->SetRefPattern(RefpObst);
            gel->Divide(sons1);

        }else{
            DebugStop();
            TPZAutoPointer<TPZRefPattern> uniform = CreateGMeshQuadRef(nreft,h_el);
            gel->SetRefPattern(uniform);
            gel->Divide(sons1);
        }
    }
    //}

    //  }

    //Set MatIds for Vugs:
    int nelem = gmesh->NElements();
        for (int iel = 0; iel < nelem; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];

            if (gel->HasSubElement()) {
                continue;
            }

            TPZManVector<REAL> xcenter(3, 0.);
            TPZManVector<REAL, 3> xicenter(gel->Dimension(), 0.);
            gel->CenterPoint(gel->NSides() - 1, xicenter);
            gel->X(xicenter, xcenter);
            if((xcenter[0]-1.)*(xcenter[0]-1.)+(xcenter[1]-1)*(xcenter[1]-1)>0.5*0.5){
                gel->SetMaterialId(2);
            }

        }

    int matSkeleton = 4;
    // BC and Hole elements
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if (gel->HasSubElement()) {
            continue;
        }

        // External boundary conditions
        if (gel->MaterialId()==2){

            TPZManVector<int64_t> TopolPlate(4);

            int64_t totalnodes = gel->NNodes();
            for (int i=0; i<totalnodes; i++){
                TopolPlate[i] = gel->NodeIndex(i);
            }

            //Set BC elements:
            TPZManVector <TPZGeoNode> Nodefinder(totalnodes);
            TPZManVector <REAL,3> nodecoord(3);

            TPZVec<int64_t> ncoordzbottVec(0); int64_t sizeOfbottVec = 0;
            TPZVec<int64_t> ncoordztopVec(0); int64_t sizeOftopVec = 0;
            TPZVec<int64_t> ncoordzleftVec(0); int64_t sizeOfleftVec = 0;
            TPZVec<int64_t> ncoordzrightVec(0); int64_t sizeOfrightVec = 0;

            for (int64_t i = 0; i < totalnodes; i++)
            {
                Nodefinder[i] = gmesh->NodeVec()[TopolPlate[i]];
                Nodefinder[i].GetCoordinates(nodecoord);

                if (nodecoord[1] == 0.)
                {
                    sizeOfbottVec++;
                    ncoordzbottVec.Resize(sizeOfbottVec);
                    ncoordzbottVec[sizeOfbottVec-1] = TopolPlate[i];
                }
                if (nodecoord[1] == h_s[1])
                {
                    sizeOftopVec++;
                    ncoordztopVec.Resize(sizeOftopVec);
                    ncoordztopVec[sizeOftopVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == 0.)
                {
                    sizeOfleftVec++;
                    ncoordzleftVec.Resize(sizeOfleftVec);
                    ncoordzleftVec[sizeOfleftVec-1] = TopolPlate[i];
                }
                if (nodecoord[0] == h_s[0])
                {
                    sizeOfrightVec++;
                    ncoordzrightVec.Resize(sizeOfrightVec);
                    ncoordzrightVec[sizeOfrightVec-1] = TopolPlate[i];
                }
            }

            if (sizeOfbottVec == 2) {
                int sidesbott = gel->WhichSide(ncoordzbottVec);
                TPZGeoElSide platesidebott(gel, sidesbott);
                TPZGeoElBC(platesidebott,fmatBCbott);
            }
            if (sizeOftopVec == 2) {
                int sidestop = gel->WhichSide(ncoordztopVec);
                TPZGeoElSide platesidetop(gel, sidestop);
                TPZGeoElBC(platesidetop,fmatBCtop);
            }
            if (sizeOfleftVec == 2) {
                int sidesleft = gel->WhichSide(ncoordzleftVec);
                TPZGeoElSide platesideleft(gel, sidesleft);
                TPZGeoElBC(platesideleft,fmatBCleft);
            }
            if (sizeOfrightVec == 2) {
                int sidesright = gel->WhichSide(ncoordzrightVec);
                TPZGeoElSide platesideright(gel, sidesright);
                TPZGeoElBC(platesideright,fmatBCright);
            }
        }
    }


    TPZManVector<REAL,3> coordFather(3,0.);
    //Remove filiation:
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(gel->HasSubElement()) {
            int nsubel = gel->NSubElements();
            for (int isub=0; isub<nsubel; isub++) {
                //set center coord of father in subelements
                int sub_index = gel->SubElement(isub)->Index();
                TPZManVector<REAL> xcenter(3,0.);
                TPZManVector<REAL,3> xicenter(gel->Dimension(),0.);
                gel->CenterPoint(gel->NSides()-1, xicenter);
                gel->X(xicenter,xcenter);
                f_HoleCoord[sub_index] = xcenter;

                int lowest_index = gel->LowestFather()->Index();

                //Verifiy if father of father index was created
                if (f_HoleCoord.find(lowest_index)!=f_HoleCoord.end()) {
                    f_HoleCoord[sub_index] = f_HoleCoord[lowest_index];
                }


                gel->SetSubElement(isub, 0);
                gel->SetFatherIndex(-1);
            }
            gmesh->DeleteElement(gel);
        }
    }

    //Verify filiation removal
    nel = gmesh->NElements();
    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            continue;
        }
        if (gel->HasSubElement()) {
            DebugStop();
        }
    }



    gmesh->ResetConnectivities();

    //Seting Blend quadrilateral elements:


    for (int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            delete gmesh->ElementVec()[iel];
            continue;
        }
        TPZManVector<int64_t> nodeindices;
        MElementType elType = gel->Type();
        int matID = gel->MaterialId();
        if(elType==EQuadrilateral){
            if (gel->HasSubElement()) {
                DebugStop();
            }
            gel->GetNodeIndices(nodeindices);
            elementid = gel->Id();
            gmesh->ElementVec()[elementid]->SetFatherIndex(-1);
            delete gmesh->ElementVec()[elementid];
            gmesh->ElementVec()[elementid] = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad >> (elementid, nodeindices, matID,*gmesh);

        }
    }


    gmesh->BuildConnectivity();


    int nref2 = 0; //refinamentos internos
    TPZVec<TPZGeoEl *> sons3;
    for (int iref = 0; iref < nref2; iref++) {
        int nel = gmesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = gmesh->ElementVec()[iel];
            if (!gel) {
                continue;
            }
            if (gel->HasSubElement()) {
                continue;
            }
            //if(gel->MaterialId()==1){
            gel->Divide(sons3);
            //}
        }
    }

    gmesh->BuildConnectivity();
    //InsertLowerDimMaterial(gmesh);
    //SetOriginalMesh(gmesh);
    //gmesh->BuildConnectivity();
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();



    {
        std::ofstream Dummyfile("GeometricMesh2d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    }

    return gmesh;


}


TPZAutoPointer<TPZRefPattern> MHMNavierStokesTest::CreateGMeshObstacle(int nrefs, TPZManVector<REAL,6> &FirstCoord, TPZManVector<REAL,6> &h_el)
{

    TPZGeoMesh geomesh;
    geomesh.SetDimension(2);

    int nrefs_c = 16;

    int nquadnods=pow(2.,1+nrefs)+1;
    int n_ext_nds = 8*pow(2,nrefs);
    int n_cir_nds = 16 * pow(2, nrefs);
    int nodes = n_ext_nds+n_cir_nds*(1+nrefs_c);

    REAL radius = h_el[0]/4.;
    geomesh.SetMaxNodeId(nodes-1);
    geomesh.NodeVec().Resize(nodes);
    //TPZManVector<TPZGeoNode,7> Node(nodes);

    TPZManVector<int64_t,8> TopolQuadrilateral(4);

    TPZManVector<REAL,3> coord(3,0.);

    TPZVec<REAL> xc(3,0.);

    REAL hx = h_el[0];
    REAL hy = h_el[1];

    int64_t nodeindex_r = 0,nodeindex = 0,  index_j1 = 0, index_j2 = 0;


    // External nodes for father
    TopolQuadrilateral[0] = n_ext_nds-nquadnods-(nquadnods-3)/2;
    TopolQuadrilateral[1] = TopolQuadrilateral[0]+nquadnods-1;
    TopolQuadrilateral[3] = n_ext_nds-1-2*(nquadnods-1)-(nquadnods-3)/2;
    TopolQuadrilateral[2] = TopolQuadrilateral[3]+1-nquadnods;

    //Exterior coordinates
    for(int i = 0; i < nquadnods; i++){
        for(int j = 0; j < nquadnods; j++){
            coord[0] = (j)*hx/(nquadnods - 1);
            coord[1] = (i)*hy/(nquadnods - 1);
            if((coord[0]>0&&coord[0]<hx)&&(coord[1]>0&&coord[1]<hy)){
                continue;
            }
            if(i==0){
                nodeindex_r=TopolQuadrilateral[0]+j;
            }else if(i>0&&i<nquadnods-1&&j==0){
                nodeindex_r=TopolQuadrilateral[0]-i;
            }else if(i>0&&i<nquadnods-1&&j==nquadnods-1&&coord[1]<hy/2.){
                index_j2++;
                nodeindex_r=TopolQuadrilateral[1]+index_j2;
            }else if(i>0&&i<nquadnods-1&&j==nquadnods-1&&coord[1]>hy/2.){
                index_j1++;
                nodeindex_r=index_j1;
            }else if(i==nquadnods-1){
                nodeindex_r=TopolQuadrilateral[3]-j;
            }else{
                nodeindex_r=0;
            }

            geomesh.NodeVec()[nodeindex_r].SetCoord(coord);
            geomesh.NodeVec()[nodeindex_r].SetNodeId(nodeindex);
        }
    }
    nodeindex = n_ext_nds;
    //Circunference coordinates
    REAL refradius = 0.;
    REAL circle_pos = hx/2.;
    REAL hipotenusa =0.;
    REAL theta_c = 0.,aux_theta =0.;
    //nrefs = 6;

    for (int iref = 0; iref <= nrefs_c; iref++) {
        for (int inode = 0; inode < n_cir_nds ; inode++) {
            // i node
            refradius = radius+iref*((circle_pos-radius)/(1+nrefs_c));
            theta_c = inode *2.*M_PI/n_cir_nds;
            if(iref>0) {
                aux_theta = theta_c;
                if(theta_c>M_PI/4.&&theta_c<=M_PI/2.) {
                    aux_theta = M_PI / 2. - theta_c;
                }
                if(theta_c>M_PI/2&&theta_c<=3.*M_PI/4.) {
                    aux_theta = theta_c - M_PI / 2.;
                }
                if(theta_c>3.*M_PI/4.&&theta_c<=M_PI) {
                    aux_theta = M_PI - theta_c;
                }
                if(theta_c>M_PI&&theta_c<=5.*M_PI/4.) {
                    aux_theta = theta_c - M_PI ;
                }
                if(theta_c>5.*M_PI/4.&&theta_c<=3.*M_PI/2.) {
                    aux_theta = 3.*M_PI/2.-theta_c;
                }
                if(theta_c>3.*M_PI/2.&&theta_c<=7.*M_PI/4) {
                    aux_theta = theta_c - 3.*M_PI/2;
                }
                if(theta_c>7.*M_PI/4&&theta_c<=2.*M_PI) {
                    aux_theta = 2.*M_PI-theta_c;
                }

                REAL ratio1 = (REAL)iref/nrefs_c;
                REAL ratio2 = (REAL)(nrefs_c-iref)/nrefs_c;

                hipotenusa = (refradius/cos(aux_theta))*ratio1+refradius*ratio2;

            }else{

                hipotenusa = refradius;
            }

            coord = ParametricCircle(hipotenusa, theta_c);
            coord[0] += circle_pos;
            coord[1] += circle_pos;
            geomesh.NodeVec()[nodeindex].SetCoord(coord);
            geomesh.NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }


    int64_t elementid = 0;

    //father element
    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
    elementid++;

    for (int iref = 0; iref <= nrefs_c; ++iref) {
        for (int inode = 0; inode < n_ext_nds; inode++) {
            if(iref==nrefs_c){ //Pegar os nós externos
                TopolQuadrilateral[0] = n_ext_nds+n_cir_nds*nrefs_c+inode*2;
                TopolQuadrilateral[1] = inode;
                TopolQuadrilateral[2] = inode+1;
                if(TopolQuadrilateral[2]==n_ext_nds){
                    TopolQuadrilateral[2] = 0;
                }
                TopolQuadrilateral[3] = TopolQuadrilateral[0]+2;
                if(TopolQuadrilateral[3]==nodes){
                    TopolQuadrilateral[3] = n_ext_nds+n_cir_nds*nrefs_c;
                }

                TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_ext = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                son_ext->SetFather(father);
                son_ext->SetFatherIndex(father->Index());
                elementid++;

            }else{

                TopolQuadrilateral[0] = n_ext_nds+n_cir_nds*iref+inode*2;
                TopolQuadrilateral[1] = n_ext_nds+n_cir_nds*(iref+1)+inode*2;
                TopolQuadrilateral[2] = n_ext_nds+n_cir_nds*(iref+1)+inode*2+2;
                if(TopolQuadrilateral[2]==n_ext_nds+n_cir_nds*(iref+2)){
                    TopolQuadrilateral[2] = n_ext_nds+n_cir_nds*(iref+1);
                }
                TopolQuadrilateral[3] = n_ext_nds+n_cir_nds*iref+inode*2+2;
                if(TopolQuadrilateral[3]==n_ext_nds+n_cir_nds*(iref+1)){
                    TopolQuadrilateral[3] = n_ext_nds+n_cir_nds*(iref);
                }

                TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_int = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                son_int->SetFather(father);
                son_int->SetFatherIndex(father->Index());
                elementid++;
            }
        }
    }

//    TPZCheckGeom check(geomesh);
//    check.CheckUniqueId();

    geomesh.BuildConnectivity();
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&geomesh, out, true);

    TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(geomesh);
    TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
    if(!Found)
    {
        gRefDBase.InsertRefPattern(refPattern);
    }

    return refPattern;

}

TPZAutoPointer<TPZRefPattern> MHMNavierStokesTest::CreateGMeshVugs(int nrefs, TPZManVector<REAL,6> &FirstCoord, TPZManVector<REAL,6> &h_el)
{

    TPZGeoMesh geomesh;
    geomesh.SetDimension(2);

    int nrefs_c = 12; //number of refinements (radial direction)

    int nquadnods=pow(2.,1+nrefs)+1;
    int n_ext_nds = 8*pow(2,nrefs);
    int n_cir_nds = 16 * pow(2, nrefs);
    int n_vug_nds = 1;
    int n_vug_rings = nrefs+1;
    for (int k = 0; k < n_vug_rings ; ++k) {
        n_vug_nds += pow(2,k)*8;
    }
    int nodes_porous = n_ext_nds+n_cir_nds*(1+nrefs_c);
    int nodes = n_ext_nds+n_cir_nds*(1+nrefs_c)+n_vug_nds;


    REAL radius = h_el[0]/4.;
    geomesh.SetMaxNodeId(nodes-1);
    geomesh.NodeVec().Resize(nodes);
    //TPZManVector<TPZGeoNode,7> Node(nodes);

    TPZManVector<int64_t,8> TopolQuadrilateral(4);

    TPZManVector<REAL,3> coord(3,0.);

    TPZVec<REAL> xc(3,0.);

    REAL hx = h_el[0];
    REAL hy = h_el[1];

    int64_t nodeindex_r = 0,nodeindex = 0,  index_j1 = 0, index_j2 = 0;


    // External nodes for father
    TopolQuadrilateral[0] = n_ext_nds-nquadnods-(nquadnods-3)/2;
    TopolQuadrilateral[1] = TopolQuadrilateral[0]+nquadnods-1;
    TopolQuadrilateral[3] = n_ext_nds-1-2*(nquadnods-1)-(nquadnods-3)/2;
    TopolQuadrilateral[2] = TopolQuadrilateral[3]+1-nquadnods;

    //Exterior coordinates
    for(int i = 0; i < nquadnods; i++){
        for(int j = 0; j < nquadnods; j++){
            coord[0] = (j)*hx/(nquadnods - 1);
            coord[1] = (i)*hy/(nquadnods - 1);
            if((coord[0]>0&&coord[0]<hx)&&(coord[1]>0&&coord[1]<hy)){
                continue;
            }
            if(i==0){
                nodeindex_r=TopolQuadrilateral[0]+j;
            }else if(i>0&&i<nquadnods-1&&j==0){
                nodeindex_r=TopolQuadrilateral[0]-i;
            }else if(i>0&&i<nquadnods-1&&j==nquadnods-1&&coord[1]<hy/2.){
                index_j2++;
                nodeindex_r=TopolQuadrilateral[1]+index_j2;
            }else if(i>0&&i<nquadnods-1&&j==nquadnods-1&&coord[1]>hy/2.){
                index_j1++;
                nodeindex_r=index_j1;
            }else if(i==nquadnods-1){
                nodeindex_r=TopolQuadrilateral[3]-j;
            }else{
                nodeindex_r=0;
            }

            geomesh.NodeVec()[nodeindex_r].SetCoord(coord);
            geomesh.NodeVec()[nodeindex_r].SetNodeId(nodeindex);
        }
    }
    nodeindex = n_ext_nds;

    //Circunferences coordinates and nodes:
    REAL refradius = 0.;
    REAL circle_pos = hx/2.;
    REAL hipotenusa =0.;
    REAL theta_c = 0.,aux_theta =0.;
    //nrefs = 6;

    for (int iref = 0; iref <= nrefs_c; iref++) {
        for (int inode = 0; inode < n_cir_nds ; inode++) {
            // i node
            refradius = radius+iref*((circle_pos-radius)/(1+nrefs_c));
            theta_c = inode *2.*M_PI/n_cir_nds;
            if(iref>0) {
                aux_theta = theta_c;
                if(theta_c>M_PI/4.&&theta_c<=M_PI/2.) {
                    aux_theta = M_PI / 2. - theta_c;
                }
                if(theta_c>M_PI/2&&theta_c<=3.*M_PI/4.) {
                    aux_theta = theta_c - M_PI / 2.;
                }
                if(theta_c>3.*M_PI/4.&&theta_c<=M_PI) {
                    aux_theta = M_PI - theta_c;
                }
                if(theta_c>M_PI&&theta_c<=5.*M_PI/4.) {
                    aux_theta = theta_c - M_PI ;
                }
                if(theta_c>5.*M_PI/4.&&theta_c<=3.*M_PI/2.) {
                    aux_theta = 3.*M_PI/2.-theta_c;
                }
                if(theta_c>3.*M_PI/2.&&theta_c<=7.*M_PI/4) {
                    aux_theta = theta_c - 3.*M_PI/2;
                }
                if(theta_c>7.*M_PI/4&&theta_c<=2.*M_PI) {
                    aux_theta = 2.*M_PI-theta_c;
                }

                REAL ratio1 = (REAL)iref/nrefs_c;
                REAL ratio2 = (REAL)(nrefs_c-iref)/nrefs_c;

                hipotenusa = (refradius/cos(aux_theta))*ratio1+refradius*ratio2;

            }else{

                hipotenusa = refradius;
            }

            coord = ParametricCircle(hipotenusa, theta_c);
            coord[0] += circle_pos;
            coord[1] += circle_pos;
            geomesh.NodeVec()[nodeindex].SetCoord(coord);
            geomesh.NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }

    //Create vug nodes

    coord[0] = circle_pos;
    coord[1] = circle_pos;
    geomesh.NodeVec()[nodeindex].SetCoord(coord);
    geomesh.NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;


    for (int iref = 1; iref <= n_vug_rings; iref++) {
        int nnode_vug_int = pow(2,iref)*4;
        for (int inode = 0; inode < nnode_vug_int ; inode++) {
            // i node
            refradius = iref*(radius/(1+n_vug_rings));
            theta_c = inode *2.*M_PI/nnode_vug_int;

                aux_theta = theta_c;
                if(theta_c>M_PI/4.&&theta_c<=M_PI/2.) {
                    aux_theta = M_PI / 2. - theta_c;
                }
                if(theta_c>M_PI/2&&theta_c<=3.*M_PI/4.) {
                    aux_theta = theta_c - M_PI / 2.;
                }
                if(theta_c>3.*M_PI/4.&&theta_c<=M_PI) {
                    aux_theta = M_PI - theta_c;
                }
                if(theta_c>M_PI&&theta_c<=5.*M_PI/4.) {
                    aux_theta = theta_c - M_PI ;
                }
                if(theta_c>5.*M_PI/4.&&theta_c<=3.*M_PI/2.) {
                    aux_theta = 3.*M_PI/2.-theta_c;
                }
                if(theta_c>3.*M_PI/2.&&theta_c<=7.*M_PI/4) {
                    aux_theta = theta_c - 3.*M_PI/2;
                }
                if(theta_c>7.*M_PI/4&&theta_c<=2.*M_PI) {
                    aux_theta = 2.*M_PI-theta_c;
                }

                REAL ratio1 = (REAL)iref/n_vug_rings;
                REAL ratio2 = (REAL)(n_vug_rings-iref)/n_vug_rings;

                //hipotenusa = (refradius/cos(aux_theta))*ratio1+refradius*ratio2;
                hipotenusa = refradius;

            coord = ParametricCircle(hipotenusa, theta_c);
            coord[0] += circle_pos;
            coord[1] += circle_pos;
            geomesh.NodeVec()[nodeindex].SetCoord(coord);
            geomesh.NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }


    int64_t elementid = 0;

    //father element
    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
    elementid++;

    // Porous Media Elements:
    for (int iref = 0; iref <= nrefs_c; ++iref) {
        for (int inode = 0; inode < n_ext_nds; inode++) {
            if(iref==nrefs_c){ //Pegar os nós externos
                TopolQuadrilateral[0] = n_ext_nds+n_cir_nds*nrefs_c+inode*2;
                TopolQuadrilateral[1] = inode;
                TopolQuadrilateral[2] = inode+1;
                if(TopolQuadrilateral[2]==n_ext_nds){
                    TopolQuadrilateral[2] = 0;
                }
                TopolQuadrilateral[3] = TopolQuadrilateral[0]+2;
                if(TopolQuadrilateral[3]==nodes_porous){
                    TopolQuadrilateral[3] = n_ext_nds+n_cir_nds*nrefs_c;
                }

                TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_ext = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                son_ext->SetFather(father);
                son_ext->SetFatherIndex(father->Index());
                elementid++;

            }else{

                TopolQuadrilateral[0] = n_ext_nds+n_cir_nds*iref+inode*2;
                TopolQuadrilateral[1] = n_ext_nds+n_cir_nds*(iref+1)+inode*2;
                TopolQuadrilateral[2] = n_ext_nds+n_cir_nds*(iref+1)+inode*2+2;
                if(TopolQuadrilateral[2]==n_ext_nds+n_cir_nds*(iref+2)){
                    TopolQuadrilateral[2] = n_ext_nds+n_cir_nds*(iref+1);
                }
                TopolQuadrilateral[3] = n_ext_nds+n_cir_nds*iref+inode*2+2;
                if(TopolQuadrilateral[3]==n_ext_nds+n_cir_nds*(iref+1)){
                    TopolQuadrilateral[3] = n_ext_nds+n_cir_nds*(iref);
                }

                TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_int = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                son_int->SetFather(father);
                son_int->SetFatherIndex(father->Index());
                elementid++;
            }
        }
    }

    // Vug Elements (nucleus):



    if(0){

        for (int j = 0; j < 4; ++j) {
            TopolQuadrilateral[0] = nodes_porous;
            TopolQuadrilateral[1] = nodes_porous+1+2*j;
            TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
            TopolQuadrilateral[3] = TopolQuadrilateral[2]+1;
            if(TopolQuadrilateral[3]==nodes_porous+9){
                TopolQuadrilateral[3]=nodes_porous+1;
            }
            TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug_2 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
            son_vug_2->SetFather(father);
            son_vug_2->SetFatherIndex(father->Index());
            elementid++;
        }

        int node_pivot = nodes_porous;
        for (int iref = 1; iref <= n_vug_rings; ++iref) {
            int nnode_vug_int = pow(2,iref)*4;
            for (int inode = 0; inode < nnode_vug_int; inode++) {
                node_pivot++;

                if(iref==n_vug_rings){ //Pegar os nós externos
                    if(inode % 2 == 0){
                        TopolQuadrilateral[0] = node_pivot;
                        TopolQuadrilateral[1] = n_ext_nds+inode*4;
                        TopolQuadrilateral[2] = n_ext_nds+inode*4+2;
                        TopolQuadrilateral[3] = node_pivot+1;

                        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_ext = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                        son_ext->SetFather(father);
                        son_ext->SetFatherIndex(father->Index());
                        elementid++;

                    }else{

                        TopolQuadrilateral[0] = node_pivot;
                        TopolQuadrilateral[1] = n_ext_nds+(inode-1)*4+2;
                        TopolQuadrilateral[2] = TopolQuadrilateral[1]+2;
                        TopolQuadrilateral[3] = TopolQuadrilateral[2]+2;

                        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                        son_vug->SetFather(father);
                        son_vug->SetFatherIndex(father->Index());
                        elementid++;

                        TopolQuadrilateral[0] = node_pivot;
                        TopolQuadrilateral[1] = TopolQuadrilateral[3];
                        TopolQuadrilateral[2] = TopolQuadrilateral[1]+2;
                        TopolQuadrilateral[3] = node_pivot+1;

                        if(inode==nnode_vug_int-1){
                            TopolQuadrilateral[2] = TopolQuadrilateral[2]-n_cir_nds;
                            TopolQuadrilateral[3] = TopolQuadrilateral[3]-nnode_vug_int;
                        }

                        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug_2 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                        son_vug_2->SetFather(father);
                        son_vug_2->SetFatherIndex(father->Index());
                        elementid++;

                    }

                }else{

                    if(inode % 2 == 0){
                        TopolQuadrilateral[0] = node_pivot;
                        TopolQuadrilateral[1] = node_pivot+nnode_vug_int+inode;
                        TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
                        TopolQuadrilateral[3] = node_pivot+1;

                        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_ext = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                        son_ext->SetFather(father);
                        son_ext->SetFatherIndex(father->Index());
                        elementid++;

                    }else{

                        TopolQuadrilateral[0] = node_pivot;
                        TopolQuadrilateral[1] = node_pivot+nnode_vug_int+inode-1;
                        TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
                        TopolQuadrilateral[3] = TopolQuadrilateral[2]+1;

                        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                        son_vug->SetFather(father);
                        son_vug->SetFatherIndex(father->Index());
                        elementid++;

                        TopolQuadrilateral[0] = node_pivot;
                        TopolQuadrilateral[1] = TopolQuadrilateral[3];
                        TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
                        TopolQuadrilateral[3] = node_pivot+1;

                        if(inode==nnode_vug_int-1){
                            TopolQuadrilateral[2] = TopolQuadrilateral[2]-nnode_vug_int*2;
                            TopolQuadrilateral[3] = TopolQuadrilateral[3]-nnode_vug_int;
                        }

                        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug_2 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                        son_vug_2->SetFather(father);
                        son_vug_2->SetFatherIndex(father->Index());
                        elementid++;

                    }
                }
            }
        }
    }

    for (int j = 0; j < 4; ++j) {
        TopolQuadrilateral[0] = nodes_porous;
        TopolQuadrilateral[1] = nodes_porous+1+2*j;
        TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
        TopolQuadrilateral[3] = TopolQuadrilateral[2]+1;
        if(TopolQuadrilateral[3]==nodes_porous+9){
            TopolQuadrilateral[3]=nodes_porous+1;
        }
        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug_2 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
        son_vug_2->SetFather(father);
        son_vug_2->SetFatherIndex(father->Index());
        elementid++;
    }

    int node_pivot = nodes_porous;

    for (int iref = 1; iref <= nrefs; ++iref) {
        int nnode_vug_int = pow(2, iref) * 4;
        for (int inode = 0; inode < nnode_vug_int; inode++) {
            node_pivot++;
            if(inode % 2 == 0){
                TopolQuadrilateral[0] = node_pivot;
                TopolQuadrilateral[1] = node_pivot+nnode_vug_int+inode;
                TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
                TopolQuadrilateral[3] = node_pivot+1;

                TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_ext = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                son_ext->SetFather(father);
                son_ext->SetFatherIndex(father->Index());
                elementid++;

            }else{

                TopolQuadrilateral[0] = node_pivot;
                TopolQuadrilateral[1] = node_pivot+nnode_vug_int+inode-1;
                TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
                TopolQuadrilateral[3] = TopolQuadrilateral[2]+1;

                TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                son_vug->SetFather(father);
                son_vug->SetFatherIndex(father->Index());
                elementid++;

                TopolQuadrilateral[0] = node_pivot;
                TopolQuadrilateral[1] = TopolQuadrilateral[3];
                TopolQuadrilateral[2] = TopolQuadrilateral[1]+1;
                TopolQuadrilateral[3] = node_pivot+1;

                if(inode==nnode_vug_int-1){
                    TopolQuadrilateral[2] = TopolQuadrilateral[2]-nnode_vug_int*2;
                    TopolQuadrilateral[3] = TopolQuadrilateral[3]-nnode_vug_int;
                }

                TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug_2 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
                son_vug_2->SetFather(father);
                son_vug_2->SetFatherIndex(father->Index());
                elementid++;

            }
        }
    }

    for (int inode = 0; inode < n_ext_nds; inode++) {
        node_pivot++;
        TopolQuadrilateral[0] = node_pivot;
        TopolQuadrilateral[1] = n_ext_nds+2*inode;
        TopolQuadrilateral[2] = TopolQuadrilateral[1]+2;
        TopolQuadrilateral[3] = TopolQuadrilateral[0]+1;
        if(TopolQuadrilateral[2]==n_cir_nds+n_ext_nds){
            TopolQuadrilateral[2] = n_ext_nds;
            TopolQuadrilateral[3] = node_pivot-n_ext_nds+1;
        }
        TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son_vug_2 = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
        son_vug_2->SetFather(father);
        son_vug_2->SetFatherIndex(father->Index());
        elementid++;
    }




//    TPZCheckGeom check(geomesh);
//    check.CheckUniqueId();

    geomesh.BuildConnectivity();
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&geomesh, out, true);

    TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(geomesh);
    TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
    if(!Found)
    {
        gRefDBase.InsertRefPattern(refPattern);
    }

    return refPattern;

}

TPZAutoPointer<TPZRefPattern> MHMNavierStokesTest::CreateGMeshQuadRef(int nrefs, TPZManVector<REAL,6> &h_el)
{

    TPZGeoMesh geomesh;
    geomesh.SetDimension(2);

    int nnodes_x = pow(2.,1+nrefs)+1;
    int nnodes = nnodes_x*nnodes_x;
    geomesh.SetMaxNodeId(nnodes-1);
    geomesh.NodeVec().Resize(nnodes);
    TPZManVector<TPZGeoNode,7> Node(nnodes);

    TPZManVector<int64_t,8> TopolQuadrilateral(4);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);

    REAL hx = h_el[0];
    REAL hy = h_el[1];

    int64_t nodeindex = 0;

    //Exterior coordinates
    for(int i = 0; i < nnodes_x; i++){
        for(int j = 0; j < nnodes_x; j++){
            coord[0] = (j)*hx/(nnodes_x - 1);
            coord[1] = (i)*hy/(nnodes_x - 1);
            geomesh.NodeVec()[nodeindex].SetCoord(coord);
            geomesh.NodeVec()[nodeindex].SetNodeId(nodeindex);
            nodeindex++;
        }
    }

    int64_t elementid = 0;

    TopolQuadrilateral[0] = 0;
    TopolQuadrilateral[1] = nnodes_x-1;
    TopolQuadrilateral[2] = nnodes_x*nnodes_x-1;
    TopolQuadrilateral[3] = nnodes_x*nnodes_x-nnodes_x;

    TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * father = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
    elementid++;

    for(int i = 0; i < nnodes_x-1; i++){
        for(int j = 0; j < nnodes_x-1; j++){
            TopolQuadrilateral[0] = nnodes_x*i+j;
            TopolQuadrilateral[1] = TopolQuadrilateral[0]+1;
            TopolQuadrilateral[3] = nnodes_x*(i+1)+j;
            TopolQuadrilateral[2] = TopolQuadrilateral[3]+1;

            TPZGeoElRefPattern< pzgeom::TPZGeoQuad > * son = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,TopolQuadrilateral, fmatID,geomesh);
            son->SetFather(father);
            son->SetFatherIndex(father->Index());
            elementid++;
        }
    }

//    for(int iel  = 1; iel < nel; iel++){
//        geomesh.Element(iel)->SetFather(geomesh.Element(0));
//    }
    geomesh.BuildConnectivity();
    std::ofstream out("QuadGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&geomesh, out, true);

    TPZAutoPointer<TPZRefPattern> refPattern = new TPZRefPattern(geomesh);
    TPZAutoPointer<TPZRefPattern> Found = gRefDBase.FindRefPattern(refPattern);
    if(!Found)
    {
        gRefDBase.InsertRefPattern(refPattern);
    }

    return refPattern;

}




TPZGeoMesh *MHMNavierStokesTest::CreateGMeshCurve()
{
    
    TPZGeoMesh * geomesh = new TPZGeoMesh;
    geomesh->SetDimension(2);
    
    int nodes = 6;
    REAL radius = 1.0;
    REAL innerradius = radius/2.0;
    geomesh->SetMaxNodeId(nodes-1);
    geomesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,7> Node(nodes);
    
    TPZManVector<int64_t,6> TopolQQuadrilateral(6);
    TPZManVector<int64_t,8> TopolQuadrilateral(4);
    TPZManVector<int64_t,6> TopolQTriangle(6);
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<int64_t,3> TopolArc(3);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    
    int64_t nodeindex = 0;
    
    for (int inode = 0; inode < 3 ; inode++) {
        // i node
        coord = ParametricCircle(radius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    for (int inode = 0; inode < 3 ; inode++) {
        // i node
        coord = ParametricCircle(innerradius, inode * M_PI/4.0);
        geomesh->NodeVec()[nodeindex].SetCoord(coord);
        geomesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
        nodeindex++;
    }
    
    //Ponto 1
    int64_t elementid = 0;
    
    TopolQuadrilateral[0] = 3;
    TopolQuadrilateral[1] = 0;
    TopolQuadrilateral[2] = 2;
    TopolQuadrilateral[3] = 5;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad > > (elementid,TopolQuadrilateral, fmatID,*geomesh);
    elementid++;
    
    // outer arcs bc's
    
    TopolLine[0] = 3;
    TopolLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCbott,*geomesh);
    elementid++;
    
    TopolLine[0] = 2;
    TopolLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,TopolLine, fmatBCtop,*geomesh);
    elementid++;
    
    TopolArc[0] = 0;
    TopolArc[1] = 2;
    TopolArc[2] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCright,*geomesh);
    elementid++;
    
    TopolArc[0] = 5;
    TopolArc[1] = 3;
    TopolArc[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, fmatBCleft,*geomesh);
    elementid++;
    
    
    
    geomesh->BuildConnectivity();
    
    int nref = 0;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++) {
        int nel = geomesh->NElements();
        for (int iel = 0; iel < nel; iel++) {
            TPZGeoEl *gel = geomesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
    
    TPZCheckGeom check(geomesh);
    check.CheckUniqueId();
    geomesh->BuildConnectivity();
    
    
    std::ofstream out("CurvedGeometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(geomesh, out, true);
    
    return geomesh;
    
}

TPZManVector<REAL,3>  MHMNavierStokesTest::ParametricCircle(REAL radius,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = radius * cos(theta);
    xcoor[1] = radius * sin(theta);
    xcoor[2] = 0.0 ;
    return xcoor;
}

void MHMNavierStokesTest::TetrahedralMeshCubo(TPZVec<int> &n_s){
    
//    TPZGeoMesh *gmesh = new TPZGeoMesh;
//    GenerateNodes(gmesh,nelem);
//
//
//
//    return gmesh;
    
}


void MHMNavierStokesTest::SubdomainRefine(int nrefine, TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindices){
    
    int ncoarse = coarseindices.size();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int i_ref = 0; i_ref < nrefine; i_ref++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            if(!gel){
                 continue;
            }
            // BC elements
//            if(gel->MaterialId()<0){
//                    if(!gel->HasSubElement()) gel->Divide(filhos);
//            }
            
            if(gel->Dimension()!=gmesh->Dimension()){
                continue;
            }
            
            // Coarse volumetric elements
            TPZGeoEl * higher_el = gel->LowestFather();

            for (int i_coarse = 0; i_coarse<ncoarse; i_coarse++) {
                if(higher_el->Index()!=coarseindices[i_coarse]) continue;
            }
            
            if(!gel) continue;

            if(f_domaintype==TStokesAnalytic::ECavity){
                TPZVec<REAL> centerMaster(2,0.), centerEuclid(2,0.);;
                TPZGeoEl * higher_el = gel->LowestFather();
                int nsides = gel->NSides();
                int side = nsides - 1;
                higher_el->CenterPoint(side, centerMaster);
                higher_el->X(centerMaster,centerEuclid);

                if ((centerEuclid[0] < 0.1 && centerEuclid[1] > 0.9)||(centerEuclid[0] > 0.9 && centerEuclid[1] > 0.9)) {
                    if(!gel->HasSubElement()) gel->Divide(filhos);
                }
            }else{
                if(!gel->HasSubElement()) gel->Divide(filhos);
            }

            
        }
        
    }
    
//    int64_t nel = coarseindices.size();
//    for (int64_t el=0; el<nel; el++) {
//        TPZGeoEl *gel = gmesh->Element(coarseindices[el]);
//
//        int nsubel = gel->NSubElements();
//
//        for (int iref = 0; iref<nrefine; iref++)
//        {
//            TPZManVector<TPZGeoEl *,10> gelsub;
//            gel->Divide(gelsub);
//        }
//    }
    
}




void MHMNavierStokesTest::UniformRefine4(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            higher_el->CenterPoint(side, centerMaster);
            higher_el->X(centerMaster,centerEuclid);
            
            
            if (fabs(centerCo[0]-centerEuclid[0]) > 1.e-9 &&  restriction == true) {
                continue;
            }
            if (fabs(centerCo[1]-centerEuclid[1]) > 1.e-9 && restriction == true) {
                continue;
            }
            
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        break;
             //           neighbour.Element()->Divide(filhos);
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}


void MHMNavierStokesTest::UniformRefine3(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        int count =0.;
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            if(higher_el->Dimension()!=2) continue;
            
            int intdiv = (higher_el->Index()/(n_div[0]*2))%2;
            if (intdiv==0) {
                count=0;
            }else{
                count=1;
            }

            count =0; //papapapa
            if((higher_el->Index()+count)%2!=0) continue;
            
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        
                        if(f_allrefine==false){
                            break;
                        }
                        neighbour.Element()->Divide(filhos);
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}


void MHMNavierStokesTest::UniformRefine2(int nDiv, TPZGeoMesh *gmesh, TPZVec<int> &n_div)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        int count =0.;
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            if(higher_el->Dimension()!=2) continue;
            
            int intdiv = (higher_el->Index()/n_div[0])%2;
            if (intdiv==0) {
                count=0;
            }else{
                count=1;
            }
            
            //if((higher_el->Index()+count)%2==0) continue;
          
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        
                        if (f_allrefine==false) {
                            break;
                        }
                        
                        neighbour.Element()->Divide(filhos);
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}



void MHMNavierStokesTest::UniformRefine(int nDiv, TPZGeoMesh *gmesh, TPZVec<REAL> centerCo, bool restriction)
{
    
    int dim = gmesh->Dimension();
    TPZManVector< TPZGeoEl *,20 > filhos;
    for(int D = 0; D < nDiv; D++)
    {
        TPZAdmChunkVector<TPZGeoEl *> gelvec = gmesh->ElementVec();
        int nels = gmesh->NElements();
        
        for(int elem = 0; elem < nels; elem++)
        {
            
            TPZGeoEl * gel = gelvec[elem];
            
            TPZGeoEl * higher_el = gel->LowestFather();
            TPZVec<REAL> centerMaster(3,0.), centerEuclid(3,0.);;
            int nsides = gel->NSides();
            int side = nsides - 1;
            higher_el->CenterPoint(side, centerMaster);
            higher_el->X(centerMaster,centerEuclid);
            
            if (fabs(centerCo[0]-centerEuclid[0]) > 1.e-9 &&  restriction == true) {
                continue;
            }
            if (fabs(centerCo[1]-centerEuclid[1]) > 1.e-9 && restriction == true) {
                continue;
            }
            
            unsigned int n_corner_sides = gel->NCornerNodes();
            
            for (int i_s=n_corner_sides; i_s<nsides; i_s++) {
                TPZGeoElSide gelside(gel,i_s);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while(neighbour != gelside)
                {
                    if(neighbour.Element()->Dimension()== dim-1){
                        if (f_allrefine==false) {
                            break;
                        }
                        neighbour.Element()->Divide(filhos);
                        
                    }
                    neighbour = neighbour.Neighbour();
                }
            }
            
            if(!gel) continue;
            if(!gel->HasSubElement()) gel->Divide(filhos);
            
        }
    }
    
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

TPZCompEl *MHMNavierStokesTest::CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index) {
    if(!gel->Reference() && gel->NumInterfaces() == 0)
        return new TPZInterfaceElement(mesh,gel,index);
    
    return NULL;
}

void MHMNavierStokesTest::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
    
//        dsol.Resize(3,3);
//        sol.Resize(4);
//
//        REAL x1 = x[0];
//        REAL x2 = x[1];
//
//        TPZVec<REAL> v_Dirichlet(3,0.);
//
////        v_Dirichlet[0] = -0.1*x2*x2+0.2*x2;
//        v_Dirichlet[0] = -1.+x2;
////        v_Dirichlet[0] = 1.;
//        v_Dirichlet[1] = 0.;
//        v_Dirichlet[2] = 0.;
//
////        STATE pressure = 1.-0.2*x1;
//        STATE pressure = 0.;
//
//        sol[0]=v_Dirichlet[0];
//        sol[1]=v_Dirichlet[1];
//        sol[2]=v_Dirichlet[2];
//        sol[3]=pressure;
//
//        // vx direction
//        dsol(0,0)= 0.;
////        dsol(0,1)= 0.2-0.2*x2;
//        dsol(0,1)= 1.;
//        //dsol(0,1)= 0.;
//        dsol(0,2)= 0.;
//
//        // vy direction
//        dsol(1,0)= 0.;
//        dsol(1,1)= 0.;
//        dsol(1,2)= 0.;
//
//        // vz direction
//        dsol(2,0)= 0.;
//        dsol(2,1)= 0.;
//        dsol(2,2)= 0.;
    
    // General form : : Artigo Botti, Di Pietro, Droniou
    
    //    dsol.Resize(3,3);
    //    sol.Resize(3);
    //
    //    REAL x1 = x[0];
    //    REAL x2 = x[1];
    //
    //    REAL m_v= 1., m_u= 1.0;
    //
    //    REAL Cf=m_v/m_u;
    //
    //    STATE v_1 = -exp(-Cf)*sin(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*sin(x1)*sin(x2);
    //    STATE v_2 = -exp(-Cf)*cos(x1)*cos(x2)-(1./m_v)*(1.-exp(-Cf))*cos(x1)*cos(x2);
    //    STATE pressure= cos(x1)*sin(x2);
    //
    //    sol[0]=v_1;
    //    sol[1]=v_2;
    //    sol[2]=pressure;
    //
    //    // vx direction
    //    dsol(0,0)= -exp(-Cf)*cos(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
    //    dsol(0,1)= exp(-Cf)*cos(x2)*sin(x1)+(1./m_v)*(1.-exp(-Cf))*cos(x2)*sin(x1);
    //
    //    // vy direction
    //    dsol(1,0)= -exp(-Cf)*cos(x2)*sin(x1)+(1./m_v)*(1.-exp(-Cf))*cos(x2)*sin(x1);
    //    dsol(1,1)= exp(-Cf)*cos(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
    //

    
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou
    
    //    dsol.Resize(3,3);
    //    sol.Resize(3);
    //
    //    REAL x1 = x[0];
    //    REAL x2 = x[1];
    //
    //    REAL e = exp(1.);
    //
    //    STATE v_1 = (1.-2./e)*sin(x1)*sin(x2);
    //    STATE v_2 = -1.*cos(x1)*cos(x2);
    //    STATE pressure= cos(x1)*sin(x2);
    //
    //    sol[0]=v_1;
    //    sol[1]=v_2;
    //    sol[2]=pressure;
    //
    //    // vx direction
    //    dsol(0,0)= (1.-2./e)*cos(x1)*sin(x2);
    //    dsol(0,1)= cos(x2)*sin(x1);
    //
    //    // vy direction
    //    dsol(1,0)= (1.-2./e)*cos(x2)*sin(x1);
    //    dsol(1,1)= cos(x1)*sin(x2);
    //

    
    // Stokes : : Artigo Botti, Di Pietro, Droniou
    
//    dsol.Resize(3,3);
//    sol.Resize(4);
//
//
//    //Applying rotation:
//    TPZVec<REAL> x_in = x;
//    TPZVec<REAL> x_rot(3,0.);
//
//    f_InvT.Apply(x_in,x_rot);
//    x[0] = x_rot[0];
//    x[1] = x_rot[1];
//
//    REAL x1 = x[0];
//    REAL x2 = x[1];
//
//    REAL e = exp(1.);
//
//    TPZVec<REAL> v_Dirichlet(3,0.), vbc_rot(3,0.);
//
//    v_Dirichlet[0] = -1.*sin(x1)*sin(x2);
//    v_Dirichlet[1] = -1.*cos(x1)*cos(x2);
//    STATE pressure= cos(x1)*sin(x2);
//
//    f_T.Apply(v_Dirichlet, vbc_rot);
//    v_Dirichlet = vbc_rot;
//
//    sol[0]=v_Dirichlet[0];
//    sol[1]=v_Dirichlet[1];
//    sol[2]=v_Dirichlet[2];
//    sol[3]=pressure;
//
//
//    // GradU * Rt
//    TPZFMatrix<STATE> GradU(3,3,0.), GradURt(3,3,0.), RGradURt(3,3,0.);
//
//    // vx direction
//    GradU(0,0)= -1.*cos(x1)*sin(x2);
//    GradU(0,1)= cos(x2)*sin(x1);
//
//    // vy direction
//    GradU(1,0)= -1.*cos(x2)*sin(x1);
//    GradU(1,1)= cos(x1)*sin(x2);
//
//    TPZFMatrix<STATE> R = f_T.Mult();
//    TPZFMatrix<STATE> Rt(3,3,0.);
//    R.Transpose(&Rt);
//
////    GradU.Print("GradU = ");
////    R.Print("R = ");
////    Rt.Print("Rt = ");
//
//    GradU.Multiply(Rt,GradURt);
////    GradURt.Print("GradURt = ");
//
//    R.Multiply(GradURt,RGradURt);
////    RGradURt.Print("RGradURt = ");
//
//    // vx direction
//    dsol(0,0)= RGradURt(0,0);
//    dsol(0,1)= RGradURt(0,1);
//    dsol(0,2)= RGradURt(0,2);
//
//    // vy direction
//    dsol(1,0)= RGradURt(1,0);
//    dsol(1,1)= RGradURt(1,1);
//    dsol(1,2)= RGradURt(1,2);
//
//    // vz direction
//    dsol(2,0)= RGradURt(2,0);
//    dsol(2,1)= RGradURt(2,1);
//    dsol(2,2)= RGradURt(2,2);
    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
    //        dsol.Resize(3,3);
    //        sol.Resize(3);
    //
    //        REAL x1 = x[0];
    //        REAL x2 = x[1];
    //
    //        STATE v_1 = sin(x1)*sin(x2);
    //        STATE v_2 = -1.*cos(x1)*cos(x2);
    //        STATE pressure= cos(x1)*sin(x2);
    //
    //        sol[0]=v_1;
    //        sol[1]=v_2;
    //        sol[2]=pressure;
    //
    //        // vx direction
    //        dsol(0,0)= cos(x1)*sin(x2);
    //        dsol(0,1)= cos(x2)*sin(x1);
    //
    //        // vy direction
    //        dsol(1,0)= cos(x2)*sin(x1);
    //        dsol(1,1)= cos(x1)*sin(x2);
    

    
    // Stokes 3D : Artigo Botti, Di Pietro, Droniou
    
//        dsol.Resize(3,3);
//        sol.Resize(4);
//
//        //Applying rotation:
//        TPZVec<REAL> x_in = x;
//        TPZVec<REAL> x_rot(3,0.);
//
//        f_InvT.Apply(x_in,x_rot);
//        x[0] = x_rot[0];
//        x[1] = x_rot[1];
//
//        REAL x1 = x[0];
//        REAL x2 = x[1];
//        REAL x3 = x[2];
//
//        TPZVec<REAL> v_Dirichlet(3,0.), vbc_rot(3,0.);
//
//        v_Dirichlet[0] = cos(x1)*cos(x3) -1.*sin(x1)*sin(x2);
//        v_Dirichlet[1] = -1.*cos(x1)*cos(x2);
//        v_Dirichlet[2] = sin(x1)*sin(x3);
//        STATE pressure= cos(x1)*sin(x2)*cos(x3);
//
//        f_T.Apply(v_Dirichlet, vbc_rot);
//        v_Dirichlet = vbc_rot;
//
//        sol[0]=v_Dirichlet[0];
//        sol[1]=v_Dirichlet[1];
//        sol[2]=v_Dirichlet[2];
//        sol[3]=pressure;
//
//
//        // GradU * Rt
//        TPZFMatrix<STATE> GradU(3,3,0.), GradURt(3,3,0.), RGradURt(3,3,0.);
//
//        // vx direction
//        GradU(0,0)= -cos(x3)*sin(x1)-1.*cos(x1)*sin(x2);
//        GradU(0,1)= cos(x2)*sin(x1);
//        GradU(0,2)= cos(x1)*sin(x3);
//
//        // vy direction
//        GradU(1,0)= -1.*cos(x2)*sin(x1);
//        GradU(1,1)= cos(x1)*sin(x2);
//        GradU(1,2)= 0.;
//
//        // vz direction
//        GradU(2,0)= -1.*cos(x1)*sin(x3);
//        GradU(2,1)= 0.;
//        GradU(2,2)= cos(x3)*sin(x1);
//
//
//        TPZFMatrix<STATE> R = f_T.Mult();
//        TPZFMatrix<STATE> Rt(3,3,0.);
//        R.Transpose(&Rt);
//
//    //    GradU.Print("GradU = ");
//    //    R.Print("R = ");
//    //    Rt.Print("Rt = ");
//
//        GradU.Multiply(Rt,GradURt);
//    //    GradURt.Print("GradURt = ");
//
//        R.Multiply(GradURt,RGradURt);
//    //    RGradURt.Print("RGradURt = ");
//
//        // vx direction
//        dsol(0,0)= RGradURt(0,0);
//        dsol(0,1)= RGradURt(0,1);
//        dsol(0,2)= RGradURt(0,2);
//
//        // vy direction
//        dsol(1,0)= RGradURt(1,0);
//        dsol(1,1)= RGradURt(1,1);
//        dsol(1,2)= RGradURt(1,2);
//
//        // vz direction
//        dsol(2,0)= RGradURt(2,0);
//        dsol(2,1)= RGradURt(2,1);
//        dsol(2,2)= RGradURt(2,2);

        dsol.Resize(3,3);
        sol.Resize(4);

        REAL x1 = x[0];
        REAL x2 = x[1];
        REAL x3 = x[2];

        STATE v_1 = 0.;
        if(x2*x2+x3*x3<0.5*0.5){
            v_1 = 10.;
        }
        STATE v_2 = 0.;

        sol[0]=v_1;
        sol[1]=0.;
        sol[2]=0.;
        sol[3]=0.;


    
    
}

void MHMNavierStokesTest::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu){
    
    //Applying rotation:
    TPZVec<REAL> x_in = x;
    TPZVec<REAL> x_rot(3,0.);
    
    f_InvT.Apply(x_in,x_rot);
    x[0] = x_rot[0];
    x[1] = x_rot[1];
    
    f.resize(3);
    REAL x1 = x[0];
    REAL x2 = x[1];
    REAL x3 = x[2];
    
    f[0] =0.;
    f[1] =0.;
    f[2] =0.;
    
    TPZVec<REAL> f_s(3,0), f_rot(3,0);
    
    // General form : : Artigo Botti, Di Pietro, Droniou
    
    //    REAL m_v= 1., m_u= 1.0;
    //
    //    REAL Cf=m_v/m_u;
    //
    //        f_1 = -sin(x1)*sin(x2)-exp(-Cf)*sin(x1)*sin(x2)+(1./m_v)*(1.-exp(-Cf))*sin(x1)*sin(x2)-m_u*(2.*exp(-Cf)*sin(x1)*sin(x2)-(1./m_v)*4.*(1-exp(-Cf))*sin(x1)*sin(x2));
    //
    //        f_2 = cos(x1)*cos(x2)-exp(-Cf)*cos(x1)*cos(x2)-(1./m_v)*(1.-exp(-Cf))*cos(x1)*cos(x2)-m_u*(2.*exp(-Cf)*cos(x1)*cos(x2)+(1./m_v)*4.*(1-exp(-Cf))*cos(x1)*cos(x2));
    //        STATE g_1 = (2./m_v)*(1.-exp(-Cf))*cos(x1)*sin(x2);
    //
    //        f[0] = f_1; // x direction
    //        f[1] = f_2; // y direction
    //
    //        f[2] = g_1; // g source
    
    
    // Brinkman : : Artigo Botti, Di Pietro, Droniou
    
    //    REAL e = exp(1.);
    //
    //    f_1 = (-8./e+ 4.)*sin(x1)*sin(x2);
    //    f_2 = (2./e- 4.)*cos(x1)*cos(x2);
    //    STATE g_1 = 2.*(1.-1./e)*cos(x1)*sin(x2);
    //
    //    f[0] = f_1; // x direction
    //    f[1] = f_2; // y direction
    //
    //    f[2] = g_1; // g source
    
    // Stokes : : Artigo Botti, Di Pietro, Droniou
    
    
    f_s[0] = -3.*sin(x1)*sin(x2);
    f_s[1] = -1.*cos(x1)*cos(x2);

    f_T.Apply(f_s, f_rot);
    f_s = f_rot;


    f[0] = f_s[0]; // x direction
    f[1] = f_s[1]; // y direction
    f[2] = f_s[2];
    
    
    // Darcy : : Artigo Botti, Di Pietro, Droniou
    
    //        f_1 = 0.;
    //        f_2 = 0.;
    //
    //        f[0] = f_1; // x direction
    //        f[1] = f_2; // y direction
    //        f[2] = 2.*cos(x1)*sin(x2);
    
    
    // Stokes 3D : Artigo Botti, Di Pietro, Droniou
    
    
//        f_s[0] = 2.*cos(x1)*cos(x3) - 1.*(2. + cos(x3))*sin(x1)*sin(x2);
//        f_s[1] = cos(x1)*cos(x2)*(-2. + cos(x3));
//        f_s[2] = (2.*sin(x1) - cos(x1)*sin(x2))*sin(x3);
//
//        f_T.Apply(f_s, f_rot);
//        f_s = f_rot;
//
//
//        f[0] = f_s[0]; // x direction
//        f[1] = f_s[1]; // y direction
//        f[2] = f_s[2];
    
}

void MHMNavierStokesTest::ChangeExternalOrderConnects(TPZCompMesh *mesh, int addToOrder){
    
    int nEl= mesh-> NElements();
    int dim = mesh->Dimension();
    
    for (int iel=0; iel<nEl; iel++) {
        TPZCompEl *cel = mesh->ElementVec()[iel];
        if (!cel) continue;
        int ncon = cel->NConnects();
        int corder = 0;
        int nshape = 0;
        int nshape2 = 0;
        
        if(cel->Dimension()== dim)
        {
            TPZConnect &conel = cel->Connect(ncon-1);
            corder = conel.Order();
            nshape = conel.NShape();
            
            int neworder = corder + addToOrder;//Aqui = +1
            int64_t cindex = cel->ConnectIndex(ncon-1);
            conel.SetOrder(neworder,cindex);
            
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            intel->SetPreferredOrder(neworder);
            nshape = intel->NConnectShapeF(ncon-1,neworder);
            
            if(dim==2 && addToOrder==1)
            {
                if(feltype==ETriangle){
                    nshape2 = (corder + 2)*(corder + 2)-1;
                }else{//Quadrilateral
                    nshape2 = 2*(corder + 1)*(corder + 2);
                }
                if(nshape2!=nshape)
                {
                    DebugStop();
                }
            }
            
            conel.SetNShape(nshape);
            mesh->Block().Set(conel.SequenceNumber(),nshape);
        }
    }
    mesh->CleanUpUnconnectedNodes();
    mesh->ExpandSolution();
}

void MHMNavierStokesTest::InsertMaterialObjects(TPZMHMeshControl *control)
{

    TPZCompMesh &cmesh = control->CMesh();

    TPZMHMNavierStokesMaterial *mat1 = new TPZMHMNavierStokesMaterial(fmatID,fdim);
    mat1->SetSimulationData(f_sim_data);

    TPZMHMNavierStokesMaterial *mat2 = NULL;
    if(f_domaintype==TStokesAnalytic::ECouplingSD||f_domaintype==TStokesAnalytic::ECouplingNSD||f_domaintype==TStokesAnalytic::EVugs2D){
        //Insert a Darcy material domain;
        mat2 = new TPZMHMNavierStokesMaterial(2,fdim);
        mat2->SetSimulationData(f_sim_data);
        mat2->SetViscosity(0.);
        mat2->SetBrinkman(1./f_sim_data->GetPermeability());
    }

    TPZAutoPointer<TPZFunction<STATE> > fp = f_ExactSol.ForcingFunction();
    TPZAutoPointer<TPZFunction<STATE> > solp = f_ExactSol.Exact();

    if(f_domaintype==TStokesAnalytic::ECavity||f_domaintype==TStokesAnalytic::EObstacles||f_domaintype==TStokesAnalytic::EVugs2D) {

    }else{
        mat1->SetForcingFunction(fp); //Caso simples sem termo fonte
        mat1->SetForcingFunctionExact(solp);
    }

    if (f_sim_data->GetShapeTest()) {
        mat1->SetForcingFunction(NULL);
        mat1->SetForcingFunctionExact(NULL);
    }

    //TPZMaterial * mat1(material);
    cmesh.InsertMaterialObject(mat1);

    if(f_domaintype==TStokesAnalytic::ECouplingSD||f_domaintype==TStokesAnalytic::ECouplingNSD){
        TPZAutoPointer<TPZFunction<STATE> > fp2 = f_ExactSol_2.ForcingFunction();
        f_ExactSol_2.fcBrinkman=1.;
        f_ExactSol_2.fvisco=0.;
        f_ExactSol_2.fExactSol = f_domaintype;
        mat2->SetProblemType(f_ExactSol_2.fProblemType); //Material 2 -> always EBrinkman (Darcy equation)
        mat2->SetForcingFunction(fp2);
        mat2->SetForcingFunctionExact(solp);
        cmesh.InsertMaterialObject(mat2);
    }
    if(f_domaintype==TStokesAnalytic::EVugs2D){
        f_ExactSol_2.fcBrinkman=1.;
        f_ExactSol_2.fvisco=0.;
        f_ExactSol_2.fExactSol = f_domaintype;
        mat2->SetProblemType(f_ExactSol_2.fProblemType); //Material 2 -> always EBrinkman (Darcy equation)

        solp = new TPZDummyFunction<STATE> (Sol_exact, 7); //oioioio
        mat2->SetForcingFunctionExact(solp); //oioioioi

        cmesh.InsertMaterialObject(mat2);
    }

    int matSkeleton = 4;

    control->SwitchLagrangeMultiplierSign(false);
    ///Inserir condicao de contorno
    TPZFMatrix<STATE> val1(3,3,0.), val2(3,1,0.);

    switch(f_domaintype) {

        case TStokesAnalytic::ECavity: //Pressure
        {
            val2(0, 0) = 0.0;
            TPZBndCond *BC_bott = mat1->CreateBC(mat1, fmatBCbott, fdirichlet_v, val1, val2);
            cmesh.InsertMaterialObject(BC_bott);

            val2(0, 0) = 0.0; // vx -> 1
            TPZBndCond *BC_top = mat1->CreateBC(mat1, fmatBCtop, fdirichlet_v, val1, val2);
            cmesh.InsertMaterialObject(BC_top);

            val2(0, 0) = 0.0;
            TPZBndCond *BC_left = mat1->CreateBC(mat1, fmatBCleft, fdirichlet_v, val1, val2);
            cmesh.InsertMaterialObject(BC_left);

            val2(0, 0) = 0.0;
            TPZBndCond *BC_right = mat1->CreateBC(mat1, fmatBCright, fdirichlet_v, val1, val2);
            cmesh.InsertMaterialObject(BC_right);

        }
            break;

        case TStokesAnalytic::EVugs2D: //Pressure
        {
            val2(0, 0) = 0.0;
            TPZBndCond *BC_bott = mat2->CreateBC(mat2, fmatBCbott, fneumann_v, val1, val2);
            cmesh.InsertMaterialObject(BC_bott);

            val2(0, 0) = 0.0; // vx -> 1
            TPZBndCond *BC_top = mat2->CreateBC(mat2, fmatBCtop, fneumann_v, val1, val2);
            cmesh.InsertMaterialObject(BC_top);

            val2(0, 0) = 1.0;
            TPZBndCond *BC_left = mat2->CreateBC(mat2, fmatBCleft, fdirichlet_v, val1, val2);
            BC_left->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_left);

            val2(0, 0) = 0.0;
            TPZBndCond *BC_right = mat2->CreateBC(mat2, fmatBCright, fneumann_v, val1, val2);
            cmesh.InsertMaterialObject(BC_right);

            TPZBndCond *BC_BJS = mat1->CreateBC(mat1, fmatBChole, f_BJS_condition, val1, val2);
            cmesh.InsertMaterialObject(BC_BJS);

            if (f_3Dmesh) {
                TPZBndCond * BCondD5 = mat2->CreateBC(mat2, fmatBCtop_z, fneumann_v, val1, val2);
                cmesh.InsertMaterialObject(BCondD5);

                TPZBndCond * BCondD6 = mat2->CreateBC(mat2, fmatBCbott_z, fneumann_v, val1, val2);
                cmesh.InsertMaterialObject(BCondD6);
            }

        }
            break;

        case TStokesAnalytic::EObstacles:
        {
            val2(0, 0) = 0.0;
            TPZBndCond *BC_bott = mat1->CreateBC(mat1, fmatBCbott, fdirichlet_v, val1, val2);
            cmesh.InsertMaterialObject(BC_bott);

            val2(0, 0) = 0.0; // vx -> 1
            TPZBndCond *BC_top = mat1->CreateBC(mat1, fmatBCtop, fdirichlet_v, val1, val2);
            cmesh.InsertMaterialObject(BC_top);

            val2(0, 0) = 1.0;
            TPZBndCond *BC_left = mat1->CreateBC(mat1, fmatBCleft, fdirichlet_v, val1, val2);
            BC_left->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_left);

            val2(0, 0) = 0.0;
            TPZBndCond *BC_right = mat1->CreateBC(mat1, fmatBCright, fneumann_v, val1, val2);
            cmesh.InsertMaterialObject(BC_right);

            val2(0,0) = 0.0;
            TPZBndCond * BC_hole = mat1->CreateBC(mat1, fmatBChole, fdirichlet_v, val1, val2);
            cmesh.InsertMaterialObject(BC_hole);

        }
            break;

        case TStokesAnalytic::ECouplingSD:
        case TStokesAnalytic::ECouplingNSD:
        {

            TPZBndCond *BC_bott = mat2->CreateBC(mat2, fmatBCbott, fdirichlet_v, val1, val2);
            BC_bott->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_bott);

            TPZBndCond *BC_top = mat1->CreateBC(mat1, fmatBCtop, fneumann_v, val1, val2);
            BC_top->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_top);

            TPZBndCond *BC_left = mat1->CreateBC(mat1, fmatBCleft, fneumann_v, val1, val2);
            BC_left->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_left);

            TPZBndCond *BC_left_2 = mat2->CreateBC(mat2, fmatBCleft_2, fneumann_v, val1, val2);
            BC_left_2->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_left_2);

            TPZBndCond *BC_right = mat1->CreateBC(mat1, fmatBCright, fneumann_v, val1, val2);
            BC_right->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_right);

            TPZBndCond *BC_right_2 = mat2->CreateBC(mat2, fmatBCright_2, fneumann_v, val1, val2);
            BC_right_2->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BC_right_2);
        }
            break;

        default: {

            TPZBndCond * BCondD1 = mat1->CreateBC(mat1, fmatBCbott, fneumann_v, val1, val2);
            BCondD1->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BCondD1);

            val1.Zero();
            TPZBndCond * BCondD2 = mat1->CreateBC(mat1, fmatBCtop, fneumann_v, val1, val2);
            BCondD2->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BCondD2);

            val1.Zero();
            val2.Zero();
            TPZBndCond * BCondD3 = mat1->CreateBC(mat1, fmatBCleft, fneumann_v, val1, val2);
            BCondD3->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BCondD3);

            TPZBndCond * BCondD4 = mat1->CreateBC(mat1, fmatBCright, fneumann_v, val1, val2);
            BCondD4->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(BCondD4);
//qwqwqwqw
            if (f_3Dmesh) {
                TPZBndCond * BCondD5 = mat1->CreateBC(mat1, fmatBCtop_z, fdirichlet_v, val1, val2);
                BCondD5->SetBCForcingFunction(0, solp);
                cmesh.InsertMaterialObject(BCondD5);

                TPZBndCond * BCondD6 = mat1->CreateBC(mat1, fmatBCbott_z, fdirichlet_v, val1, val2);
                BCondD6->SetBCForcingFunction(0, solp);
                cmesh.InsertMaterialObject(BCondD6);
            }

        }
    }

    //Skeleton::

    TPZBndCond * bcFlux = mat1->CreateBC(mat1, matSkeleton, fneumann_v, val1, val2);
    //bcFlux->SetBCForcingFunction(0, solp);
    cmesh.InsertMaterialObject(bcFlux);
    
    // 2.1 - Material para tração tangencial 1D (Interior)
    TPZNullMaterial *matLambda = new TPZNullMaterial(fmatLambda);
    matLambda->SetDimension(fdim-1);
    matLambda->SetNStateVariables(fdim-1);
    control->CMesh()->InsertMaterialObject(matLambda);

//    if(f_domaintype==TStokesAnalytic::EVugs2D){
//        TPZNullMaterial *matBJS = new TPZNullMaterial(fmatBChole);
//        matBJS->SetDimension(fdim-1);
//        matBJS->SetNStateVariables(fdim-1);
//        control->CMesh()->InsertMaterialObject(matBJS);
//    }

    // 2.2 - Material for interfaces (Interior)
    TPZMHMNavierStokesMaterial *matInterfaceLeft = new TPZMHMNavierStokesMaterial(control->fLagrangeMatIdLeft,fdim);
    matInterfaceLeft->SetSimulationData(f_sim_data);
    matInterfaceLeft->SetMultiplier(1.);
    cmesh.InsertMaterialObject(matInterfaceLeft);
    
    TPZMHMNavierStokesMaterial *matInterfaceRight = new TPZMHMNavierStokesMaterial(control->fLagrangeMatIdRight,fdim);
    matInterfaceRight->SetSimulationData(f_sim_data);
    matInterfaceRight->SetMultiplier(-1.);
    cmesh.InsertMaterialObject(matInterfaceRight);
    
    // 3.1 - Material para tração tangencial 1D nos contornos

    switch(f_domaintype) {

        case TStokesAnalytic::ECavity: {

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_bott = mat1->CreateBC(mat1, fmatLambdaBC_bott, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_bott);

            val2(0,0) = 1.;
            TPZBndCond *matLambdaBC_top = mat1->CreateBC(mat1, fmatLambdaBC_top, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_top);

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_left = mat1->CreateBC(mat1, fmatLambdaBC_left, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_left);

            TPZBndCond *matLambdaBC_right = mat1->CreateBC(mat1, fmatLambdaBC_right, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_right);

        }
            break;


        case TStokesAnalytic::EVugs2D: {

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_bott = mat2->CreateBC(mat2, fmatLambdaBC_bott, fneumann_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_bott);

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_top = mat2->CreateBC(mat2, fmatLambdaBC_top, fneumann_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_top);

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_left = mat2->CreateBC(mat2, fmatLambdaBC_left, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_left);

            TPZBndCond *matLambdaBC_right = mat2->CreateBC(mat2, fmatLambdaBC_right, fneumann_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_right);

            if (f_3Dmesh) {
                TPZBndCond *matLambdaBC_top_z = mat2->CreateBC(mat2, fmatLambdaBC_top_z, fneumann_sigma, val1, val2);
                cmesh.InsertMaterialObject(matLambdaBC_top_z);

                TPZBndCond *matLambdaBC_bott_z = mat2->CreateBC(mat2, fmatLambdaBC_bott_z, fneumann_sigma, val1, val2);
                cmesh.InsertMaterialObject(matLambdaBC_bott_z);
            }


        }
            break;

        case TStokesAnalytic::EObstacles: {

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_bott = mat1->CreateBC(mat1, fmatLambdaBC_bott, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_bott);

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_top = mat1->CreateBC(mat1, fmatLambdaBC_top, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_top);

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_left = mat1->CreateBC(mat1, fmatLambdaBC_left, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_left);

            TPZBndCond *matLambdaBC_right = mat1->CreateBC(mat1, fmatLambdaBC_right, fneumann_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_right);

            val2(0,0) = 0.0;
            TPZBndCond *matLambdaBC_hole = mat1->CreateBC(mat1, fmatLambdaBC_hole, fdirichlet_sigma, val1, val2);
            cmesh.InsertMaterialObject(matLambdaBC_hole);

        }
            break;

        case TStokesAnalytic::ECouplingSD:
        case TStokesAnalytic::ECouplingNSD:
        {

            TPZBndCond *matLambdaBC_bott = mat2->CreateBC(mat2, fmatLambdaBC_bott, fdirichlet_sigma, val1, val2);
            matLambdaBC_bott->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_bott);

            TPZBndCond *matLambdaBC_top = mat1->CreateBC(mat1, fmatLambdaBC_top, fdirichlet_sigma, val1, val2);
            matLambdaBC_top->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_top);

            TPZBndCond *matLambdaBC_left = mat1->CreateBC(mat1, fmatLambdaBC_left, fdirichlet_sigma, val1, val2);
            matLambdaBC_left->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_left);

            TPZBndCond *matLambdaBC_left_2 = mat2->CreateBC(mat2, fmatLambdaBC_left_2, fdirichlet_sigma, val1, val2);
            matLambdaBC_left_2->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_left_2);

            TPZBndCond *matLambdaBC_right = mat1->CreateBC(mat1, fmatLambdaBC_right, fdirichlet_sigma, val1, val2);
            matLambdaBC_right->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_right);

            TPZBndCond *matLambdaBC_right_2 = mat2->CreateBC(mat2, fmatLambdaBC_right_2, fdirichlet_sigma, val1, val2);
            matLambdaBC_right_2->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_right_2);


        }
            break;

        default: {

            TPZBndCond *matLambdaBC_bott = mat1->CreateBC(mat1, fmatLambdaBC_bott, fdirichlet_sigma, val1, val2);
            matLambdaBC_bott->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_bott);

            TPZBndCond *matLambdaBC_top = mat1->CreateBC(mat1, fmatLambdaBC_top, fdirichlet_sigma, val1, val2);
            matLambdaBC_top->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_top);

            TPZBndCond *matLambdaBC_left = mat1->CreateBC(mat1, fmatLambdaBC_left, fdirichlet_sigma, val1, val2);
            matLambdaBC_left->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_left);

            TPZBndCond *matLambdaBC_right = mat1->CreateBC(mat1, fmatLambdaBC_right, fdirichlet_sigma, val1, val2);
            matLambdaBC_right->SetBCForcingFunction(0, solp);
            cmesh.InsertMaterialObject(matLambdaBC_right);
//qwqwqwqw
            if (f_3Dmesh) {
                TPZBndCond *matLambdaBC_top_z = mat1->CreateBC(mat1, fmatLambdaBC_top_z, fneumann_sigma, val1, val2);
                matLambdaBC_top_z->SetBCForcingFunction(0, solp);
                cmesh.InsertMaterialObject(matLambdaBC_top_z);

                TPZBndCond *matLambdaBC_bott_z = mat1->CreateBC(mat1, fmatLambdaBC_bott_z, fneumann_sigma, val1, val2);
                matLambdaBC_bott_z->SetBCForcingFunction(0, solp);
                cmesh.InsertMaterialObject(matLambdaBC_bott_z);
            }

        }
    }


}


void MHMNavierStokesTest::InsertInterfaces(TPZMultiphysicsCompMesh *cmesh_m){
    
    std::set<int> boundaries_ids;
    boundaries_ids.insert(fmatBCbott);
    boundaries_ids.insert(fmatBCleft);
    boundaries_ids.insert(fmatBCtop);
    boundaries_ids.insert(fmatBCright);
    if (f_3Dmesh) {
        boundaries_ids.insert(fmatBCbott_z);
        boundaries_ids.insert(fmatBCtop_z);
    }
    
    
    TPZInterfaceInsertion InterfaceInsertion(cmesh_m, fmatLambda, boundaries_ids, feltype);
    TPZManVector<int64_t,3> Interfaces(2,0);
    Interfaces[0] = fmatInterfaceLeft;
    Interfaces[1] = fmatInterfaceRight;
    InterfaceInsertion.SetInterfaceVectorId(Interfaces);
    
    if (f_allrefine) {
        InterfaceInsertion.AddMultiphysicsInterfacesLeftNRight(fmatLambda);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_bott,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_top,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_left,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface(fmatLambdaBC_right,fmatInterfaceLeft);
    }else{
        InterfaceInsertion.AddMultiphysicsInterfacesLeftNRight2(fmatLambda);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_bott,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_top,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_left,fmatInterfaceLeft);
        InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_right,fmatInterfaceLeft);
        if (f_3Dmesh) {
            InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_bott_z,fmatInterfaceLeft);
            InterfaceInsertion.AddMultiphysicsBCInterface2(fmatLambdaBC_top_z,fmatInterfaceLeft);
        }
        
    }
    
    
    
}

void MHMNavierStokesTest::ComputeSkelNeighbours(){
    
    if (!f_mesh0) {
        DebugStop();
    }
    
    int64_t nel = f_mesh0->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = f_mesh0->Element(el);
//        if(gel->HasSubElement()&&f_allrefine)
//        {
//            continue;
//        }
        if (gel->MaterialId() != fmatLambda) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();

        if (neighbour == gelside) {
            continue;
        }
        
        while (neighbour != gelside) {
            int neigh_matID = neighbour.Element()->MaterialId();
            if (neighbour.Element()->Dimension() == f_mesh0->Dimension() && neigh_matID==fmatID) {
                f_skellNeighs.Push(neighbour);
            }
            if(neighbour.Element()->HasSubElement()){
                break;
            }
            neighbour = neighbour.Neighbour();
        }
    }
    
}


void MHMNavierStokesTest::GroupAndCondense(TPZMultiphysicsCompMesh *cmesh_m){
   
    
    //Criando apropamento de elementos
    
    int64_t ncompel = cmesh_m->ElementVec().NElements();
    int dim = cmesh_m->Reference()->Dimension();

    std::vector<int64_t> GroupIndex;
    TPZStack<TPZElementGroup *> elgroups;
    int count = 0;
    int64_t index =0;
    
    for(int64_t el = 0; el < ncompel; el++){
        
        TPZCompEl *cel = cmesh_m->Element(el);
        if (cel->Dimension()!=dim) {
            continue;
        }
        //GroupIndex[el] = cel->Index();
        count++;
        GroupIndex.resize(count);
        GroupIndex[count-1]=cel->Index();
        TPZElementGroup *GroupEl = new TPZElementGroup(*cmesh_m,index);
        elgroups.Push(GroupEl);
        elgroups[count-1]->AddElement(cel);
    }
    

    //Inserindo as respectivas interfaces e condições de contorno
    
    for(int64_t el = 0; el < ncompel; el++){
        TPZCompEl *cel = cmesh_m->Element(el);
        
        TPZMultiphysicsInterfaceElement *interel = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (interel) {
            TPZCompEl *Leftel = interel->LeftElement();
            
            if (Leftel->Dimension()!=dim) {
                continue;
            }
            int leftindex = Leftel->Index();
            
            for(int64_t iel = 0; iel < GroupIndex.size(); iel++){
                if (leftindex==GroupIndex[iel]) {
                    elgroups[iel]->AddElement(cel);
                }
            }
        }
        
        if (!cel) {
            continue;
        }
        
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension()==dim-1) {
            
            TPZBndCond *elBC = dynamic_cast<TPZBndCond *>(cel->Material());
            if (!elBC) {
                continue;
            }
            
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel, gel->NSides() - 1);
            
            gelside.EqualLevelCompElementList(celstack, 0, 0);
            
            for (auto &celstackindex : celstack) {
                if (celstackindex.Reference().Element()->Dimension() == dim) {
                    int bcindex = celstackindex.Element()->Index();
                    
                    for(int64_t iel = 0; iel < GroupIndex.size(); iel++){
                        if (bcindex==GroupIndex[iel]) {
                            elgroups[iel]->AddElement(cel);
                        }
                    }
                }
            }
        }
    }
    
    cmesh_m->ComputeNodElCon();
    // create condensed elements
    // increase the NumElConnected of one pressure connects in order to prevent condensation
//    for (int64_t ienv=0; ienv<nenvel; ienv++) {
//        TPZElementGroup *elgr = elgroups[ienv];
    int nenvel = elgroups.NElements();
    for (int64_t ienv=0; ienv<nenvel; ienv++) {
        TPZElementGroup *elgr = elgroups[ienv];
        
        int nc = elgroups[ienv]->GetElGroup()[0]->NConnects();
        elgroups[ienv]->GetElGroup()[0]->Connect(nc-1).IncrementElConnected();
        
        
//        for (int ic=0; ic<nc; ic++) {
//            TPZConnect &c = elgr->Connect(ic);
//            int connectpM = elgroups[ienv]->GetElGroup()[0]->NConnects();
//            int nc = elgr->NConnects();
//            TPZConnect &c = elgr->Connect(nc-1);
//            if (c.LagrangeMultiplier() > 0) {
//                c.IncrementElConnected();
//                break;
//            }
//        }
        new TPZCondensedCompEl(elgr);
    }

    
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
    
}

void MHMNavierStokesTest::InsertArcInterface(TPZAutoPointer<TPZGeoMesh> gmesh) {

    //Build Arcs in BC holes

    int nnodes = 0, matSkeleton = 4;
    int64_t nodeindex = 0, elementid = 0;;
    TPZManVector<REAL,3> coord0(3,0.),coord1(3,0.),coordHole(3,0),coordNode(3,0);
    int nel = gmesh->NElements();
    for(int iel = 0; iel < nel; iel++) {
        TPZGeoEl *gel = gmesh->ElementVec()[iel];
        if(!gel){
            continue;
        }
        if(gel->MaterialId()!=matSkeleton){
            continue;
        }

        TPZStack<TPZCompElSide> neigh;
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gelside.Neighbour();

        int matIdNeighs = 0;
        while (neighbour != gelside) {
            if(neighbour.Element()->Dimension()==2){
                matIdNeighs+= neighbour.Element()->MaterialId();
            }
            neighbour = neighbour.Neighbour();
        }

        if(matIdNeighs!=3){ // Ids 1 + 2 (Stokes - Darcy)
            continue;
        }

        gel->ResetReference();

        TPZManVector<int64_t> nodeindices;
        TPZManVector<int64_t,3> TopolArc(3);
        gel->GetNodeIndices(nodeindices);
        TopolArc[0]=nodeindices[0];
        TopolArc[1]=nodeindices[1];


        //Add new node:
        nnodes = gmesh->NNodes();
        gmesh->NodeVec().Resize(nnodes+1);
        nodeindex = nnodes;
        gmesh->SetNodeIdUsed(nodeindex);
        //REAL radius = h_s[0]/(4.*n_div[0]);
        REAL radius = 0.5;

        while (neighbour != gelside) {
            if(neighbour.Element()->Dimension()==2){
                int index_neigh = neighbour.Element()->Index();
                coordHole = f_HoleCoord[index_neigh];
                f_ArcCentralNode[index_neigh] = nodeindex;
                break;
            }
            neighbour = neighbour.Neighbour();
        }

        coordHole[0]=1.;
        coordHole[1]=1.;

        gel->Node(0).GetCoordinates(coord0);
        gel->Node(1).GetCoordinates(coord1);

        REAL rel = (coord0[1]-coordHole[1])/radius;
        REAL theta0H = asin(rel);
        rel = (coord1[1]-coordHole[1])/radius;
        REAL theta1H = asin(rel);
        REAL theta = theta0H - (theta0H-theta1H)/2;

        if (fabs(coord1[0]-coordHole[0])<1.e-6) {
            coord1[0] = coordHole[0];
        }
        if (fabs(coord0[1]-coordHole[1])<1.e-6) {
            coord0[1] = coordHole[1];
        }

        if(coord0[0]<coordHole[0]||coord1[0]<coordHole[0]){
            theta = Pi - theta;
        }

        coordNode = ParametricCircle(radius, theta);
        coordNode[0] +=coordHole[0];
        coordNode[1] +=coordHole[1];
        gmesh->NodeVec()[nodeindex].SetCoord(coordNode);
        gmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);

        TopolArc[2]=nodeindex;

        elementid = gel->Id();

        //delete gmesh->ElementVec()[iel];
        gmesh->DeleteElement(gel);
        gmesh->ElementVec()[iel] = new TPZGeoElRefPattern< pzgeom::TPZArc3D > (elementid,TopolArc, matSkeleton,*gmesh);
        // gmesh->ElementVec()[elementid]->IncrementNumInterfaces();

    }
    gmesh->BuildConnectivity();


}

void MHMNavierStokesTest::GetElIndexCoarseMesh(TPZGeoMesh *gmesh, TPZVec<int64_t> &coarseindex)
{
    int nel = gmesh->NElements();
    int iel;
    int hassubel=0;
    int dim = gmesh->Dimension();
    int eldim;
    int count =0.;
    
    for(iel = 0; iel<nel; iel++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[iel];
        if(!gel) continue;
        
        hassubel = gel->HasSubElement();
        eldim = gel->Dimension();
        if(!hassubel && eldim ==dim)
        {
            count++;
            coarseindex.resize(count);
            coarseindex[count-1] = gel->Index();
        }
    }
    
}

