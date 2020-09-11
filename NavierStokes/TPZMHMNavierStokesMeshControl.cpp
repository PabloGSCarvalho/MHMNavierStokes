//
//  TPZMHMixedMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMNavierStokesMeshControl.h"
#include "TPZMHMeshControl.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
#include "TPZMatLaplacian.h"
#include "TPZLagrangeMultiplier.h"

#include "pzsmanal.h"
#include "pzanalysis.h"
#include "TPZSpStructMatrix.h"
#include "pzstepsolver.h"

#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzcompelwithmem.h"
#include "TPZVTKGeoMesh.h"
#include "TPZNullMaterial.h"

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr logger(Logger::getLogger("pz.mhmNavierStokesmeshcontrol"));
#endif


TPZMHMNavierStokesMeshControl::TPZMHMNavierStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMixedMeshControl(gmesh,coarseindices), fBCTractionMatIds()
{
    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCoarseAveragePressMesh = new TPZCompMesh(fGMesh);
    fCoarseDistrFluxMesh = new TPZCompMesh(fGMesh);
    fsetCoarseAverageMultipliers = false;
    fsetStaticCondensation = true;
    
    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
    
}

TPZMHMNavierStokesMeshControl::TPZMHMNavierStokesMeshControl(int dimension) : TPZMHMixedMeshControl(dimension), fBCTractionMatIds(){

    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCoarseAveragePressMesh = new TPZCompMesh(fGMesh);
    fCoarseDistrFluxMesh = new TPZCompMesh(fGMesh);
    fsetCoarseAverageMultipliers = false;
    fsetStaticCondensation = true;
    
    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
}

TPZMHMNavierStokesMeshControl::TPZMHMNavierStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMixedMeshControl(gmesh), fBCTractionMatIds(){

    fAveragePressMesh = new TPZCompMesh(fGMesh);
    fDistrFluxMesh = new TPZCompMesh(fGMesh);
    fCoarseAveragePressMesh = new TPZCompMesh(fGMesh);
    fCoarseDistrFluxMesh = new TPZCompMesh(fGMesh);
    fsetCoarseAverageMultipliers = false;
    fsetStaticCondensation = true;

    fBCTractionMatIds.clear();
    for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
    {
        fBCTractionMatIds[*it]=*it-10;
    }
}


TPZMHMNavierStokesMeshControl::TPZMHMNavierStokesMeshControl(const TPZMHMNavierStokesMeshControl &copy) : TPZMHMixedMeshControl(copy){
    
    this->operator=(copy);
}

TPZMHMNavierStokesMeshControl &TPZMHMNavierStokesMeshControl::operator=(const TPZMHMNavierStokesMeshControl &cp){

    fDistrFluxMesh = cp.fDistrFluxMesh;
    fAveragePressMesh = cp.fAveragePressMesh;
    fCoarseAveragePressMesh = cp.fCoarseAveragePressMesh;
    fCoarseDistrFluxMesh = cp.fCoarseDistrFluxMesh;
    fBCTractionMatIds = cp.fBCTractionMatIds;
    fsetCoarseAverageMultipliers = cp.fsetCoarseAverageMultipliers;
    fsetStaticCondensation = cp.fsetStaticCondensation;
    return *this;
}


TPZMHMNavierStokesMeshControl::~TPZMHMNavierStokesMeshControl()
{

}

void TPZMHMNavierStokesMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    BuildWrapMesh(fGMesh->Dimension());
    fGeoToMHMDomain.Resize(fGMesh->NElements(), -1);

    //InsertPeriferalMaterialObjects();  //Skeleton com dirichlet v=0
    CreateHDivMHMMesh();
    
    InsertBCSkeleton();
    InsertInternalSkeleton();

#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    
    if(1){
        std::ofstream fileg1("MalhaGeo_test.txt"); //Impressão da malha geométrica (formato txt)
        std::ofstream filegvtk1("MalhaGeo_test.vtk"); //Impressão da malha geométrica (formato vtk)
        fGMesh->Print(fileg1);
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh, filegvtk1,true);
    }

    InsertPeriferalPressureMaterialObjects();
    CreatePressureAndTractionMHMMesh();
    
    InsertDistributedFluxMaterialObjects();
    CreateDistributedFluxMHMMesh();
   
    InsertPeriferalAveragePressMaterialObjects();
    CreateAveragePressMHMMesh();

    if (fsetCoarseAverageMultipliers) {
        CreateLagrangeMultiplierMesh();

        int ncL1 = fCMeshLagrange.operator->()->NConnects();
        for(int i_con1 =0; i_con1 < ncL1; i_con1++){
            fCMeshLagrange.operator->()->ConnectVec()[i_con1].SetLagrangeMultiplier(3);
        }
        int ncL2 = fCMeshConstantPressure.operator->()->NConnects();
        for(int i_con2 =0; i_con2 < ncL2; i_con2++) {
            fCMeshConstantPressure.operator->()->ConnectVec()[i_con2].SetLagrangeMultiplier(3);
        }
    }

  //  CreateCoarseAveragePressMHMMesh();
  //  CreateCoarseDistributedFluxMHMMesh();
    
    CreateMultiPhysicsMHMMesh();
//    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
//    std::cout << "Number of connects (obs) " << fGlobalSystemWithLocalCondensationSize << std::endl;
    fGlobalSystemSize = fCMesh->Solution().Rows();
    
    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG
    if(0){
        std::ofstream out("Friendly.txt");
        PrintFriendly(out);
    }
    CheckMeshConsistency();
#endif

    //if (usersubstructure) {
    //HideTheElements();
    //GroupAndCondense(fCMesh.operator->());
    //}
    fNumeq = fCMesh->NEquations();
    
   
    
#ifdef PZDEBUG
    if(0){
        int64_t nel = fCMesh->NElements();
        for(int64_t el = 0; el<nel; el++)
        {
            TPZCompEl *cel = fCMesh->Element(el);
            TPZSubCompMesh *sub = dynamic_cast<TPZSubCompMesh *>(cel);
            if(sub)
            {
                std::stringstream sout;
                sout << "submesh_" << el << ".vtk";
                std::ofstream file(sout.str());
                TPZVTKGeoMesh::PrintCMeshVTK(sub, file,true);
            }
        }
        
    }
#endif
    
}

void TPZMHMNavierStokesMeshControl::CreatePressureAndTractionMHMMesh(){
    
    CreatePressureMHMMesh();
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    
    int64_t nskeletonconnects = fPressureFineMesh->NConnects();
    int porder = fpOrderInternal;
    TPZCompMesh * cmeshTraction = fPressureFineMesh.operator->();
    gmesh->ResetReference();
    cmeshTraction->SetName("PressureAndTractionMesh");
    cmeshTraction->SetDimModel(gmesh->Dimension()-1);
    cmeshTraction->SetAllCreateFunctionsDiscontinuous();
    cmeshTraction->SetDefaultOrder(porder-1);
    int meshdim = cmeshTraction->Dimension();
    
    std::set<int> matids;
    TPZMaterial *mat = cmeshTraction->FindMaterial(fTractionMatId);
    if (mat && mat->Dimension() == meshdim) {
        matids.insert(fTractionMatId);
    }
    
    for (auto it:fMaterialBCIds) {
        int dsmatid = fBCTractionMatIds[it];
        TPZMaterial *mat = cmeshTraction->FindMaterial(fBCTractionMatIds[it]);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(fBCTractionMatIds[it]);
        }
    }
    
    cmeshTraction->AutoBuild(matids);
    fPressureFineMesh->ExpandSolution();

    

#ifdef PZDEBUG
    {
        int64_t nel = fGMesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (gel && gel->MaterialId() == 1) {
                if (gel->Dimension() != fGMesh->Dimension()) {
                    DebugStop();
                }
            }
        }
    }
#endif
    
    int64_t nc = cmeshTraction->NConnects();
    if(nskeletonconnects != 0){
  //      DebugStop();
    }
    for (int64_t ic=nskeletonconnects; ic<nc; ic++) {
        cmeshTraction->ConnectVec()[ic].SetLagrangeMultiplier(2);
    }
    gmesh->ResetReference();
    
    
    // associate the connects with the proper subdomain aoaoaoaoa
    int64_t nel = cmeshTraction->NElements();

    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshTraction->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();

        int dim = gmesh->Dimension();

        if (gel->Dimension()==dim) {
            continue;
        }
        int gelindex = gel->Index();

        int nsides = gel->NSides();

        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gel->Neighbour(nsides-1);

        int notintern = 0;

        while(neighbour != gelside)
        {
            if(neighbour.Element()->MaterialId()==fInternalWrapMatId) {
                notintern = 1;
                break;
            }
            neighbour = neighbour.Neighbour();
        }

        if (notintern==1) {

            neighbour = gel->Neighbour(nsides-1);

            while(neighbour != gelside)
            {
                if (neighbour.Element()->Dimension() == dim) {
                    break;
                }
                neighbour = neighbour.Neighbour();
            }

            int domain = fGeoToMHMDomain[neighbour.Element()->Index()];
#ifdef PZDEBUG
            if (domain == -1) {
                DebugStop();
            }
#endif

            SetSubdomain(cel, domain);
        }else{
            SetSubdomain(cel, -1);
        }

    }

    if(0){
        std::ofstream out("PressureAndTractionFineMesh.txt");
        fPressureFineMesh->Print(out);
    }

}

void TPZMHMNavierStokesMeshControl::InsertInternalSkeleton(){
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (gel->Dimension()!=fGMesh->Dimension()-1) {
            continue;
        }
        if (gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        
        if (gel->MaterialId()==fSkeletonMatId||gel->MaterialId()==fInternalWrapMatId) {
            TPZGeoElBC(gelside, fTractionMatId);
        }
        
    }

}

void TPZMHMNavierStokesMeshControl::InsertBCSkeleton(){
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        
        if(gel->HasSubElement()){
            continue;
        }
        
        for(auto it = fMaterialBCIds.begin(); it != fMaterialBCIds.end(); it++)
        {
            if (gel->MaterialId()==*it) {
                TPZGeoElBC(gelside, fBCTractionMatIds[*it]);
            }
        }
    }
    
}


void TPZMHMNavierStokesMeshControl::InsertPeriferalPressureMaterialObjects()
{
    TPZCompMesh * cmeshPressure = fPressureFineMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshPressure->MaterialVec().find(matid) == cmeshPressure->MaterialVec().end())
        {
            TPZNullMaterial *matl2 = new TPZNullMaterial((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshPressure->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
    // Material for interior traction:
    
    TPZVecL2 *matTraction = new TPZVecL2(fTractionMatId);
    matTraction->SetDimension(fGMesh->Dimension()-1);
    cmeshPressure->InsertMaterialObject(matTraction);

    for (auto it:fMaterialBCIds)
    {
        if (fBCTractionMatIds.size()!=fMaterialBCIds.size()) {
            DebugStop();
        }
        
        int matid= fBCTractionMatIds[it];
        if (cmeshPressure->MaterialVec().find(matid) == cmeshPressure->MaterialVec().end())
        {
            TPZVecL2 *matBCTraction = new TPZVecL2(matid);
            matBCTraction->SetDimension(fGMesh->Dimension()-1);
            cmeshPressure->InsertMaterialObject(matBCTraction);
        }
        
    }
    
}

void TPZMHMNavierStokesMeshControl::InsertPeriferalAveragePressMaterialObjects(){
    
    TPZCompMesh *cmeshAverPressure = fAveragePressMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshAverPressure->MaterialVec().find(matid) == cmeshAverPressure->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshAverPressure->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
    
}

void TPZMHMNavierStokesMeshControl::CreateAveragePressMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshAverPressute = fAveragePressMesh.operator->();
    gmesh->ResetReference();
    cmeshAverPressute->SetName("AveragePressureMesh");
    cmeshAverPressute->SetDimModel(gmesh->Dimension());
    cmeshAverPressute->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshAverPressute->SetDefaultOrder(porder);
    
    int meshdim = cmeshAverPressute->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshAverPressute->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshAverPressute->AutoBuild(matids);
    fAveragePressMesh->ExpandSolution();

    
    int64_t nel = fAveragePressMesh->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmeshAverPressute->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) DebugStop();
        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
        {
            DebugStop();
        }
        celdisc->SetTotalOrderShape();
        celdisc->SetFalseUseQsiEta();
    }
    
    int64_t nc = cmeshAverPressute->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshAverPressute->ConnectVec()[ic].SetLagrangeMultiplier(3);
    }
    gmesh->ResetReference();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshAverPressute->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            continue;
        }
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
    }
    
    if(0)
    {
        std::ofstream out("AveragePressureMesh.txt");
        fAveragePressMesh->Print(out);
    }
    
    
}

void TPZMHMNavierStokesMeshControl::InsertPeriferalCoarseAveragePressMaterialObjects(){
    
    TPZCompMesh *cmeshCoarseAverPressure = fCoarseAveragePressMesh.operator->();
    
    TPZVecL2 *matl2 = new TPZVecL2(fSkeletonMatId);
    matl2->SetNStateVariables(fNState);
    matl2->SetDimension(fGMesh->Dimension());
    cmeshCoarseAverPressure->InsertMaterialObject(matl2);
    
}

void TPZMHMNavierStokesMeshControl::CreateCoarseAveragePressMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshAverPressute = fCMeshConstantPressure.operator->();
    gmesh->ResetReference();
    cmeshAverPressute->SetName("CoarseAveragePressureMesh");
    cmeshAverPressute->SetDimModel(gmesh->Dimension());
    cmeshAverPressute->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshAverPressute->SetDefaultOrder(porder);
    
    int meshdim = cmeshAverPressute->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshAverPressute->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshAverPressute->AutoBuild(matids);
    fAveragePressMesh->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("CoarseAveragePressureMesh.txt");
        fAveragePressMesh->Print(out);
    }
    
    
//    int64_t nel = fAveragePressMesh->NElements();
//    for(int64_t i=0; i<nel; i++){
//        TPZCompEl *cel = cmeshAverPressute->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(!celdisc) DebugStop();
//        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
//        {
//            DebugStop();
//        }
//        celdisc->SetTotalOrderShape();
//        celdisc->SetFalseUseQsiEta();
//    }
    
    int64_t nc = cmeshAverPressute->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshAverPressute->ConnectVec()[ic].SetLagrangeMultiplier(2);
    }
    gmesh->ResetReference();
    
//    for (int64_t el=0; el<nel; el++)
//    {
//        TPZCompEl *cel = cmeshAverPressute->Element(el);
//#ifdef PZDEBUG
//        if (! cel) {
//            DebugStop();
//        }
//#endif
//        TPZGeoEl *gel = cel->Reference();
//        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
//        {
//            continue;
//        }
//#ifdef PZDEBUG
//        if (fGeoToMHMDomain[gel->Index()] == -1) {
//            DebugStop();
//        }
//#endif
//
//        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
//    }
    
    
    return;
    
}



void TPZMHMNavierStokesMeshControl::InsertDistributedFluxMaterialObjects(){
    
    TPZCompMesh *cmeshDistributedFlux = fDistrFluxMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshDistributedFlux->MaterialVec().find(matid) == cmeshDistributedFlux->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshDistributedFlux->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
}


void TPZMHMNavierStokesMeshControl::CreateDistributedFluxMHMMesh(){

    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshDistributedFlux = fDistrFluxMesh.operator->();
    gmesh->ResetReference();
    cmeshDistributedFlux->SetName("DistributedFluxMesh");
    cmeshDistributedFlux->SetDimModel(gmesh->Dimension());
    cmeshDistributedFlux->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshDistributedFlux->SetDefaultOrder(porder);
    
    int meshdim = cmeshDistributedFlux->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshDistributedFlux->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshDistributedFlux->AutoBuild(matids);
    fDistrFluxMesh->ExpandSolution();
    
    
    int64_t nel = fDistrFluxMesh->NElements();
    for(int64_t i=0; i<nel; i++){
        TPZCompEl *cel = cmeshDistributedFlux->ElementVec()[i];
        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
        if(!celdisc) DebugStop();
        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
        {
            DebugStop();
        }
        celdisc->SetTotalOrderShape();
        celdisc->SetFalseUseQsiEta();
    }
    
    int64_t nc = cmeshDistributedFlux->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshDistributedFlux->ConnectVec()[ic].SetLagrangeMultiplier(2);
    }
    gmesh->ResetReference();
    
    for (int64_t el=0; el<nel; el++)
    {
        TPZCompEl *cel = cmeshDistributedFlux->Element(el);
#ifdef PZDEBUG
        if (! cel) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
        {
            continue;
        }
#ifdef PZDEBUG
        if (fGeoToMHMDomain[gel->Index()] == -1) {
            DebugStop();
        }
#endif
        
        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
    }
    
    
    if(0)
    {
        std::ofstream out("DistributedFluxMesh.txt");
        fDistrFluxMesh->Print(out);
    }
    
}

void TPZMHMNavierStokesMeshControl::InsertCoarseDistributedFluxMaterialObjects(){
    
    TPZCompMesh *cmeshDistributedFlux = fDistrFluxMesh.operator->();
    
    for (auto it = fMaterialIds.begin(); it != fMaterialIds.end(); it++)
    {
        int matid = *it;
        if (cmeshDistributedFlux->MaterialVec().find(matid) == cmeshDistributedFlux->MaterialVec().end())
        {
            TPZVecL2 *matl2 = new TPZVecL2((matid));
            matl2->SetNStateVariables(fNState);
            matl2->SetDimension(fGMesh->Dimension());
            cmeshDistributedFlux->InsertMaterialObject(matl2);
        }
        else
        {
            DebugStop();
        }
    }
    
}

void TPZMHMNavierStokesMeshControl::CreateCoarseDistributedFluxMHMMesh(){
    
    TPZGeoMesh * gmesh = fGMesh.operator->();
    gmesh->ResetReference();
    int porder = 0; //constante
    TPZCompMesh * cmeshDistributedFlux = fCMeshLagrange.operator->();
    gmesh->ResetReference();
    cmeshDistributedFlux->SetName("CoarseDistributedFluxMesh");
    cmeshDistributedFlux->SetDimModel(gmesh->Dimension());
    cmeshDistributedFlux->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    cmeshDistributedFlux->SetDefaultOrder(porder);
    
    int meshdim = cmeshDistributedFlux->Dimension();
    std::set<int> matids;
    for (auto it:fMaterialIds) {
        TPZMaterial *mat = cmeshDistributedFlux->FindMaterial(it);
        if (mat && mat->Dimension() == meshdim) {
            matids.insert(it);
        }
    }
    cmeshDistributedFlux->AutoBuild(matids);
    fDistrFluxMesh->ExpandSolution();
    
    if(0)
    {
        std::ofstream out("CoarseDistributedFluxMesh.txt");
        fDistrFluxMesh->Print(out);
    }
    
    
//    int64_t nel = fDistrFluxMesh->NElements();
//    for(int64_t i=0; i<nel; i++){
//        TPZCompEl *cel = cmeshDistributedFlux->ElementVec()[i];
//        TPZCompElDisc *celdisc = dynamic_cast<TPZCompElDisc *>(cel);
//        if(!celdisc) DebugStop();
//        if(celdisc && celdisc->Reference()->Dimension() != meshdim)
//        {
//            DebugStop();
//        }
//        celdisc->SetTotalOrderShape();
//        celdisc->SetFalseUseQsiEta();
//    }
    
    int64_t nc = cmeshDistributedFlux->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
        cmeshDistributedFlux->ConnectVec()[ic].SetLagrangeMultiplier(2);
    }
    gmesh->ResetReference();
    
//    for (int64_t el=0; el<nel; el++)
//    {
//        TPZCompEl *cel = cmeshDistributedFlux->Element(el);
//#ifdef PZDEBUG
//        if (! cel) {
//            DebugStop();
//        }
//#endif
//        TPZGeoEl *gel = cel->Reference();
//        if(fMaterialIds.find (gel->MaterialId()) == fMaterialIds.end())
//        {
//            continue;
//        }
//#ifdef PZDEBUG
//        if (fGeoToMHMDomain[gel->Index()] == -1) {
//            DebugStop();
//        }
//#endif
//
//        SetSubdomain(cel, fGeoToMHMDomain[gel->Index()]);
//    }
    
    return;
    
}


void TPZMHMNavierStokesMeshControl::CreateMultiPhysicsMHMMesh()
{
    TPZManVector<TPZCompMesh *,6 > cmeshes(4);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    cmeshes[2] = fDistrFluxMesh.operator->();
    cmeshes[3] = fAveragePressMesh.operator->();
    if (fsetCoarseAverageMultipliers) {
        cmeshes.Resize(6);
        cmeshes[4] = fCMeshLagrange.operator->();
        cmeshes[5] = fCMeshConstantPressure.operator->();
    }
    
    TPZGeoMesh *gmesh = cmeshes[0]->Reference();
    if(!gmesh)
    {
        std::cout<< "Geometric mesh doesn't exist" << std::endl;
        DebugStop();
    }
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    
    // Malha computacional
    TPZCompMesh * MixedFluxPressureCmesh = fCMesh.operator->();
    MixedFluxPressureCmesh->SetDimModel(dim);
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElem();
    MixedFluxPressureCmesh->SetAllCreateFunctionsMultiphysicElemWithMem();
    //gSinglePointMemory = true;
    
    BuildMultiPhysicsMesh();

    {
        int dim = MixedFluxPressureCmesh->Dimension();
        int64_t nelem = MixedFluxPressureCmesh->NElements();
        for(int64_t el=0; el<nelem; el++)
        {
            TPZCompEl *cel = MixedFluxPressureCmesh->Element(el);
            if(!cel) DebugStop();
            TPZGeoEl *gel = cel->Reference();
            if(!gel) DebugStop();
//            if(gel->Dimension() == dim)
//            {
                cel->PrepareIntPtIndices();
//            }
        }
    }

    TPZManVector<TPZCompMesh * ,6> meshvector;
    
    if(0)
    {
        std::ofstream out2("gmesh.txt");
        gmesh->Print(out2);
        std::ofstream out3("HDivMesh.txt");
        fFluxMesh->Print(out3);
        std::ofstream out4("PressureMesh.txt");
        fPressureFineMesh->Print(out4);
    }
    
    meshvector = cmeshes;
    
    JoinSubdomains(meshvector, MixedFluxPressureCmesh);
    
    // Transferindo para a multifisica
    //TPZBuildMultiphysicsMesh::AddElements(meshvector, MixedFluxPressureCmesh);
    //TPZBuildMultiphysicsMesh::AddConnects(meshvector, MixedFluxPressureCmesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvector, MixedFluxPressureCmesh);
    
    MixedFluxPressureCmesh->LoadReferences();
    
    std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
    //CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-1);
    //CreateMultiPhysicsInterfaceElements(fGMesh->Dimension()-2);
    
    CreateMultiPhysicsInterfaceElements();
    CreateMultiPhysicsBCInterfaceElements();

//    std::map<int64_t,int64_t> submeshindices;
//    TPZCompMeshTools::PutinSubmeshes(mixed_cmesh, ElementGroups, submeshindices, KeepOneLagrangian);

#ifdef PZDEBUG
    if(0) {
        std::ofstream file("multiphysics.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(MixedFluxPressureCmesh, file, true);
        std::ofstream out("multiphysics_without_condensation.txt");
        MixedFluxPressureCmesh->Print(out);
    }
#endif

    if(0) {
        std::ofstream out2("gmesh_elementos.txt");
        gmesh->Print(out2);
    }


    std::cout << "Number of DoFs (without condesation) = " << fCMesh->NEquations() << std::endl;

    if(fsetStaticCondensation){
        BuildSubMeshes(); //aoaoaoa
        std::cout << "Number of DoFs (after condesation) = " << fCMesh->NEquations() << std::endl;
        fGlobalSystemSize = fCMesh->Solution().Rows();

#ifdef PZDEBUG
        if(0)
        {
            std::ofstream out("multiphysics_after_condensation.txt");
            MixedFluxPressureCmesh->Print(out);
        }
#endif
    }

    //GroupAndCondense(MixedFluxPressureCmesh);
    MixedFluxPressureCmesh->CleanUpUnconnectedNodes();
    //MixedFluxPressureCmesh->ExpandSolution();
    //GroupAndCondense(MixedFluxPressureCmesh);

    
    return;
    
}

void TPZMHMNavierStokesMeshControl::BuildSubMeshes(){
    
    bool KeepOneLagrangian = true;
    if (fHybridize) {
        KeepOneLagrangian = false;
    }
    typedef std::set<int64_t> TCompIndexes;
    std::map<int64_t, TCompIndexes> ElementGroups;
    TPZGeoMesh *gmesh = fCMesh->Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    fCMesh->LoadReferences();
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int el_index = gel->Index();

        if (cel->Reference()->HasSubElement()) {
            continue;
        }
        TPZBndCond *elBC = dynamic_cast<TPZBndCond *>(cel->Material());
        if (elBC) {
            continue;
        }

        int nsides=cel->Reference()->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZGeoElSide neighbour = gel->Neighbour(nsides-1);
        
        int notintern = 0;
        
//        while(neighbour != gelside)
//        {
//            if(neighbour.Element()->MaterialId()!=fInternalWrapMatId) {
//                notintern = 1;
//                break;
//            }
//            neighbour = neighbour.Neighbour();
//        }
//        
//        if (notintern == 1) {
//        //    continue;
//        }

        
        int64_t domain = WhichSubdomain(cel);
        if (domain == -1) {
            continue;
        }
        ElementGroups[domain].insert(el);
    }
    
    
#ifdef PZDEBUG
    if (ElementGroups.size() <= 5) {
        std::map<int64_t, TCompIndexes>::iterator it;
        for (it = ElementGroups.begin(); it != ElementGroups.end(); it++) {
            std::cout << "Group " << it->first << " with size " << it->second.size() << std::endl;
            std::cout << " elements ";
            std::set<int64_t>::iterator its;
            for (its = it->second.begin(); its != it->second.end(); its++) {
                std::cout << *its << " ";
            }
            std::cout << std::endl;
        }
    }
#endif
    KeepOneLagrangian = 3;
    std::map<int64_t,int64_t> submeshindices;


//    int64_t nc = fCMesh.operator->()->NConnects();
//    fCMesh.operator->()->ConnectVec()[nc-2].SetLagrangeMultiplier(3);
//    fCMesh.operator->()->ConnectVec()[nc-1].SetLagrangeMultiplier(3);

    PutinSubmeshes(fCMesh.operator->(), ElementGroups, submeshindices, 3);

    //std::cout << "Number of connects (1) " << fCMesh->NEquations() << std::endl;

    std::cout << "After putting in substructures\n";
    fMHMtoSubCMesh = submeshindices;

    fCMesh->ComputeNodElCon();
    fCMesh->CleanUpUnconnectedNodes();

    GroupandCondenseSubMeshes();
//    GroupandCondenseElements();

    std::cout << "Finished substructuring\n";
    
    
}

/// Put the element set into a subcompmesh and make the connects internal
void TPZMHMNavierStokesMeshControl::PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, int KeepOneLagrangian)
{
    for (std::map<int64_t,std::set<int64_t> >::iterator it = elindices.begin(); it != elindices.end(); it++) {
        int64_t index;
        TPZSubCompMesh *subcmesh = new TPZSubCompMesh(*cmesh,index);
        indices[it->first] = index;
        for (std::set<int64_t>::iterator itloc = it->second.begin(); itloc != it->second.end(); itloc++) {
            subcmesh->TransferElement(cmesh, *itloc);
        }
    }
    cmesh->ComputeNodElCon();
    for (std::map<int64_t,int64_t>::iterator it = indices.begin(); it != indices.end(); it++) {
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cmesh->Element(it->second));
        if (!subcmesh) {
            DebugStop();
        }
        int count = 0;
        if (KeepOneLagrangian)
        {
            int64_t nconnects = subcmesh->NConnects();
            for (int64_t ic=0; ic<nconnects; ic++) {
                TPZConnect &c = subcmesh->Connect(ic);
                if (c.LagrangeMultiplier() == KeepOneLagrangian) {
                    c.IncrementElConnected();
//                    count++;
//                    if(count == 6 && c.NState() == 1)
//                    {
//                        break;
//                    }
//                    else if(count == 2 && c.NState() == 2)
//                    {
//                        break;
//                    }
//                    else if(count == 3 && c.NState() == 3)
//                    {
//                        break;
//                    }
                }

                //    else if(c.LagrangeMultiplier() == 1){
//                    if(count >= 2 && c.NState() == 1)
//                    {
//                        continue;
//                    }
//                    c.IncrementElConnected();
//                    count++;
//                }

            }
        }


        subcmesh->MakeAllInternal();
    }



}

void TPZMHMNavierStokesMeshControl::CreateMultiPhysicsInterfaceElements(){
    
    TPZCompMesh *cmesh = fCMesh.operator->();
    TPZVec<int> m_interfaceVector_ids(2,0);
    m_interfaceVector_ids[0] = fLagrangeMatIdLeft;
    m_interfaceVector_ids[1] = fLagrangeMatIdRight;
    
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        int meshdim = fGMesh->Dimension();
        int matid = gel->MaterialId();
        
        if (matid != fTractionMatId) {
            continue;
        }
        
        if (gel->HasSubElement() == 1) {
            continue;
        }
        
        int nsides = gel->NSides();
        TPZGeoElSide gelside(gel,nsides-1);
        TPZCompElSide celside = gelside.Reference();
        
        if (!celside ) {
            DebugStop();
        }
        
        TPZStack<TPZGeoElSide> neighbourset;
       // gelside.AllNeighbours(neighbourset);
        
        //Find elements with the same mesh dimension :
        TPZGeoElSide neigh = gelside.Neighbour();
        while(neigh != gelside)
        {
            if (neigh.Element()->Dimension()!=fGMesh->Dimension()) {
                neigh = neigh.Neighbour();
                continue;
            }
            
            neighbourset.Push(neigh);
            neigh = neigh.Neighbour();
        }
        
        //gelside.LowerLevelCompElementList2(1);
        
        int nneighs = neighbourset.size();
        if (nneighs!=2) {
            DebugStop();
        }
        
        TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
        LeftElIndices[0]=0;
        RightElIndices[0]=1;
        
        for(int stack_i=0; stack_i <nneighs; stack_i++){
            TPZGeoElSide neigh = neighbourset[stack_i];
            
            if (neigh.Element()->Dimension()!=meshdim) {
                continue;
            }
            
            TPZCompElSide celneigh = neigh.Reference();

            
            if(neigh.Element()->HasSubElement()) {

                TPZStack<TPZGeoElSide> subel;
                neigh.YoungestChildren(subel);
                
                for (int i_sub =0; i_sub<subel.size(); i_sub++) {
                    
                    TPZCompElSide cel_sub_neigh = subel[i_sub].Reference();
                    
                    TPZGeoElBC gbc_sub(subel[i_sub],m_interfaceVector_ids[stack_i]);
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc_sub.CreatedElement(),index,cel_sub_neigh,celside);
                    elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                    
#ifdef PZDEBUG
                    if(0){
                        std::cout << "****Created an interface element between volumetric element " << subel[i_sub].Element()->Index() <<
                                  " side " << subel[i_sub].Side() <<
                                  " and Skeleton element " << gelside.Element()->Index() << std::endl;
                    }
#endif

                }
                
                
            }else{
                
                TPZGeoElBC gbc(gelside,m_interfaceVector_ids[stack_i]);
                int64_t index;
                
                TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
                elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);

#ifdef PZDEBUG
                if(0) {
                    std::cout << "Created an interface element between volumetric element " << neigh.Element()->Index()
                              <<
                              " side " << neigh.Side() <<
                              " and interior 1D element " << gelside.Element()->Index() << std::endl;
                }
#endif
                
            }
            
        }
        
    }
    
    
}

void TPZMHMNavierStokesMeshControl::CreateMultiPhysicsBCInterfaceElements(){
    
    int matBCinterface = fLagrangeMatIdLeft;
    for (auto it : fMaterialBCIds) {
        
        int matfrom = fBCTractionMatIds[it];
        TPZCompMesh *cmesh = fCMesh.operator->();
        
        int64_t nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            int meshdim = fGMesh->Dimension();
            int matid = gel->MaterialId();
            
            if (matid != matfrom) {
                continue;
            }
            
            int nsides = gel->NSides();
            TPZGeoElSide gelside(gel,nsides-1);
            TPZCompElSide celside = gelside.Reference();
            
            TPZStack<TPZGeoElSide> neighbourset;
            gelside.AllNeighbours(neighbourset);
            
            int nneighs = neighbourset.size();
//            if(nneighs!=2){
//                //    DebugStop();
//            }
            
            TPZManVector<int64_t,3> LeftElIndices(1,0.),RightElIndices(1,0.);
            LeftElIndices[0]=0;
            RightElIndices[0]=1;
            
            for(int stack_i=0; stack_i <nneighs; stack_i++){
                TPZGeoElSide neigh = neighbourset[stack_i];
                if (neigh.Element()->Dimension()!=meshdim) {
                    continue;
                }
                
                TPZCompElSide celneigh = neigh.Reference();
                if (!celside || !celneigh) {
                    //    DebugStop();
                }
                int64_t neigh_index = neigh.Element()->Index();
                if (neigh.Element()->Dimension()!=meshdim){
                    continue;
                }
                
                if (neigh.Element()->HasSubElement()) {
                    
                    TPZStack<TPZGeoElSide > subelside;
                    neigh.YoungestChildren(subelside);
                    
                    for (int i_sub = 0; i_sub<subelside.size(); i_sub++) {
                        
                        TPZCompElSide cel_sub_neigh = subelside[i_sub].Reference();
                        
                        TPZGeoElBC gbc_sub(subelside[i_sub],matBCinterface);
                        int64_t index;
                        
                        TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc_sub.CreatedElement(),index,cel_sub_neigh,celside);
                        elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);

#ifdef PZDEBUG
                        if(0){
                            std::cout << "****Created an BC interface element between volumetric element " << subelside[i_sub].Element()->Index() <<
                                      " side " << subelside[i_sub].Side() <<
                                      " and boundary 1D element " << gelside.Element()->Index() << std::endl;
                        }
#endif
                        
                    }
                    
                }else{
                    
                    TPZGeoElBC gbc(gelside,matBCinterface);
                    int64_t index;
                    
                    TPZMultiphysicsInterfaceElement *elem_inter = new TPZMultiphysicsInterfaceElement(*cmesh,gbc.CreatedElement(),index,celneigh,celside);
                    elem_inter->SetLeftRightElementIndices(LeftElIndices,RightElIndices);
                    
#ifdef PZDEBUG
                    if(0){
                        std::cout << "Created an BC interface element between volumetric element " << neigh.Element()->Index() <<
                                  " side " << neigh.Side() <<
                                  " and boundary 1D element " << gelside.Element()->Index() << std::endl;
                    }
#endif
                }
                
            }
        }
        
    }
    
    
}

void TPZMHMNavierStokesMeshControl::GroupandCondenseSubMeshes()
{
    for (std::map<int64_t,int64_t>::iterator it=fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
        TPZCompEl *cel = fCMesh->Element(it->second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        subcmesh->ComputeNodElCon();

//        TPZCompMeshTools::GroupElements(subcmesh);

//
//#ifdef LOG4CXX2
//        if(logger->isDebugEnabled())
//        {
//            std::stringstream sout;
//            subcmesh->Print(sout);
//            LOGPZ_DEBUG(logger, sout.str())
//        }
//#endif
//        // TODO: increment nelconnected of exterior connects

//        bool keeplagrange = true;
//        TPZCompMeshTools::CreatedCondensedElements(subcmesh, false, keeplagrange);

        if(fMHMtoSubCMesh.size()<5)
        {
            std::ofstream filesub("submesh.vtk");
            TPZVTKGeoMesh::PrintCMeshVTK(subcmesh, filesub,true);
            std::ofstream out("submesh.txt");
            subcmesh->Print(out);
        }

        // Increment nelconnected of exterior connects

        int nel = subcmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = subcmesh->Element(el);
            if (!cel) {
                continue;
            }
            int nconnects = cel->NConnects();
            for (int icon=0; icon<nconnects; icon++) {
                TPZConnect &connect = cel->Connect(icon);

                int lagrangemult = connect.LagrangeMultiplier();
                //Increment the number of connected elements for the avg pressure in order to not condense them
                if (lagrangemult==3) {
                    connect.IncrementElConnected();
                }
            }
        }
        GroupAndCondense(subcmesh);

     //   std::cout << "Number of connects (2) " << fCMesh->NEquations() << std::endl;
     //   TPZCompMeshTools::CreatedCondensedElements(subcmesh, false, true);

        subcmesh->CleanUpUnconnectedNodes();

        int numthreads = 4;
        int preconditioned = 0;
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZAutoPointer<TPZGuiInterface> guiInterface;

        int ninter_materials = subcmesh->MaterialVec().size();
//
        if(ninter_materials>3){
            subcmesh->SetAnalysisFStruct(numthreads);
        }else{
            subcmesh->SetAnalysisNonSymSparse(numthreads);
        }

        //subcmesh->SetAnalysisFStruct(numthreads);
//        subcmesh->SetAnalysisFStruct(numthreads);
        //subcmesh->SetAnalysisFStruct(numthreads);
//        subcmesh->SetAnalysisNonSymSparse(numthreads);
        //subcmesh->SetAnalysisFrontal(numthreads, guiInterface);
        //subcmesh->SetAnalysisNonSymSparse(numthreads);
        //subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);

        if(0) {
            std::ofstream filehide2("subcmeshAfter.txt");
            subcmesh->Print(filehide2);
        }

    }
    
}


void TPZMHMNavierStokesMeshControl::GroupAndCondense(TPZCompMesh *cmesh_m){
    
    //Criando agrupamento de elementos
    
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

    int nenvel = elgroups.NElements();
    for (int64_t ienv=0; ienv<nenvel; ienv++) {
        TPZElementGroup *elgr = elgroups[ienv];

        int nc = elgroups[ienv]->GetElGroup()[0]->NConnects();
    //    elgroups[ienv]->GetElGroup()[0]->Connect(nc-1).IncrementElConnected();
        if (fsetCoarseAverageMultipliers) {
        //        elgroups[ienv]->GetElGroup()[0]->Connect(nc-2).IncrementElConnected();
                int con_index = elgroups[ienv]->GetElGroup()[0]->ConnectIndex(nc-3);
                elgroups[ienv]->GetElGroup()[0]->Connect(nc-3).IncrementElConnected(); //pressão média
        }
        new TPZCondensedCompEl(elgr);
    }
    
    
    cmesh_m->CleanUpUnconnectedNodes();
    cmesh_m->ExpandSolution();
    
}

/// build the multi physics mesh (not at the finest geometric mesh level
void TPZMHMNavierStokesMeshControl::BuildMultiPhysicsMesh()
{
    if (fCMesh->NElements() != 0) {
        DebugStop();
    }
    fCMesh->SetAllCreateFunctionsMultiphysicElem();
    fCMesh->SetAllCreateFunctionsMultiphysicElemWithMem();

    TPZMultiphysicsCompMesh *mphysics = dynamic_cast<TPZMultiphysicsCompMesh *>(fCMesh.operator->());

    TPZManVector<TPZCompMesh *,6 > cmeshes(4);
    cmeshes[0] = fFluxMesh.operator->();
    cmeshes[1] = fPressureFineMesh.operator->();
    cmeshes[2] = fDistrFluxMesh.operator->();
    cmeshes[3] = fAveragePressMesh.operator->();
    if (fsetCoarseAverageMultipliers) {
        cmeshes.Resize(6);
        cmeshes[4] = fCMeshLagrange.operator->();
        cmeshes[5] = fCMeshConstantPressure.operator->();
    }

    TPZManVector<int64_t> shouldcreate(fGMesh->NElements(),0);
    std::set<int> matids;
    for (auto it : fCMesh->MaterialVec()) {
        matids.insert(it.first);
    }
    int64_t nel = fFluxMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        // this means that all geometric elements associated with flux elements will generate a computational element
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    nel = fPressureFineMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    nel = fDistrFluxMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fDistrFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }
    nel = fAveragePressMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fAveragePressMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (matids.find(gel->MaterialId()) == matids.end()) {
            DebugStop();
        }
        if (cel) {
            shouldcreate[cel->Reference()->Index()] = 1;
        }
    }

    if (fsetCoarseAverageMultipliers) {
        nel = fCMeshLagrange->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMeshLagrange->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (matids.find(gel->MaterialId()) == matids.end()) {
                DebugStop();
            }
            if (cel) {
                shouldcreate[cel->Reference()->Index()] = 1;
            }
        }

        nel = fCMeshConstantPressure->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fCMeshConstantPressure->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if (matids.find(gel->MaterialId()) == matids.end()) {
                DebugStop();
            }
            if (cel) {
                shouldcreate[cel->Reference()->Index()] = 1;
            }
        }
    }

    // define the intersection of the finest references
    nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (shouldcreate[el])
        {
            TPZGeoEl *fat = gel->Father();
            while(fat)
            {
                if(shouldcreate[fat->Index()] == 1)
                {
                    shouldcreate[fat->Index()] = 0;
                }
                fat = fat->Father();
            }
        }
    }
    TPZStack<int64_t> gelindexes;
    for (int64_t el=0; el<nel; el++) {
        if (shouldcreate[el])
        {
            gelindexes.Push(el);
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Geometric indices for which we will create multiphysics elements" << std::endl;
        sout << gelindexes;
//        std::cout << sout.str() << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
   // mphysics->BuildMultiphysicsSpace(cmeshes,gelindexes);

    TPZManVector<int,5> active_approx_spaces(6);
    active_approx_spaces[0] = 1;
    active_approx_spaces[1] = 1;
    active_approx_spaces[2] = 1;
    active_approx_spaces[3] = 1;
    active_approx_spaces[4] = 1;
    active_approx_spaces[5] = 1;

//    int dimension = fGMesh->Dimension();
    //mphysics->SetDimModel(dimension);
    //mphysics->BuildMultiphysicsSpace(cmeshes);
    TPZManVector<int> active(cmeshes.size(),1);
    mphysics->BuildMultiphysicsSpaceWithMemory(active,cmeshes);

}
