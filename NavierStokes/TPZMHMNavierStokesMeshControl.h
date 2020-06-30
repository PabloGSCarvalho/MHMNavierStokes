//
//  TPZMHMStokesMeshControl.hpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#ifndef TPZMHMNavierStokesMeshControl_hpp
#define TPZMHMNavierStokesMeshControl_hpp

#include <stdio.h>

#include "TPZMHMixedMeshControl.h"
#include "TPZMultiphysicsCompMesh.h"

/// class for creating TPZMHMM with Stokes problem Meshes
class TPZMHMNavierStokesMeshControl : public TPZMHMixedMeshControl
{
    
protected:
    
    /// Computational mesh to contain the avarege pressure elements
    TPZAutoPointer<TPZCompMesh> fAveragePressMesh;

    /// Computational mesh to contain the distributed flux elements
    TPZAutoPointer<TPZCompMesh> fDistrFluxMesh;

    /// Computational mesh to contain the coarse avarege pressure elements
    TPZAutoPointer<TPZCompMesh> fCoarseAveragePressMesh;
    
    /// Computational mesh to contain the coarse distributed flux elements
    TPZAutoPointer<TPZCompMesh> fCoarseDistrFluxMesh;
    
    /// Set coarse average pressure and distributed flux
    bool fsetCoarseAverageMultipliers;

    bool fsetStaticCondensation;


public:
    
    TPZMHMNavierStokesMeshControl() : TPZMHMixedMeshControl()
    {
        
    }
    
    TPZMHMNavierStokesMeshControl(int dimension);
    
    TPZMHMNavierStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices);
    
    TPZMHMNavierStokesMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh);
    
    TPZMHMNavierStokesMeshControl(const TPZMHMNavierStokesMeshControl &copy);
    
    TPZMHMNavierStokesMeshControl &operator=(const TPZMHMNavierStokesMeshControl &cp);
    
    virtual ~TPZMHMNavierStokesMeshControl();
    
    /// Create all data structures for the computational mesh
    virtual void BuildComputationalMesh(bool usersubstructure);
    
    void InsertInternalSkeleton();
    
    void InsertBCSkeleton();
    
    /// Put the pointers to the meshes in a vector
    void GetMeshVec(TPZVec<TPZCompMesh *> &meshvec)
    {
        meshvec.Resize(4);
        meshvec[0] = fFluxMesh.operator->();
        meshvec[1] = fPressureFineMesh.operator->();
        meshvec[2] = fDistrFluxMesh.operator->();
        meshvec[3] = fAveragePressMesh.operator->();

        if (fsetCoarseAverageMultipliers) {
            meshvec.Resize(6);
            meshvec[4] = fCMeshLagrange.operator->();
            meshvec[5] = fCMeshConstantPressure.operator->();
        }

    }
    
    TPZVec<TPZAutoPointer<TPZCompMesh> > GetMeshes()
    {
        TPZManVector<TPZAutoPointer<TPZCompMesh>,6> result(4);
        result[0] = fFluxMesh;
        result[1] = fPressureFineMesh;
        result[2] = fDistrFluxMesh;
        result[3] = fAveragePressMesh;

        if (fsetCoarseAverageMultipliers) {
            result.Resize(6);
            result[4] = fCMeshLagrange;
            result[5] = fCMeshConstantPressure;
        }
        
        return result;
    }
    
    /// Set the flag for creating Lagrange Dofs for the average pressure and distributed flux on coarse mesh
    void SetCoarseAverageMultipliers(bool flag)
    {
        fsetCoarseAverageMultipliers = flag;
    }

    void SetStaticCondensation(bool flag)
    {
        fsetStaticCondensation = flag;
    }
    
protected:

    
    /// Insert the necessary pressure material objects to create the pressure mesh
    virtual void InsertPeriferalPressureMaterialObjects();
    
    /// Create the pressure mesh which contain traction variable
    void CreatePressureAndTractionMHMMesh();

    /// Insert the necessary average pressure material objects to create the average pressure mesh
    void InsertPeriferalAveragePressMaterialObjects();

    /// Create the average pressure mesh
    void CreateAveragePressMHMMesh();

    /// Insert the necessary coarse average pressure material objects to create the average pressure mesh
    void InsertPeriferalCoarseAveragePressMaterialObjects();
    
    /// Create the coarse average pressure mesh
    void CreateCoarseAveragePressMHMMesh();
    
    /// Insert the necessary distributed flux material objects to create the distributed flux material pressure mesh
    void InsertDistributedFluxMaterialObjects();
    
    /// Create the distributed flux mesh
    void CreateDistributedFluxMHMMesh();
    
    /// Insert the necessary distributed flux material objects to create the distributed flux material pressure mesh
    void InsertCoarseDistributedFluxMaterialObjects();
    
    /// Create the distributed flux mesh
    void CreateCoarseDistributedFluxMHMMesh();
    
    /// Create multiphysics MHM mesh
    void CreateMultiPhysicsMHMMesh();
    
//    /// build the multi physics mesh (not at the finest geometric mesh level)
//    virtual void BuildMultiPhysicsMesh();

    /// Create the multiphysics interface elements between elements of traction material id
    void CreateMultiPhysicsInterfaceElements();

    /// Create the multiphysics BC interface elements between elements of BC traction material id
    void CreateMultiPhysicsBCInterfaceElements();

    void GroupandCondenseSubMeshes();
    
    void GroupAndCondense(TPZCompMesh *cmesh_m);

    void BuildSubMeshes();

    void PutinSubmeshes(TPZCompMesh *cmesh, std::map<int64_t,std::set<int64_t> >&elindices, std::map<int64_t,int64_t> &indices, int KeepOneLagrangian);

    void BuildMultiPhysicsMesh();

public:
    
    /// material id associated with the internal skeleton elements
    int64_t fTractionMatId = 3;

    /// material ids associated with the BC skeleton elements
    std::map<int64_t,int64_t> fBCTractionMatIds;
    
};

#endif /* TPZMHMixedMeshControl_hpp */


