//
//  TPZNSBuildTransfers.h
//  PZ
//
//  Created by Pablo Carvalho on 08/31/20.
//  This class storage approximation space in global integration points arrays for volumetric and boundary elements
//

#ifndef TRMBuildTransfers_h
#define TRMBuildTransfers_h

#include <stdio.h>

#include "tpzintpoints.h"
#include "TPZMatWithMem.h"
#include "TPZNSMemory.h"
#include "TPZMHMNavierStokesMaterial.h"
#include "TPZMaterial.h"
#include "pzinterpolationspace.h"
#include "pzmultiphysicselement.h"
#include "pzcondensedcompel.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"

#include "pzysmp.h"
#include "pzblockdiag.h"
#include "TPZSimulationData.h"


class TPZNSBuildTransfers{
    
    
private:
    
    
    
    /**
     * @defgroup Attributes
     * @{
     */
    
    /** @brief Autopointer of simulation data */
    TPZSimulationData * fSimulationData;
    
    
    /**
     * @defgroup Sparses Matrices to transfer information to the mixed problem on Omega
     * @{
     */
    
    /** @brief Diagonal block matrix to transfer u flux solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<STATE> fu_To_Mixed;
    
    /** @brief flux dof indexes per element */
    TPZVec< TPZVec<int64_t> > fu_dof_scatter;
    
    /** @brief Diagonal block matrix to transfer Pressure solution to integrations points of the mixed mesh */
    TPZBlockDiagonal<STATE> fp_To_Mixed;
    
    /** @brief pressure dof indexes per element */
    TPZVec< TPZVec<int64_t> > fp_dof_scatter;


public:
    
    /** @brief Default constructor */
    TPZNSBuildTransfers();
    
    /** @brief Default desconstructor */
    ~TPZNSBuildTransfers();
    
    /** @brief Copy constructor $ */
    TPZNSBuildTransfers(const TPZNSBuildTransfers &copy);
    
    /** @brief Copy assignemnt operator $ */
    TPZNSBuildTransfers &operator=(const TPZNSBuildTransfers &other);
    

    /**
     * @defgroup Apply transfers to different meshes
     * @{
     */

    /** @brief Transfer Flux to integration points of multiphysics mesh over volumetric elements */
    void u_To_Mixed_Memory(TPZCompMesh * cmesh_flux, TPZCompMesh * cmesh_multiphysics);
    
    /** @brief Transfer Pressure to integration points of multiphysics mesh over volumetric elements */
    void p_To_Mixed_Memory(TPZCompMesh * cmesh_pressure, TPZCompMesh * cmesh_multiphysics);


    /**
     * @defgroup Create, compute and get transfer matrices
     * @{
     */

    /** @brief Initializate diagonal block matrix to transfer flux to multiphysics mesh  */
    void Initialize_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Initializate diagonal block matrix to transfer flux to multiphysics mesh  */
    void Fill_u_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TPZBlockDiagonal<STATE> Transfer_u_To_Mixed(){
        return fu_To_Mixed;
    }
    
    /** @brief Initializate  diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Initialize_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Initializate diagonal block matrix to transfer Pressure to multiphysics mesh  */
    void Fill_p_To_Mixed(TPZCompMesh * cmesh_multiphysics, int mesh_index);
    
    /** @brief Get the sparse matrix to transfer Pressure to multiphysics mesh  */
    TPZBlockDiagonal<STATE> Transfer_p_To_Mixed(){
        return fp_To_Mixed;
    }
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Memory operations
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Get Global integration point indexes associaded  */
    void GlobalPointIndexes(TPZCompEl * cel, TPZManVector<int64_t,30> &int_point_indexes);

    
    /** @brief Compute element dof indexes */
    void ElementDofIndexes(TPZInterpolationSpace * &intel,  TPZVec<int64_t> &dof_indexes);


    /** @brief Set autopointer of Simulation data */
    void SetSimulationData(TPZSimulationData * SimulationData){
        fSimulationData = SimulationData;
    }
    
    /** @brief Get autopointer of Simulation data */
    TPZSimulationData * SimulationData(){
        return fSimulationData;
    }
    
};


#endif /* TPZNSBuildTransfers_h */

