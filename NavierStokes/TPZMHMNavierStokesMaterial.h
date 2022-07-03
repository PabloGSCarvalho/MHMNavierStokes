/*
 *  TPZMHMNavierStokesMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMatWithMem.h"
#include "pzfmatrix.h"
#include "TPZBndCondT.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "TPZNavierStokesMaterial.h"
#include "pztrnsform.h"
#include "TPZNSMemory.h"
#include "TPZMatInterfaceCombinedSpaces.h"

#ifndef TPZMHMNavierStokesMATERIAL
#define TPZMHMNavierStokesMATERIAL

#include "fad.h"
#include "tfad.h"

class TPZMHMNavierStokesMaterial : public TPZNavierStokesMaterial  {
    
protected:
    
    STATE fMultiplier;
    
    
public:

    /**
     * Empty Constructor
     */
    TPZMHMNavierStokesMaterial() : TPZNavierStokesMaterial(), fMultiplier(1.)
    {
    }
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZMHMNavierStokesMaterial(int matid, int dimension) : TPZNavierStokesMaterial(matid,dimension), fMultiplier(1.)
    {

    }
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZMHMNavierStokesMaterial(const TPZMHMNavierStokesMaterial &mat) : TPZNavierStokesMaterial(mat)
    {}
    
    /**
     * Destructor
     */
    ~TPZMHMNavierStokesMaterial()
    {
    }
    
    TPZMHMNavierStokesMaterial &operator=(const TPZMHMNavierStokesMaterial &copy)
    {
        DebugStop();
        return *this;
    }
    
    virtual void SetMultiplier(STATE mult)
    {
        fMultiplier = mult;
    }


    
    virtual TPZMaterial *NewMaterial() const override
    {
        return new TPZMHMNavierStokesMaterial(*this);
    }
    
    TPZFMatrix<STATE> Transpose(TPZFMatrix<STATE> &GradU );

    TPZFNMatrix<9,REAL> Tangential(TPZManVector<REAL,3> normal);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                              REAL weight, TPZFMatrix<STATE> &ek,
                              TPZFMatrix<STATE> &ef,
                              TPZBndCondT<STATE> &bc) override;
    
    /**
     * @brief This method defines which parameters need to be initialized in order to compute the contribution
     * of interface elements
     */
    virtual void
    FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data,
                                  std::map<int, TPZMaterialDataT<STATE>> &datavec_left,
                                  std::map<int, TPZMaterialDataT<STATE>> &datavec_right) override;
    

    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                          std::map<int, TPZMaterialDataT<STATE>> &datavec,
                          REAL weight,
                          TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                                                   TPZBndCondT<STATE> &bc) override;

    /**
     * @brief Computes a contribution to the stiffness matrix and load vector at one integration point
     * to multiphysics simulation
     * @param [in] data
     * @param [in] dataleft
     * @param [in] dataright
     * @param [in] weight
     * @param [out] ek is the stiffness matrix
     * @param [out] ef is the load vector
     * @since June 5, 2012
     */
    virtual void ContributeInterface(const TPZMaterialDataT<STATE> &data,
                                     std::map<int, TPZMaterialDataT<STATE>> &dataleft,
                                     std::map<int, TPZMaterialDataT<STATE>> &dataright, REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef) override;


    TPZManVector<REAL,3> ComputeNormal(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright);
    
    /** @brief Creates an associated boundary condition.
     @param[in] reference The volumetric material associated with the BC.
     @param[in] id Boundary condition identifier.
     @param[in] type Type of the boundary condition.
     @param[in] val1 Value to be set at the element matrix.
     @param[in] val2 Value to be set at the rhs vector.
    */
    virtual TPZBndCondT<STATE>* CreateBC(TPZMaterial *reference,
                                        int id, int type,
                                        const TPZFMatrix<STATE> &val1,
                                        const TPZVec<STATE> &val2) override
    {
        return new  TPZBndCondBase<STATE,TPZMatCombinedSpacesBC<STATE>, TPZMatInterfaceCombinedSpacesBC<STATE> >
        (reference,id, type,val1,val2);
    }

    /**
     * @brief Returns the solution associated with the var index based on the finite element approximation around
     * one interface element
     */
    virtual void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                                   const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                                   const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                                   int var, TPZVec<STATE> &Solout) override {
        DebugStop();
    }

    /**
     * @brief Returns the solution associated with the var index based on the finite element approximation around
     * one interface element
     */
    virtual void SolutionInterface(const TPZMaterialDataT<STATE> &data,
                                   const std::map<int, TPZMaterialDataT<STATE>> &dataleftvec,
                                   const std::map<int, TPZMaterialDataT<STATE>> &datarightvec,
                                   int var, TPZVec<STATE> &Solout,
                                   TPZCompEl *left,TPZCompEl *right) override
    {
        DebugStop();
    }

};

#endif
