//
//  TPZSimulationData.h
//  PZ
//
//  Created by Omar on 8/28/18.
//
//

#ifndef TPZSimulationData_h
#define TPZSimulationData_h

#include <stdio.h>
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include <iostream>
#include <stdio.h>
#include <string>
#include "TPZGmshReader.h"
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"
#include "TPZAnalyticSolution.h"


/** @brief Object conatining several kind of informations being used anytime and anywhere */
class TPZSimulationData
{
    
protected:
    
    /** Number of coarse elements in each axis*/
    TPZVec<int> m_n_divs;
    
    /** Number of refinements in each internal element*/
    int m_n_intrefs;
    
    /** Size of domain*/
    TPZVec<REAL> m_h_domain;
    
    /** Polynomial order for internal elements */
    int m_internal_order;
    
    /** Polynomial order for internal elements */
    int m_skeleton_order;
    
    /** Physical dimension of the domain */
    int m_dimesion;
    
    /** Number of threads */
    int m_n_threads;

    /** Viscosity coeficient */
    REAL m_visco;

    /** Permeability coeficient */
    REAL m_permeability;

    /** Brinkman coeficient */
    REAL m_cBrinkman;

    /** Density of the fluid */
    REAL m_density;

    /** Coeficient for symmetry */
    REAL m_Sym_coef;

    /** Coeficient for approximation Space*/
    int m_space;

    /** Set the formulation of the problem*/
    TStokesAnalytic::MProblemType m_problem_type;

    /** Set the domain of the problem*/;
    TStokesAnalytic::EExactSol m_domain_type;

    /** Testing shape function */
    bool m_testshape;
    
    /** @brief Number of iteration */
    int m_n_iterations;
    
    /** @brief Residue overal tolerance */
    REAL m_epsilon_res;
    
    /** @brief Correction overal tolerance */
    REAL m_epsilon_cor;
    
    bool m_is_initial_state_Q;
    
    bool m_is_current_state_Q;
    
    bool m_optimizeBandwidth_Q;
    
    /** @brief Directive that states if the current solution must be accepted inside the memory  */
    bool m_must_accept_solution_Q;

    bool m_run_post_porcessing;

    bool m_is_pardiso_Q;

    bool m_is_static_condensation_Q;

    /** @brief Vector that storages only volumetric material identifiers (higher dimension elements) */
    std::vector<int> m_volumetric_material_id;

    bool m_is_transient_Q;

    REAL m_deltaT;

    REAL m_time;

    REAL m_timeTotal;

public:
    
    /** default constructor */
    TPZSimulationData();
    
    /** default constructor */
    TPZSimulationData(const TPZSimulationData & other);
    
    /** default constructor */
    TPZSimulationData &operator=(const TPZSimulationData &other);
    
    /** destructor */
    ~TPZSimulationData();
    
    /** Print object attributes */
    void Print();

    /** Set the number of coarse elements in each axis*/
    void SetCoarseDivisions(TPZVec<int> n_divs){
        m_n_divs = n_divs;
    }

    /** Get normal stiffness for each fracture*/
    TPZVec<int> GetCoarseDivisions(){
        return m_n_divs;
    }
    
    /** Set the number of refinements in each internal element*/
    void SetNInterRefs(int n_refs){
        m_n_intrefs = n_refs;
    }
    
    /** Get the number of refinements in each internal element*/
    int GetNInterRefs(){
        return m_n_intrefs;
    }
    
    /** Set the size of domain*/
    void SetDomainSize(TPZVec<REAL> h_domain){
        m_h_domain = h_domain;
    }

    /** Get the size of domain*/
    TPZVec<REAL> GetDomainSize(){
        return m_h_domain;
    }
    
    /** Set polynomial order for internal elements */
    void SetInternalOrder(int internal_order){
        m_internal_order = internal_order;
    }

    /** Get polynomial order for internal elements */
    int GetInternalOrder(){
        return m_internal_order;
    }
    
    /** Set polynomial order for skeleton elements */
    void SetSkeletonOrder(int skeleton_order){
        m_skeleton_order = skeleton_order;
    }
    
    /** Get polynomial order for skeleton elements */
    int GetSkeletonOrder(){
        return m_skeleton_order;
    }

    /** Set the physical dimension of the domain */
    void SetDimension(int dim){
        m_dimesion = dim;
    }
    
    /** Get the physical dimension of the domain */
    int GetDimension(){
        return m_dimesion;
    }

    /** Set the number of threads */
    void SetNthreads(int nthreads){
        m_n_threads = nthreads;
    }
    
    /** Get the number of threads */
    int GetNthreads(){
        return m_n_threads;
    }

    /** Set the viscosity coeficient */
    void SetViscosity(REAL viscosity){
        m_visco = viscosity;
    }
    
    /** Get the viscosity coeficient */
    REAL GetViscosity(){
        return m_visco;
    }

    /** Set the permeability coeficient */
    void SetPermeability(REAL perm){
        m_permeability = perm;
    }

    /** Get the permeability coeficient */
    REAL GetPermeability(){
        return m_permeability;
    }

    /** Set the Brinkman coeficient */
    void SetBrinkmanCoef(REAL coefB){
        m_cBrinkman = coefB;
    }

    /** Get the viscosity coeficient */
    REAL GetBrinkmanCoef(){
        return m_cBrinkman;
    }


    /** Set the density coeficient */
    void SetDensity(REAL density){
        m_density = density;
    }

    /** Get the density coeficient */
    REAL GetDensity(){
        return m_density;
    }

    /** Set the symmetry coeficient */
    void SetSymmetryCoef(REAL sym_coef){
        m_Sym_coef = sym_coef;
    }

    /** Get the symmetry coeficient */
    REAL GetSymmetryCoef(){
        return m_Sym_coef;
    }

    /** Set the space coeficient */
    void SetSpaceCoef(int space_coef){
        m_space = space_coef;
    }

    /** Get the symmetry coeficient */
    int GetSpaceCoef(){
        return m_space;
    }

    void SetProblemType(TStokesAnalytic::MProblemType type){
        m_problem_type = type;
    }

    TStokesAnalytic::MProblemType GetProblemType(){
        return m_problem_type;
    }


    void SetDomainType(TStokesAnalytic::EExactSol type){
        m_domain_type = type;
    }

    TStokesAnalytic::EExactSol GetDomainType(){
        return m_domain_type;
    }

    /** Set shape test true */
    void SetShapeTest(){
        m_testshape = true;
    }
    
    /** Get shape test true */
    bool GetShapeTest(){
        return m_testshape;
    }
    
    /** @brief Set Newton iterations */
    void Set_n_iterations(int n_iterations){
        m_n_iterations = n_iterations;
    }
    
    /** @brief Get Newton iterations  */
    REAL Get_n_iterations(){
        return m_n_iterations;
    }
    
    /** @brief Set residue tolerance */
    void Set_epsilon_res(REAL epsilon_res){
        m_epsilon_res = epsilon_res;
    }
    
    /** @brief Get residue tolerance  */
    REAL Get_epsilon_res(){
        return m_epsilon_res;
    }
    
    /** @brief Set correction tolerance */
    void Set_epsilon_cor(REAL epsilon_cor){
        m_epsilon_cor = epsilon_cor;
    }
    
    /** @brief Get correction tolerance  */
    REAL Get_epsilon_cor(){
        return m_epsilon_cor;
    }
    
    /** @brief Set initial state */
    void SetInitialStateQ(bool state) {
        m_is_initial_state_Q = state;
    }
    
    /** @brief Get initial state */
    bool IsInitialStateQ() {
        return m_is_initial_state_Q;
    }

    /** @brief Set Pardiso option */
    void SetPardisoSolver() {
        m_is_pardiso_Q = true;
    }

    /** @brief Get Pardiso option */
    bool IsPardisoSolverQ() {
        return m_is_pardiso_Q;
    }

    /** @brief Set Pardiso option */
    void SetStaticCondensation(bool static_condense) {
        m_is_static_condensation_Q = static_condense;
    }

    /** @brief Get Pardiso option */
    bool IsStaticCondensedQ() {
        return m_is_static_condensation_Q;
    }

    /** @brief Set current time state */
    void SetCurrentStateQ(bool state) { m_is_current_state_Q = state; }
    
    /** @brief Get current time state */
    bool IsCurrentStateQ() {return m_is_current_state_Q;}

    /** @brief Set current time state */
    void ActivatePostProcessing() { m_run_post_porcessing = true; }

    bool IsPostProcessingActivatedQ() {return m_run_post_porcessing;}


    /** @brief Set the directive that states if the current solution must be accepted inside the memory  */
    void Set_must_accept_solution_Q(bool must_accept_solution_Q){
        m_must_accept_solution_Q = must_accept_solution_Q;
    }
    
    /** @brief Get the directive that states if the current solution must be accepted inside the memory  */
    bool Get_must_accept_solution_Q() { return m_must_accept_solution_Q; }
    
    /** @brief Get the vector that storages only volumetric material identifiers (higher dimension elements)  */
    std::vector<int> & Get_volumetric_material_id(){
        return m_volumetric_material_id;
    }
    
    /** @brief Set the vector that storages only volumetric material identifiers (higher dimension elements) */
    void Set_volumetric_material_id(std::vector<int> & volumetric_material_id){
        m_volumetric_material_id = volumetric_material_id;
    }

    void SetOptimizeBandwidthQ(bool set_optimze){
        m_optimizeBandwidth_Q = set_optimze;
    }
    
    bool GetOptimizeBandwidthQ(){
        return m_optimizeBandwidth_Q;
    }

    void SetTimeStep(REAL deltaT){
        m_deltaT = deltaT;
    }

    REAL GetTimeStep(){
        return m_deltaT;
    }

    void SetTransientQ(bool is_transient_Q){
        m_is_transient_Q = is_transient_Q;
    }

    bool IsTransientQ(){
        return m_is_transient_Q;
    }

    void SetTime(REAL time){
        m_time = time;
    }

    REAL GetTime(){
        return m_time;
    }

    void SetTimeTotal(REAL timeTotal){
        m_timeTotal = timeTotal;
        if(m_timeTotal>0){
            SetTransientQ(true);
        }
    }

    REAL GetTimeTotal(){
        return m_timeTotal;
    }

    void UpdateTime();

};

#endif /* TPZSimulationData_h */
