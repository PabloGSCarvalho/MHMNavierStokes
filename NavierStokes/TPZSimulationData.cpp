//
//  TPZSimulationData.cpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#include "TPZSimulationData.h"

TPZSimulationData::TPZSimulationData()
{
    m_n_divs.resize(0);
    m_h_domain.resize(0);
    m_internal_order = 0;
    m_skeleton_order = 0;
    m_dimesion = 0;
    m_n_threads = 0;
    m_visco = 0;
    m_density = 0.;
    m_Sym_coef = 0;
    m_space = 1;
    m_problem_type = TStokesAnalytic::EStokes;
    m_domain_type = TStokesAnalytic::ENone;
    m_n_intrefs = 0.;
    m_testshape = false;
    m_n_iterations = 0;
    m_epsilon_res = 0;
    m_epsilon_cor = 0;
    m_is_initial_state_Q  = false;
    m_is_current_state_Q  = false;
    m_must_accept_solution_Q = false;
    m_volumetric_material_id.resize(1);
    m_volumetric_material_id[0]=1;
    m_optimizeBandwidth_Q = false;
    m_run_post_porcessing = false;
    m_is_pardiso_Q = false;
    m_is_static_condensation_Q = true;
    m_is_transient_Q = false;
    m_deltaT = 0.;
    m_time = 0.;
    m_timeTotal = 0.;

}

TPZSimulationData::~TPZSimulationData()
{
    
}

TPZSimulationData::TPZSimulationData(const TPZSimulationData & other)
{
    m_n_divs                           = other.m_n_divs;
    m_h_domain                         = other.m_h_domain;
    m_internal_order                   = other.m_internal_order;
    m_skeleton_order                   = other.m_skeleton_order;
    m_dimesion                         = other.m_dimesion;
    m_n_threads                        = other.m_n_threads;
    m_visco                            = other.m_visco;
    m_density                          = other.m_density;
    m_Sym_coef                         = other.m_Sym_coef;
    m_space                            = other.m_space;
    m_problem_type                     = other.m_problem_type;
    m_domain_type                      = other.m_domain_type;
    m_n_intrefs                        = other.m_n_intrefs;
    m_testshape                        = other.m_testshape;
    m_n_iterations                      = other.m_n_iterations;
    m_epsilon_res                       = other.m_epsilon_res;
    m_epsilon_cor                       = other.m_epsilon_cor;
    m_is_initial_state_Q                = other.m_is_initial_state_Q;
    m_is_current_state_Q                = other.m_is_current_state_Q;
    m_must_accept_solution_Q            = other.m_must_accept_solution_Q;
    m_volumetric_material_id            = other.m_volumetric_material_id;
    m_optimizeBandwidth_Q               = other.m_optimizeBandwidth_Q;
    m_run_post_porcessing               = other.m_run_post_porcessing;
    m_is_pardiso_Q                      = other.m_is_pardiso_Q;
    m_is_static_condensation_Q          = other.m_is_static_condensation_Q;
    m_is_transient_Q                    = other.m_is_transient_Q;
    m_deltaT                            = other.m_deltaT;
    m_time                              = other.m_time;
    m_timeTotal                         = other.m_timeTotal;

}

TPZSimulationData & TPZSimulationData::operator=(const TPZSimulationData &other)
{
    if (this != & other) // prevent self-assignment
    {
        m_n_divs                           = other.m_n_divs;
        m_h_domain                         = other.m_h_domain;
        m_internal_order                   = other.m_internal_order;
        m_skeleton_order                   = other.m_skeleton_order;
        m_dimesion                         = other.m_dimesion;
        m_n_threads                        = other.m_n_threads;
        m_visco                            = other.m_visco;
        m_density                          = other.m_density;
        m_Sym_coef                         = other.m_Sym_coef;
        m_space                            = other.m_space;
        m_problem_type                     = other.m_problem_type;
        m_domain_type                      = other.m_domain_type;
        m_n_intrefs                        = other.m_n_intrefs;
        m_testshape                        = other.m_testshape;
        m_n_iterations                      = other.m_n_iterations;
        m_epsilon_res                       = other.m_epsilon_res;
        m_epsilon_cor                       = other.m_epsilon_cor;
        m_is_initial_state_Q                = other.m_is_initial_state_Q;
        m_is_current_state_Q                = other.m_is_current_state_Q;
        m_must_accept_solution_Q            = other.m_must_accept_solution_Q;
        m_volumetric_material_id            = other.m_volumetric_material_id;
        m_optimizeBandwidth_Q               = other.m_optimizeBandwidth_Q;
        m_run_post_porcessing               = other.m_run_post_porcessing;
        m_is_pardiso_Q                      = other.m_is_pardiso_Q;
        m_is_static_condensation_Q          = other.m_is_static_condensation_Q;
        m_is_transient_Q                    = other.m_is_transient_Q;
        m_deltaT                            = other.m_deltaT;
        m_time                              = other.m_time;
        m_timeTotal                         = other.m_timeTotal;
    }
    return *this;
}

void TPZSimulationData::UpdateTime()
{

    if(m_time + m_deltaT > m_timeTotal)
    {
        m_deltaT = m_timeTotal - m_time;
    }
    m_time += m_deltaT;

    std::cout << "\n----------------- Actual Time of Simulation = " << m_time << " -----------------" << std::endl;

}

void TPZSimulationData::Print()
{
    
    std::cout << " TPZSimulationData class members : " << std::endl;
    std::cout << std::endl;
    std::cout << " m_n_divs = " << m_n_divs << std::endl;
    std::cout << " m_h_domain = " << m_h_domain << std::endl;
    std::cout << " m_internal_order = " << m_internal_order << std::endl;
    std::cout << " m_skeleton_order = " << m_skeleton_order << std::endl;
    std::cout << " m_dimesion = " << m_dimesion << std::endl;
    std::cout << " m_n_threads = " << m_n_threads << std::endl;
    std::cout << " m_visco = " << m_visco << std::endl;
    std::cout << " m_density = " << m_density << std::endl;
    std::cout << " m_n_intrefs = " << m_n_intrefs << std::endl;
    std::cout << " m_testshape = " << m_testshape << std::endl;
    std::cout << " m_n_iterations = " << m_n_iterations << std::endl;
    std::cout << " m_epsilon_res = " << m_epsilon_res << std::endl;
    std::cout << " m_epsilon_cor = " << m_epsilon_cor << std::endl;
    std::cout << " m_volumetric_material_id = " << &m_volumetric_material_id << std::endl;
    std::cout << " m_deltaT = " << &m_deltaT << std::endl;
    std::cout << std::endl;
    
}

