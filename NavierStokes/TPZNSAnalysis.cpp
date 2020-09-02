//
//  TPZNSAnalysis.cpp
//  PMRS
//
//  Created by Omar Durán on 9/13/18.
//

#include "TPZNSAnalysis.h"
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("NavierStokes.Analysis"));
#endif

TPZNSAnalysis::TPZNSAnalysis() : TPZAnalysis() {
    m_simulation_data = NULL;
    m_U_Plus.Resize(0, 0);
    m_U_n.Resize(0, 0);
    m_mesh_vec.Resize(0);
    m_res_error = 0;
    m_dU_norm = 0;
    m_k_iterations = 0;
    m_post_processor = NULL;
    m_var_names.resize(0);
    m_vec_var_names.resize(0);
    m_R_Plus.Resize(0,0);
    m_R_n.Resize(0,0);
    m_LastStepRhs.Resize(0,0);

}

TPZNSAnalysis::~TPZNSAnalysis(){
    
}

TPZNSAnalysis::TPZNSAnalysis(const TPZNSAnalysis & other){
    m_simulation_data   = other.m_simulation_data;
    m_U_Plus            = other.m_U_Plus;
    m_U_n               = other.m_U_n;
    m_mesh_vec          = other.m_mesh_vec;
    m_res_error         = other.m_res_error;
    m_dU_norm           = other.m_dU_norm;
    m_k_iterations      = other.m_k_iterations;
    m_post_processor    = other.m_post_processor;
    m_var_names         = other.m_var_names;
    m_vec_var_names     = other.m_vec_var_names;
    m_R_Plus            = other.m_R_Plus;
    m_R_n               = other.m_R_n;
    m_LastStepRhs       = other.m_LastStepRhs;
}


void TPZNSAnalysis::ConfigureAnalysis(DecomposeType decomposition, TPZSimulationData * simulation_data, TPZCompMesh * cmesh_M, TPZVec<TPZCompMesh *> & mesh_vec, TPZVec<std::string> & var_names){
    
    SetSimulationData(simulation_data);
    bool mustOptimizeBandwidth = simulation_data->GetOptimizeBandwidthQ();
    this->SetCompMesh(cmesh_M,mustOptimizeBandwidth);
    TPZStepSolver<STATE> step;
    unsigned int n_threads = m_simulation_data->GetNthreads();
    
    if(!Mesh()){
        std::cout << "Call SetCompMesh method." << std::endl;
        DebugStop();
    }
    
    m_mesh_vec = mesh_vec;
    switch (decomposition) {
        case ELU:
        {
            if(m_simulation_data->IsPardisoSolverQ() == true){
                TPZSpStructMatrix struct_mat(Mesh());
                struct_mat.SetNumThreads(n_threads);
                this->SetStructuralMatrix(struct_mat);
            }else{
                TPZFStructMatrix struct_mat(Mesh());
                //TPZSkylineNSymStructMatrix struct_mat(Mesh());
                struct_mat.SetNumThreads(n_threads);
                this->SetStructuralMatrix(struct_mat);
            }
        }
            break;
        case ELDLt:
        {
            if(m_simulation_data->IsPardisoSolverQ() == true){
                TPZSymetricSpStructMatrix struct_mat(Mesh());
                struct_mat.SetNumThreads(n_threads);
                this->SetStructuralMatrix(struct_mat);
            }else{
                TPZFStructMatrix struct_mat(Mesh());
                struct_mat.SetNumThreads(n_threads);
                this->SetStructuralMatrix(struct_mat);
            }
        }
            break;
        default:
        {
            DebugStop();
        }
            break;
    }
    step.SetDirect(decomposition);
    this->SetSolver(step);
    this->Solution().Resize(Mesh()->NEquations(), 1);
    m_U_n.Resize(Mesh()->NEquations(), 1);
    m_U_Plus.Resize(Mesh()->Solution().Rows(), 1);
    
    m_post_processor = new TPZPostProcAnalysis;
    m_post_processor->SetCompMesh(Mesh());
    
    int n_vols = m_simulation_data->Get_volumetric_material_id().size();
    TPZManVector<int,10> post_mat_id(n_vols);
    for (int ivol = 0; ivol < n_vols; ivol++)
    {
        int matid = m_simulation_data->Get_volumetric_material_id()[ivol];
        post_mat_id[ivol] = matid;
    }
    
    for (auto i : var_names) {
        m_var_names.Push(i);
    }
    
    m_post_processor->SetPostProcessVariables(post_mat_id, m_var_names);
    int dim = Mesh()->Dimension();
    int div = 0;
    TPZStack< std::string> vecnames;
    std::string plotfile("NavierStokesPostProcess.vtk");
    
    m_post_processor->DefineGraphMesh(dim,m_var_names,vecnames,plotfile);
    
    TPZFStructMatrix structmatrix(m_post_processor->Mesh());
    structmatrix.SetNumThreads(n_threads);
    m_post_processor->SetStructuralMatrix(structmatrix);
    
}


void TPZNSAnalysis::ExecuteOneTimeStep(){

    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
    if (!cmesh) {
        DebugStop();
    }

    // m_X means the solution at the previous time step
    if (m_simulation_data->IsInitialStateQ()) {
        m_U_n = Solution();
    }

    //    // Initial guess
    //    m_X_n = m_X;
    // Set current state true means overwriting p_n of the memory object
    //m_simulation_data->SetCurrentStateQ(true);

    
    std::ofstream plotNavierEK("NavierStiffness.txt");
    std::ofstream plotNavierEF("NavierRhs.txt");
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled() && 0)
    {
        std::stringstream sout;
        fRhs.Print("Rhs =",sout);
        //EFormatted, EInputFormat, EMathematicaInput, EMatlabNonZeros, EMatrixMarket
        //  fSolver->Matrix()->Print("ek = ",plotDarcyEK,EMathematicaInput);
        fRhs.Print("ef = ",plotNavierEF,EMathematicaInput);

        PrintVectorByElement(sout, fRhs);
        PrintVectorByElement(sout, fSolution);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    TPZFMatrix<STATE> dU, U_n;
    bool residual_stop_criterion_Q = false;
    bool correction_stop_criterion_Q = false;
    REAL norm_res, norm_dU;
    REAL res_norm = m_simulation_data->Get_epsilon_res();
    REAL dU_norm = m_simulation_data->Get_epsilon_cor();
    int n_it = m_simulation_data->Get_n_iterations();
    //m_U_Plus = fCompMesh->Solution();
    //m_U_Plus.Redim(Solution().Rows(),Solution().Cols());
    if(0)
    {
        AssembleResidual();
        std::cout<< "residual norm 0 = " << Norm(this->Rhs()) <<std::endl;
    }

//    this->AssembleResidual();
//    m_R_n = this->Rhs();

    for (int i = 1; i <= n_it; i++) {

        TPZMHMNavierStokesMaterial *mat = dynamic_cast<TPZMHMNavierStokesMaterial *>(fCompMesh->FindMaterial(1));

//        if(i==1){
//            mat->SetProblemType(TStokesAnalytic::EOseenCDG);
//        }else{
//            mat->SetProblemType(TStokesAnalytic::ENavierStokesCDG);
//        }

        this->ExecuteNewtonIteration();

        dU = Solution();

      //  m_U_Plus += fCompMesh->Solution();
      //  m_U_Plus += dU;

  //      std::cout<<this->Rhs()<<std::endl;

        cmesh->UpdatePreviousState(1.);
        //cmesh->LoadSolutionFromMultiPhysics();
        Mesh()->TransferMultiphysicsSolution();

        AssembleResidual();

        norm_dU  = Norm(dU);
        m_R_Plus = this->Rhs();

        //std::cout<<this->Rhs()<<std::endl;

        m_res_error =  Norm(m_R_Plus); // residue error

        
        
#ifdef LOG4CXX
        if(0 && logger->isDebugEnabled())
        {
            std::stringstream sout;
            fRhs.Print("Rhs =",sout);
            {
                std::ofstream plotNavierStiff("NavierStiffness.txt");
                std::ofstream plotNavierRhs("NavierRhs.txt");
                fSolver->Matrix()->Print("ek = ",plotNavierStiff,EMathematicaInput);
                fRhs.Print("ef = ",plotNavierRhs,EMathematicaInput);
            }
            PrintVectorByElement(sout, fRhs);
            PrintVectorByElement(sout, fSolution);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        if(0)
        {
            
            std::ofstream plotNavierStiff("NavierStiffness_iteration.txt");
            fSolver->Matrix()->Print("ek = ",plotNavierStiff,EMathematicaInput);
            std::ofstream plotNavierRhs("NavierRhs_iteration.txt");;
            fRhs.Print("Rhs =",plotNavierRhs,EMathematicaInput);
        }

        norm_res = Norm(Rhs());
//        if(m_simulation_data->IsTransientQ()){
//            m_R_Plus += m_LastStepRhs;
//            norm_res = Norm(m_R_Plus);
//        }

        std::cout << "Correction norm 1 = " << norm_dU << std::endl;
        std::cout<< "residual norm 1 = " << norm_res <<std::endl;

        residual_stop_criterion_Q   = norm_res < res_norm;
        correction_stop_criterion_Q = norm_dU  < dU_norm;
        
        m_k_iterations = i;
        m_res_error = norm_res;
        m_dU_norm = norm_dU;
        
        
        if (residual_stop_criterion_Q && correction_stop_criterion_Q) {
#ifdef PZDEBUG
            std::cout << "TPZNSAnalysis:: Nonlinear process converged with residue norm = " << norm_res << std::endl;
            std::cout << "TPZNSAnalysis:: Number of iterations = " << i << std::endl;
            std::cout << "TPZNSAnalysis:: Correction norm = " << norm_dU << std::endl;
#endif
            //m_simulation_data->SetCurrentStateQ(true);

            break;
        }
    }
    
    if (residual_stop_criterion_Q == false) {
        std::cout << "TPZNSAnalysis:: Nonlinear process not converged with residue norm = " << norm_res << std::endl;
    }
    

}

void TPZNSAnalysis::PostProcessTimeStep(std::string & res_file){

  //  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, this->Mesh());

//    TPZMultiphysicsCompMesh * cmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
//    if (!cmesh) {
//        DebugStop();
//    }
    //LoadSolution();
    //cmesh->LoadSolutionFromMeshes();

    const int dim = this->Mesh()->Dimension();
    int div = 1;
    TPZStack<std::string> scalnames, vecnames;
//    if (fSimulationData->IsInitialStateQ()) {
//        plotfile =  "DualSegregatedDarcyOnBox_I.vtk";
//        return;
//    }
//    else{
//        plotfile =  "DualSegregatedDarcyOnBox.vtk";
//    }
    
    //Pós-processamento (paraview):

    scalnames.Push("P");
    vecnames.Push("V");
    vecnames.Push("f");
    vecnames.Push("V_exact");
    scalnames.Push("P_exact");
    scalnames.Push("Div");

    TPZMHMNavierStokesMaterial *mat = dynamic_cast<TPZMHMNavierStokesMaterial *>(fCompMesh->FindMaterial(1));
    if(mat->GetProblemType()==TStokesAnalytic::EOseenCDG||mat->GetProblemType()==TStokesAnalytic::ENavierStokesCDG){
        scalnames.Push("P_CDG");
        scalnames.Push("P_exact_CDG");
    }


    this->DefineGraphMesh(dim, scalnames, vecnames, res_file);
    //std::cout<<this->Mesh()->Solution()<<std::endl;
    REAL time = m_simulation_data->GetTime();
    this->SetTime(time);
    this->PostProcess(div,dim);
    

}

void TPZNSAnalysis::SolveSystem(){

    std::string file_NavierStokes("NavierStokes.vtk");
    //Testes
    std::string file_NavierStokes_test("NavierStokes_test.vtk");

    int n_max_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    int n_enforced_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    REAL res_norm = m_simulation_data->Get_epsilon_res();
    REAL dU_norm = m_simulation_data->Get_epsilon_cor();
    bool error_stop_criterion_Q = false;
    bool dU_stop_criterion_Q = false;
    this->SetInitialParameters();

    {
   //     std::ofstream filecE("MalhaC_E_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
   //     m_elastoplast_analysis->Mesh()->Print(filecE);
    }
    if(0)
    {
        std::ofstream filecM("MalhaC_M_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
        Mesh()->Print(filecM);
    }
    if(0)
    {
        int64_t neq = Mesh()->Solution().Rows();
        for (int is = 0; is<neq; is++) {
            Mesh()->Solution()(is,0) = 1.;
        }
        Mesh()->LoadSolution((Mesh()->Solution()));
        TPZMultiphysicsCompMesh *mphys = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
        mphys->LoadSolutionFromMultiPhysics();
    }


    for (int k = 1; k <= n_max_fss_iterations; k++) {
        this->ExecuteOneTimeStep();

        REAL rhsnomr = this->Get_error();
        REAL dunorm = this->Get_dU_norm();
        error_stop_criterion_Q = (rhsnomr < res_norm);
        dU_stop_criterion_Q = (dunorm < dU_norm);
        //this->PostProcessTimeStep(file_NavierStokes_test);

        if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) && dU_stop_criterion_Q) {
            //this->PostProcessTimeStep(file_NavierStokes);
            std::cout << "TPZNSAnalysis:: Iterative process converged with residue norm  = " << this->Get_error() << std::endl;
            UpdateState();
            break;
        }
        bool postprocessingQ = m_simulation_data->IsPostProcessingActivatedQ();
        if (error_stop_criterion_Q && dU_stop_criterion_Q && postprocessingQ) {
            std::cout << "Running post processing" << std::endl;
            this->PostProcessTimeStep(file_NavierStokes);
        }

//      if (postprocessingQ) {
//           std::cout << "Running post processing" << std::endl;
//                this->PostProcessTimeStep(file_NavierStokes);
//      }

    }

}

void TPZNSAnalysis::SolveSystemTransient(){

    std::string file_NavierStokes("NavierStokes.vtk");
    //Testes
    std::string file_NavierStokes_test("NavierStokes_test.vtk");

    int n_max_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    int n_enforced_fss_iterations = 1; // @TODO:: MS, please to xml file structure
    REAL time_total = m_simulation_data->GetTimeTotal();
    REAL time = 0.;
    REAL deltaT = m_simulation_data->GetTimeStep();
    int n_time_steps = time_total/deltaT;
    REAL res_norm = m_simulation_data->Get_epsilon_res();
    REAL dU_norm = m_simulation_data->Get_epsilon_cor();

    bool error_stop_criterion_Q = false;
    bool dU_stop_criterion_Q = false;
    this->SetInitialParameters();

    {
        //     std::ofstream filecE("MalhaC_E_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
        //     m_elastoplast_analysis->Mesh()->Print(filecE);
    }
    if(0)
    {
        std::ofstream filecM("MalhaC_M_AfterAdjust.txt"); //Impressão da malha computacional da velocidade (formato txt)
        Mesh()->Print(filecM);
    }
    if(0)
    {
        int64_t neq = Mesh()->Solution().Rows();
        for (int is = 0; is<neq; is++) {
            Mesh()->Solution()(is,0) = 1.;
        }
        Mesh()->LoadSolution((Mesh()->Solution()));
        TPZMultiphysicsCompMesh *mphys = dynamic_cast<TPZMultiphysicsCompMesh *>(Mesh());
        mphys->LoadSolutionFromMultiPhysics();
    }

    // Transfer object
//    TPZNSBuildTransfers * Transfer = new TPZNSBuildTransfers;
//    Transfer->SetSimulationData(m_simulation_data);
//    //Transfer->Fill_u_To_Mixed(this->Mesh(),0);
//    Transfer->BuildTransferData(this->Mesh());
//    SetTransfer(Transfer);

    //m_simulation_data->SetTransientQ(false);
    this->ExecuteOneTimeStep();
    AssembleLastStep();
    //m_simulation_data->SetTransientQ(true);

    //SetLastState();
//    this->ExecuteOneTimeStep();
//    UpdateMemory_LastStep();

    //std::cout<<this->Rhs()<<std::endl;
    for (int it = 0; it < n_time_steps; it++) {
        for (int k = 1; k <= n_max_fss_iterations; k++) {

            SetCurrentState();
            this->ExecuteOneTimeStep();
            //UpdateMemory_LastStep();


            REAL rhsnomr = this->Get_error();
            REAL dunorm = this->Get_dU_norm();
            error_stop_criterion_Q = (rhsnomr < res_norm);
            dU_stop_criterion_Q = (dunorm < dU_norm);
            //this->PostProcessTimeStep(file_NavierStokes_test);

            if ((error_stop_criterion_Q && (k > n_enforced_fss_iterations)) && dU_stop_criterion_Q) {
                //this->PostProcessTimeStep(file_NavierStokes);
                std::cout << "TPZNSAnalysis:: Iterative process converged with residue norm  = " << this->Get_error() << std::endl;
                UpdateState();
                break;
            }
            bool postprocessingQ = m_simulation_data->IsPostProcessingActivatedQ();
            if (error_stop_criterion_Q && dU_stop_criterion_Q && postprocessingQ) {
                std::cout << "Running post processing" << std::endl;
                this->PostProcessTimeStep(file_NavierStokes);
            }

//            if (postprocessingQ) {
//                std::cout << "Running post processing" << std::endl;
//                this->PostProcessTimeStep(file_NavierStokes);
//            }
        }
        m_simulation_data->UpdateTime();
        //UpdateMemory_LastStep();
        AssembleLastStep();

    }

}

void TPZNSAnalysis::SetInitialParameters(){
    
    // Updating volumetric parameters :
    
    return;
    
}



void TPZNSAnalysis::UpdateParameters(){
    
    // Updating volumetric parameters :

    return;
    
}




void TPZNSAnalysis::AdjustIntegrationOrder(TPZCompMesh * cmesh_o, TPZCompMesh * cmesh_d){
    
    // Assuming the cmesh_o as directive.

    return;
}

void TPZNSAnalysis::ExecuteNewtonIteration(){

    //std::cout<<this->Rhs()<<std::endl;
    this->Assemble();
    if(0)
    {
        std::ofstream plotNavierStiff("NavierStiffness_NewtonIteration.txt");
        fSolver->Matrix()->Print("ek = ",plotNavierStiff,EMathematicaInput);
        std::ofstream plotNavierRhs("NavierRhs_NewtonIteration.txt");;
        fRhs.Print("Rhs =",plotNavierRhs,EMathematicaInput);
    }

//    if(m_simulation_data->IsTransientQ()){
//        m_LastStepRhs.Resize(Rhs().Rows(),Rhs().Cols());
//        this->Rhs()+=m_LastStepRhs;
//    }

    this->Rhs() *= 1.0;
    {
        std::cout<< "residual norm 0.5 = " << Norm(this->Rhs()) <<std::endl;
    }
    this->Solve();
    {
        std::cout<< "solution norm 0.5 = " << Norm(this->Solution()) <<std::endl;
    }
}

void TPZNSAnalysis::LoadCurrentState(){
     fCompMesh->Solution() = m_U_Plus;
    // LoadSolution(m_U_Plus);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());

}

void TPZNSAnalysis::LoadLastState(){
    LoadSolution(m_U_n);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(m_mesh_vec, Mesh());
}

void TPZNSAnalysis::SetCurrentState(){

    TPZCompMesh * mesh = this->Mesh();
    std::map<int, TPZMaterial * >::iterator matit;
    for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
    {
        if(!matit->second) continue;
        TPZMHMNavierStokesMaterial * MHMNavierStokeMaterial = dynamic_cast< TPZMHMNavierStokesMaterial *>(matit->second);
        if (MHMNavierStokeMaterial){
            MHMNavierStokeMaterial->SetCurrentState();
        }
    }

}

void TPZNSAnalysis::SetLastState(){

    TPZCompMesh * mesh = this->Mesh();
    std::map<int, TPZMaterial * >::iterator matit;
    for(matit = mesh->MaterialVec().begin(); matit != mesh->MaterialVec().end(); matit++)
    {
        if(!matit->second) continue;
        TPZMHMNavierStokesMaterial * MHMNavierStokeMaterial = dynamic_cast< TPZMHMNavierStokesMaterial *>(matit->second);
        if (MHMNavierStokeMaterial){
            MHMNavierStokeMaterial->SetLastState();
        }
    }

}

void TPZNSAnalysis::AssembleLastStep(){

    SetLastState();
    AssembleResidual();
    m_LastStepRhs = this->Rhs();

}

void TPZNSAnalysis::AcceptTimeStepSolution(){
    
    bool state = m_simulation_data->IsCurrentStateQ();
    if (state) {
        // must accept solution changes a global data structure shared by the material objects
        // which indicates the solution should be overwritten in memory
        m_simulation_data->Set_must_accept_solution_Q(true);
        // load current state copies m_X_n into the solution vector
        LoadCurrentState();
        // puts the solution vector into a variable depending on yet another global variable
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }else{
        // m_simulation_data is pointer shared by the material object
        // this call forces the solution to be loaded into the memory object
        m_simulation_data->Set_must_accept_solution_Q(true);
        // put m_X in the mesh solution
        LoadLastState();
        // copy the state vector into the memory because must_accept_solution_Q in the m_simulation_data is true
        AssembleResidual();
        m_simulation_data->Set_must_accept_solution_Q(false);
    }
}



/** @brief Update memory using the Transfer object at state n */
void TPZNSAnalysis::UpdateMemory_LastStep(){

    int matid = m_simulation_data->Get_volumetric_material_id()[0];

    TPZMaterial * material = this->Mesh()->FindMaterial(matid);
    TPZMatWithMem<TPZNSMemory,TPZDiscontinuousGalerkin> * associated_material = dynamic_cast<TPZMatWithMem<TPZNSMemory,TPZDiscontinuousGalerkin> *>(material);

    if (!associated_material) {
        DebugStop();
    }

    long N_ipoints = associated_material->GetMemory().get()->NElements();
    for (int ip_index = 0 ; ip_index < N_ipoints; ip_index++) {
        TPZNSMemory & memory = associated_material->GetMemory().get()->operator[](ip_index);
        TPZManVector<REAL,3> u_last(3,0.);
        u_last = memory.u();
        memory.Set_u_last(u_last);
    }

}
