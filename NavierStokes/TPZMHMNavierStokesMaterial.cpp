/*
 *  TPZMHMNavierStokesMaterial.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 10/05/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMHMNavierStokesMaterial.h"
#include "TPZBndCondT.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"

using namespace std;

void TPZMHMNavierStokesMaterial::FillDataRequirementsInterface(TPZMaterialDataT<STATE> &data, std::map<int,TPZMaterialDataT<STATE> > &datavec_left, std::map<int, TPZMaterialDataT<STATE> > &datavec_right)
{
    //TPZMaterial::FillDataRequirementsInterface(data, datavec_left, datavec_right);
    int nref_left = datavec_left.size();
    datavec_left[0].fNeedsNormal = true;
    datavec_right[1].fNeedsNormal = true;
    datavec_left[0].fNeedsDeformedDirectionsFad = NeedsNormalVecFad;
    
}

void TPZMHMNavierStokesMaterial::ContributeInterface(const TPZMaterialDataT<STATE> &data,
                                 std::map<int, TPZMaterialDataT<STATE>> &datavecleft,
                                 std::map<int, TPZMaterialDataT<STATE>> &datavecright, REAL weight,
                                 TPZFMatrix<STATE> &ek,
                                                     TPZFMatrix<STATE> &ef) {


    //2 = 1 Vel space + 1 Press space for datavecleft

    if(fState==ELastState){
        return;
    }

    int nrefleft =  datavecleft.size();
    if (nrefleft != 4 ) {
    //    std::cout << " Erro. The size of the datavec is different from 2 \n";
    //    DebugStop();
    }

    //2 = 1 Vel space + 1 Press space for datavecright
    int nrefright =  datavecright.size();
    if (nrefright != 4 ) {
   //     std::cout << " Erro. The size of the datavec is different from 2 \n";
   //     DebugStop();
    }

    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if(datavecleft.find(vindex) == datavecleft.end()) DebugStop();
    if(datavecleft.find(pindex) == datavecleft.end()) DebugStop();
    if(datavecright.find(vindex) == datavecleft.end()) DebugStop();
    if(datavecright.find(pindex) == datavecleft.end()) DebugStop();

    TPZMaterialDataT<STATE> &vdataleft = datavecleft.find(vindex)->second;
    TPZMaterialDataT<STATE> &pdataleft = datavecleft.find(pindex)->second;
    TPZMaterialDataT<STATE> &vdataright = datavecright.find(vindex)->second;
    TPZMaterialDataT<STATE> &pdataright = datavecright.find(pindex)->second;

    if (vdataleft.fVecShapeIndex.size() == 0) {
        DebugStop();
//        FillVecShapeIndex(vdataleft);
    }

    if (pdataleft.fVecShapeIndex.size() == 0) {
        DebugStop();
//         FillVecShapeIndex(datavecleft[pindex]);
    }

    // Setting the phis
    // V - left
    const TPZFNMatrix<9,REAL>  &tan = pdataleft.axes;

    TPZManVector<STATE,3> u_n = vdataleft.sol[0];

    STATE sLambda_n = pdataright.sol[0][0];

    TPZManVector<STATE,3> Lambda_n = pdataright.sol[0];


    int nshapeV , nshapeP , nshapeLambda, nstateVariablesL;
    nshapeV = vdataleft.fVecShapeIndex.NElements();
    nshapeP = datavecleft[pindex].phi.Rows();
    nshapeLambda = pdataright.phi.Rows();
    int dim = this->Dimension();
    nstateVariablesL = dim-1;
//    if(nstateVariablesL>1){
//        tan = Tangential(datavecright[pindex].normal);
//    }

    int normvecRows = datavecleft[vindex].fDeformedDirections.Rows();
    int normvecCols = datavecleft[vindex].fDeformedDirections.Cols();
    TPZFNMatrix<3,REAL> Normalvec(normvecRows,normvecCols,0.);

    if(datavecleft[vindex].fNeedsDeformedDirectionsFad){
#ifdef _AUTODIFF
        for (int e = 0; e < normvecRows; e++) {
            for (int s = 0; s < normvecCols; s++) {
                Normalvec(e,s)=datavecleft[vindex].fDeformedDirectionsFad(e,s).val();
            }
        }
#else
        DebugStop();
#endif
    }else{
        Normalvec=datavecleft[vindex].fDeformedDirections;
    }
    //qwqw
    for(int i1 = 0; i1 < nshapeV; i1++)
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;

        TPZFNMatrix<9, STATE> phiVi(3,1,0.),lambda_n(3,1,0.);
        for (int e=0; e< 3 ; e++) {
            phiVi(e,0)=Normalvec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
        }

        STATE Lambda_dot_phiV = 0.; //rhs term
        for(int f =0; f< nstateVariablesL ; f++){
            for (int e=0; e< 3 ; e++) {
                lambda_n(e,0) += sLambda_n*tan(f,e);
            }
        }

        Lambda_dot_phiV = InnerVec(lambda_n, phiVi);
        ef(i1) += -fMultiplier * weight*Lambda_dot_phiV;
//        TPZManVector<REAL,3>  normal = datavecleft[vindex].normal;
//        TPZManVector<REAL,3>  cood = data.x;
//        if(normal[1]=-1&&cood[1]==0.){
//            for(int j1 = 0; j1 < nshapeV; j1++) {
//                int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
//                int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
//
//                TPZFNMatrix<9, STATE> phiVj(3,1);
//                for (int e=0; e< 3 ; e++) {
//                    phiVj(e,0)=Normalvec(e,jvec1)*datavecleft[vindex].phi(jphi1,0);
//                }
//                STATE fact = fMultiplier * weight * InnerVec(phiVi,phiVj);
//                ek(i1,j1) += fact;
//            }
//        }

        // K12 e K21 - (test V left) * (trial Lambda right)
        for(int j1 = 0; j1 < nshapeLambda; j1++)
        {
            // Var. Sigma, Sn :
            TPZFNMatrix<9, STATE> lambda_j(3,1,0.);
            TPZFNMatrix<9, STATE> phiLamb = datavecright[pindex].phi;

            for(int f =0; f< nstateVariablesL ; f++){
                lambda_j.Zero();
                for (int e=0; e< 3 ; e++) {
                    lambda_j(e, 0) += phiLamb(j1, 0) * tan(f, e);
                }

                STATE fact = fMultiplier * weight * InnerVec(phiVi,lambda_j);
                ek(i1,j1*nstateVariablesL+f+nshapeV) += fact;
                ek(j1*nstateVariablesL+f+nshapeV,i1) += fact;

            }

        }
    }

    for(int i1 = 0; i1 < nshapeLambda; i1++)
    {
        // Var. Sigma, Sn :
        TPZFNMatrix<9, STATE> lambda_i(3,1,0.);
        TPZFNMatrix<9, STATE> phiLamb = datavecright[pindex].phi;

        // Tangencial comp. vector (t x t)Sn :
        for(int f =0; f< nstateVariablesL ; f++) {
            lambda_i.Zero();
            for (int e=0; e< 3 ; e++) {
                lambda_i(e, 0) += phiLamb(i1, 0) * tan(f, e);
            }
        }

        STATE phiLambda_dot_U = 0.; //rhs term
        for (int e=0; e< 3 ; e++) {
                phiLambda_dot_U += lambda_i(e, 0) * u_n[e];
        }
        ef(i1+nshapeV) += -fMultiplier * weight*phiLambda_dot_U;
    }

    //lalala
    //ef.Zero();

}


void TPZMHMNavierStokesMaterial::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                          REAL weight, TPZFMatrix<STATE> &ek,
                          TPZFMatrix<STATE> &ef,
                                              TPZBndCondT<STATE> &bc) {

    if(fState==ELastState){
        return;
    }
    int nshapeV , nshapeLambda;
    
    nshapeV = datavec[0].phi.Rows();
    nshapeLambda = datavec[1].phi.Rows();
    int nstateBC = 1;

    TPZManVector<REAL,3> &normal = datavec[0].normal;
    if(nshapeLambda){
        nstateBC = ek.Rows()/nshapeLambda;
        normal = datavec[1].normal;
    }

    TPZFNMatrix<9,REAL>  &tan = datavec[1].axes;

    TPZManVector<REAL,3> u_n  = datavec[0].sol[0];
    TPZManVector<REAL,3> l_n  = datavec[1].sol[0];


    if(nshapeV!=0&&nshapeLambda!=0){
        DebugStop();
    }
    
    //TPZFNMatrix<9, STATE> u_D(dim,1);
    TPZVec<STATE> v_Dirichlet(3,0.);
    TPZVec<STATE> v_Neumann(3,0.);
    STATE p_D =0.;
    
    TPZFMatrix<STATE> * phi;
    
    TPZManVector<REAL> x;
    
    if (nshapeV!=0) {
        x=datavec[0].x;
    }else{
        x=datavec[1].x;
    }
    
    TPZFMatrix<STATE> gradu(3,3,0.);
    TPZManVector<STATE> vbc(4,0.);
    TPZFMatrix<STATE> Du(3,3,0.),Dun(3,1,0.),u_x_beta(3,3,0.), u_x_beta_n(3,1,0.);
    if(bc.HasForcingFunctionBC())
    {
        
        bc.ForcingFunctionBC()(x,vbc,gradu);
        //std::cout<<gradu<<std::endl;
        v_Dirichlet[0] = vbc[0];
        v_Dirichlet[1] = vbc[1];
        v_Dirichlet[2] = vbc[2];
        p_D = vbc[3];

        REAL time = f_sim_data->GetTime();
//        v_Dirichlet[0] = v_Dirichlet[0]*(sin(M_PI*(time+0.1-0.5))*0.5+0.5);
//        v_Dirichlet[1] = v_Dirichlet[1]*(sin(M_PI*(time+0.1-0.5))*0.5+0.5);
//        if(time<0.9){
//            v_Dirichlet[0] = v_Dirichlet[0]*(time+0.1);
//            v_Dirichlet[1] = v_Dirichlet[1]*(time+0.1);
//        }

        // Calculo do vetor (Du)n :
        if (nshapeV!=0) {
            gradu.Resize(3,3);
            //std::cout<<gradu<<std::endl;
            STATE visco = GetViscosity();
            for (int i = 0; i<3; i++) {
                for (int j = 0; j<3; j++) {
                    Du(i,j)=  visco * (gradu(i,j) + gradu(j, i));
                }
            }

            //std::cout<<Du<<std::endl;
            for (int i = 0; i<3; i++) {
                for (int j = 0; j<3; j++) {
                    Dun(i,0) +=  Du(i,j) * datavec[0].normal[j];
                }
            }
            //std::cout<<datavec[0].normal<<std::endl;
            //std::cout<<Dun<<std::endl;

            // Neumann vector

            for (int i = 0; i<3; i++) {
                v_Neumann[i] = Dun(i,0) - p_D * datavec[0].normal[i];
            }
            
        }
        //qwqw
        
    }else{


        for (int i = 0; i<3; i++) {
            v_Dirichlet[i] = bc.Val2()[i];
        }

        for (int i = 0; i<3; i++) {
            for (int j = 0; j<3; j++) {
                v_Neumann[i] +=  bc.Val1()(i,j) * normal[j];
            }
        }

    }
    
    // Calculo do vetor (Du)n :
    if (nshapeLambda!=0) {
        gradu.Resize(3,3);
        STATE visco = GetViscosity();
        for (int i = 0; i<3; i++) {
            for (int j = 0; j<3; j++) {
                Du(i,j)=  visco * (gradu(i,j) + gradu(j, i));
            }
        }
        Dun.Zero();
        for (int i = 0; i<3; i++) {
            for (int j = 0; j<3; j++) {
                Dun(i,0) +=  Du(i,j) * datavec[1].normal[j];
            }
        }

        // Neumann vector
        
        for (int i = 0; i<3; i++) {
            v_Neumann[i] = Dun(i,0) - p_D * datavec[1].normal[i] ;
        }
  
    }

    //STATE value = 0.;
    TPZVec<STATE> v_value(nstateBC,0.);

    if(nshapeV!=0){
        
        if (bc.Type()==0) {

            for (int i = 0; i<3; i++) {
                v_value[0] += v_Dirichlet[i]*datavec[0].normal[i];
            }
            
//            std::cout << "v_d = "<< v_Dirichlet << std::endl;
//            std::cout << "u_n = "<< u_n << std::endl;
//            std::cout << "normal = "<< datavec[0].normal << std::endl;
//            std::cout << "val 1 = "<< value << std::endl;
//
            v_value[0] = v_value[0] - u_n[0];

//            std::cout << "val 2 = "<< value << std::endl;
//            std::cout << "________" << std::endl;
            

        }else if (bc.Type()==1){
            
            for (int j = 0; j<3; j++) {
                v_value[0] += v_Neumann[j] * datavec[0].normal[j];
            }
            
            
        }
        
        phi = &datavec[0].phi;
        
    }else if (nshapeLambda!=0&&(bc.Type()!=5)){

        if(bc.Type()==0){
            for(int f =0; f< nstateBC ; f++) {
                for (int i = 0; i<3; i++) {
                    v_value[f] += v_Neumann[i] * tan(f, i);
                }
            }
                //                        std::cout << "v_n = "<< v_Neumann << std::endl;
                //                        std::cout << "l_n = "<< l_n << std::endl;
                //                        std::cout << "normal = "<< datavec[1].axes << std::endl;
                //                        std::cout << "val 1 = "<< value << std::endl;
                //
            //value = value - l_n[0];
            for(int f =0; f< nstateBC ; f++) {
                v_value[f] = v_value[f] - l_n[f];
            }


                //                        std::cout << "val 2 = "<< value << std::endl;
                //                        std::cout << "________" << std::endl;

        }else if(bc.Type()==1){

            for(int f =0; f< nstateBC ; f++) {
                for (int i = 0; i<3; i++) {
                    v_value[f] += v_Dirichlet[i] * tan(f, i);
                }
            }
        }


        phi = &datavec[1].phi;
        
    }else if (bc.Type()==5){
        
        phi = &datavec[1].phi;
        
    }else {
        DebugStop();
    }
    
    switch (bc.Type()) {
        case 0: //Dirichlet for continuous formulation
        {
            for(int istate = 0; istate < nstateBC; istate++) {
                for (int j1 = 0; j1 < phi->Rows(); j1++) {
                    ef(j1*nstateBC+istate, 0) += fBigNumber * v_value[istate] * (*phi)(j1, 0) * weight;

                    for (int i1 = 0; i1 < phi->Rows(); i1++) {
                        ek(i1*nstateBC+istate, j1*nstateBC+istate) += fBigNumber * weight * (*phi)(j1, 0) * (*phi)(i1, 0);
                    }
                }
            }

        }
            break;
            
            
        case 1: //Neumann for continuous formulation
        {
            //Oscillations
            REAL time = f_sim_data->GetTime();
            REAL valbc1 = bc.Val1()(2,2);
            if(nshapeLambda!=0&&valbc1!=0){
                if((time>=0.1)&&(time<=1.4)){
                    v_value[0] = 0.5;
                }
                if((time>=1.5)&&(time<=2.9)){
                    v_value[0] = -0.6;
                }
            }

            for(int j1 = 0; j1 < phi->Rows(); j1++)
            {
                for(int istate = 0; istate < nstateBC; istate++) {
                    ef(j1*nstateBC+istate, 0) += v_value[istate] * (*phi)(j1, 0) * weight;
                }
            }
        }
            break;
            
        case 5: //Ponto pressao
        {

            //return;
            p_D = bc.Val2()[0];
            
            
            for(int i = 0; i < phi->Rows(); i++ )
            {
                ef(i) += 1.0 * p_D * (*phi)(i,0);
                
                for(int j = 0; j < phi->Rows(); j++){
    
                    ek(i,j) += 1.0 * ((*phi)(i,0) * (*phi)(j,0));
    
                }
                
            }
            
        }
            break;

        case 7:
        {
            for(int j1 = 0; j1 < phi->Rows(); j1++)
            {
                for(int istate = 0; istate < nstateBC; istate++) {
                    for (int j1 = 0; j1 < phi->Rows(); j1++) {
                       // ef(j1*nstateBC+istate, 0) += -0. * (*phi)(j1, 0) * weight;

                        for (int i1 = 0; i1 < phi->Rows(); i1++) {
                            ek(i1*nstateBC+istate, j1*nstateBC+istate) +=  -weight * (*phi)(j1, 0) * (*phi)(i1, 0);
                        }
                    }
                }
            }
        }
            break;

        default:
        {
            std::cout << "Boundary not implemented " << std::endl;
            DebugStop();
        }
            break;
    }

}

void TPZMHMNavierStokesMaterial::ContributeBCInterface(const TPZMaterialDataT<STATE> &data,
                      std::map<int, TPZMaterialDataT<STATE>> &datavec,
                      REAL weight,
                      TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
                                                       TPZBndCondT<STATE> &bc) {
    
    
    
    STATE rhsnorm = Norm(ef);
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavec[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[vindex]);
    }
    
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    //Normal
    const TPZManVector<REAL,3> &normal = data.normal;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    // nshapeV = phiV.Rows()*NStateVariables();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    //Dirichlet
    
    auto v_2=bc.Val2();
    TPZFMatrix<STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    
    switch (bc.Type()) {
            
            
        case 0: //Dirichlet for continuous formulation
        {
            
            if(bc.HasForcingFunctionBC())
            {
                TPZManVector<STATE> vbc(3);
                TPZFMatrix<STATE> gradu;
                bc.ForcingFunctionBC()(datavec[vindex].x,vbc,gradu);
                v_2[0] = vbc[0];
                v_2[1] = vbc[1];
                p_D=vbc[2];
                
            }
            
            //Componente tangencial -> imposta fracamente:
            
            for(int i = 0; i < nshapeV; i++ )
            {
                int iphi = datavec[vindex].fVecShapeIndex[i].second;
                int ivec = datavec[vindex].fVecShapeIndex[i].first;
                TPZFNMatrix<9,STATE> GradVni(fDimension,1,0.),phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                GradVni.Zero();
                
                TPZFNMatrix<4,STATE> GradVi(fDimension,fDimension,0.),GradVit(fDimension,fDimension,0.),Dui(fDimension,fDimension,0.),Duni(fDimension,1,0.);
                
                for (int e=0; e<fDimension; e++) {
                    
                    for (int f=0; f<fDimension; f++) {
                        GradVi(e,f) = datavec[vindex].fDeformedDirections(e,ivec)*dphiVx(f,iphi);
                        //termo transposto:
                        GradVit(f,e) = datavec[vindex].fDeformedDirections(e,ivec)*dphiVx(f,iphi);
                        
                    }
                }
                
                //Du = 0.5(GradU+GradU^T)
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        Dui(e,f)= (1./2.) * (GradVi(e,f) + GradVit(e,f));
                    }
                }
                
                //Duni
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        Duni(e,0) += Dui(e,f)*normal[f] ;
                    }
                }
                
                //GradVni
                for (int e=0; e<fDimension; e++) {
                    for (int f=0; f<fDimension; f++) {
                        GradVni(e,0) += GradVi(e,f)*normal[f] ;
                    }
                }
                
                
                for (int e=0; e<fDimension; e++) {
                    phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*datavec[vindex].phi(iphi,0);
                    phiVni(0,0)+=phiVi(e,0)*normal[e];
                    
                }
                
                TPZManVector<REAL> n = data.normal;
                TPZManVector<REAL> t(2);
                t[0]=n[1];
                t[1]=n[0];
                
                
                
                phiVti(0,0)= t[0] * phiVi(0,0) + t[1] * phiVi(1,0);
                TPZFNMatrix<9,STATE> phiVtit(fDimension,1,0.);
                phiVtit(0,0)=phiVti(0,0)*t[0];
                phiVtit(1,0)=phiVti(0,0)*t[1];
                
                TPZFNMatrix<9,STATE> phiVnin(fDimension,1,0.);
                phiVnin(0,0)=phiVni(0,0)*n[0];
                phiVnin(1,0)=phiVni(0,0)*n[1];
                
                
                if(fSpace==1||fSpace==3){
                    
                    REAL vh_t = 0.;
                    if(v_h.size()>1){
                        vh_t = v_h[1];
                    }
                    
                    REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
                    
                    TPZManVector<REAL> v_tt(2);
                    v_tt[0]=v_t*t[0];
                    v_tt[1]=v_t*t[1];
                    
                    TPZManVector<REAL> vh_tt(2);
                    vh_tt[0]=vh_t*t[0];
                    vh_tt[1]=vh_t*t[1];
                    
                    TPZFNMatrix<9,STATE> diffvt(fDimension,1,0.);
                    diffvt(0,0)=v_tt[0];
                    diffvt(1,0)=v_tt[1];
                    
                    
                    STATE factef = weight * fSigma * v_t * phiVti(0,0) * fViscosity;
                    
                    ef(i,0) += factef;
                    
                    
                    STATE fact= 2. * weight * fViscosity * InnerVec(diffvt, Duni) ;
                    
                    ef(i,0) += fTheta*fact;
                    
                    
                    for(int j = 0; j < nshapeV; j++){
                        int jphi = datavec[vindex].fVecShapeIndex[j].second;
                        int jvec = datavec[vindex].fVecShapeIndex[j].first;
                        
                        TPZFNMatrix<9,STATE> GradVnj(fDimension,1),phiVtj(1,1,0.),phiVj(fDimension,1);
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVj(e,0)=datavec[vindex].fDeformedDirections(e,jvec)*datavec[vindex].phi(jphi,0);
                        }
                        
                        
                        phiVtj(0,0)= t[0] * phiVj(0,0) + t[1] * phiVj(1,0);
                        
                        
                        
                        TPZFNMatrix<4,STATE> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);
                        
                        for (int e=0; e<fDimension; e++) {
                            
                            for (int f=0; f<fDimension; f++) {
                                GradVj(e,f) = datavec[vindex].fDeformedDirections(e,jvec)*dphiVx(f,jphi);
                                //termo transposto:
                                GradVjt(f,e) = datavec[vindex].fDeformedDirections(e,jvec)*dphiVx(f,jphi);
                                
                            }
                        }
                        
                        //Du = 0.5(GradU+GradU^T)
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Duj(e,f)= (1./2.) * (GradVj(e,f) + GradVjt(e,f));
                            }
                        }
                        
                        //Du2nj
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Dunj(e,0) += Duj(e,f)*normal[f] ;
                            }
                        }
                        
                        STATE factek = weight * fSigma * phiVtj(0,0)* phiVti(0,0) * fViscosity;
                        ek(i,j) +=  factek;
                        
                        
                        STATE fact =(-1.) * weight * 2. * fViscosity * InnerVec(phiVtit, Dunj) ;
                        ek(i,j) += fact ;
                        ek(j,i) += -fTheta*fact;
                        
                        
                    }
                    
                    
                    //Componente normal -> imposta fortemente:
                    if(fSpace==1||fSpace==3){
                        for(int i = 0; i < nshapeV; i++ )
                        {
                            
                            int iphi = datavec[vindex].fVecShapeIndex[i].second;
                            int ivec = datavec[vindex].fVecShapeIndex[i].first;
                            TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                            
                            
                            for (int e=0; e<fDimension; e++) {
                                phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*datavec[vindex].phi(iphi,0);
                                phiVni(0,0)+=phiVi(e,0)*n[e];
                                phiVti(0,0)+=phiVi(e,0)*t[e];
                            }
                            
                            
                            REAL vh_n = v_h[0];
                            REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                            
                            ef(i,0) += -weight * fBigNumber * (vh_n-v_n) * (phiVni(0,0));
                            
                            
                            for(int j = 0; j < nshapeV; j++){
                                
                                int jphi = datavec[vindex].fVecShapeIndex[j].second;
                                int jvec = datavec[vindex].fVecShapeIndex[j].first;
                                
                                TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
                                
                                for (int e=0; e<fDimension; e++) {
                                    phiVj(e,0)=datavec[vindex].fDeformedDirections(e,jvec)*datavec[vindex].phi(jphi,0);
                                    phiVnj(0,0)+=phiVj(e,0)*n[e];
                                    phiVtj(0,0)+=phiVj(e,0)*t[e];
                                    
                                }
                                
                                ek(i,j) += weight * fBigNumber * (phiVni(0,0)) * (phiVnj(0,0)) ;
                                
                            }
                            
                        }
                        
                        
                    }
                    
                }
                
                
            }
            break;
            
        case 1: //Neumann for continuous formulation
            {
                
                
                if(bc.HasForcingFunctionBC())
                {
                    TPZManVector<STATE> vbc(3);
                    TPZFMatrix<STATE> gradu;
                    bc.ForcingFunctionBC()(datavec[vindex].x,vbc,gradu);
                    v_2[0] = vbc[0];
                    v_2[1] = vbc[1];
                    p_D=vbc[2];
                    
                }
                
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
                    TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*phiV(iphi,0);
                    }
                    
                    TPZManVector<REAL> n = data.normal;
                    
                    TPZFNMatrix<9,STATE> pn(fDimension,1);
                    
                    
                    for (int f=0; f<fDimension; f++) {
                        pn(f,0)=n[f]*v_1(0,0);
                    }
                    
                    //Adaptação para Hdiv
                    
                    STATE factef=0.0;
                    
                    factef += InnerVec(pn, phiVi) ;
                    
                    
                    ef(i,0) += weight * factef;
                    
                }
                
                
            }
            
            
            
            break;
            
        case 2: //Penetração com slip for continuous formulation
            {
                
                if(bc.HasForcingFunctionBC())
                {
                    TPZManVector<STATE> vbc(3);
                    TPZFMatrix<STATE> gradu;
                    bc.ForcingFunctionBC()(datavec[vindex].x,vbc,gradu);
                    v_2[0] = vbc[0];
                    v_2[1] = vbc[1];
                    p_D=vbc[2];
                    
                }
                
                
                if(fSpace==1||fSpace==2){
                    TPZManVector<REAL> n = data.normal;
                    TPZManVector<REAL> t(2);
                    t[0]=-n[1];
                    t[1]=n[0];
                    
                    //Componente normal -> imposta fortemente:
                    
                    for(int i = 0; i < nshapeV; i++ )
                    {
                        
                        int iphi = datavec[vindex].fVecShapeIndex[i].second;
                        int ivec = datavec[vindex].fVecShapeIndex[i].first;
                        TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                        
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*datavec[vindex].phi(iphi,0);
                            phiVni(0,0)+=phiVi(e,0)*n[e];
                            phiVti(0,0)+=phiVi(e,0)*t[e];
                        }
                        
                        REAL vh_n = v_h[0];
                        REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                        
                        ef(i,0) += -weight * fBigNumber * (vh_n-v_n) * (phiVni(0,0));
                        
                        
                        for(int j = 0; j < nshapeV; j++){
                            
                            int jphi = datavec[vindex].fVecShapeIndex[j].second;
                            int jvec = datavec[vindex].fVecShapeIndex[j].first;
                            
                            TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
                            
                            for (int e=0; e<fDimension; e++) {
                                phiVj(e,0)=datavec[vindex].fDeformedDirections(e,jvec)*datavec[vindex].phi(jphi,0);
                                phiVnj(0,0)+=phiVj(e,0)*n[e];
                                phiVtj(0,0)+=phiVj(e,0)*t[e];
                                
                            }
                            
                            ek(i,j) += weight * fBigNumber * (phiVni(0,0)) * (phiVnj(0,0)) ;
                            
                        }
                        
                    }
                }
                
                
                if(fSpace==3){
                    
                    for(int i = 0; i < nshapeV; i++ )
                    {
                        int iphi = datavec[vindex].fVecShapeIndex[i].second;
                        int ivec = datavec[vindex].fVecShapeIndex[i].first;
                        TPZFNMatrix<9,STATE> GradVni(fDimension,1,0.),phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                        GradVni.Zero();
                        
                        TPZFNMatrix<4,STATE> GradVi(fDimension,fDimension,0.),GradVit(fDimension,fDimension,0.),Dui(fDimension,fDimension,0.),Duni(fDimension,1,0.);
                        
                        for (int e=0; e<fDimension; e++) {
                            
                            for (int f=0; f<fDimension; f++) {
                                GradVi(e,f) = datavec[vindex].fDeformedDirections(e,ivec)*dphiVx(f,iphi);
                                //termo transposto:
                                GradVit(f,e) = datavec[vindex].fDeformedDirections(e,ivec)*dphiVx(f,iphi);
                                
                            }
                        }
                        
                        //Du = 0.5(GradU+GradU^T)
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Dui(e,f)= (1./2.) * (GradVi(e,f) + GradVit(e,f));
                            }
                        }
                        
                        //Duni
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Duni(e,0) += Dui(e,f)*normal[f] ;
                            }
                        }
                        
                        //GradVni
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                GradVni(e,0) += GradVi(e,f)*normal[f] ;
                            }
                        }
                        
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*datavec[vindex].phi(iphi,0);
                            phiVni(0,0)+=phiVi(e,0)*normal[e];
                        }
                        
                        TPZManVector<REAL> n = data.normal;
                        TPZManVector<REAL> t(2);
                        t[0]=n[1];
                        t[1]=n[0];
                        
                        phiVti(0,0)= t[0] * phiVi(0,0) + t[1] * phiVi(1,0);
                        TPZFNMatrix<9,STATE> phiVtit(fDimension,1,0.);
                        phiVtit(0,0)=phiVti(0,0)*t[0];
                        phiVtit(1,0)=phiVti(0,0)*t[1];
                        
                        TPZFNMatrix<9,STATE> phiVnin(fDimension,1,0.);
                        phiVnin(0,0)=phiVni(0,0)*n[0];
                        phiVnin(1,0)=phiVni(0,0)*n[1];
                        
                        
                        //Componente normal -> imposta fortemente:
                        
                        for(int i = 0; i < nshapeV; i++ )
                        {
                            
                            int iphi = datavec[vindex].fVecShapeIndex[i].second;
                            int ivec = datavec[vindex].fVecShapeIndex[i].first;
                            TPZFNMatrix<9,STATE> phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                            
                            
                            for (int e=0; e<fDimension; e++) {
                                phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*datavec[vindex].phi(iphi,0);
                                phiVni(0,0)+=phiVi(e,0)*n[e];
                                phiVti(0,0)+=phiVi(e,0)*t[e];
                            }
                            
                            
                            REAL vh_n = v_h[0];
                            REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                            
                            ef(i,0) += -weight * fBigNumber * (vh_n-v_n) * (phiVni(0,0));
                            
                            
                            for(int j = 0; j < nshapeV; j++){
                                
                                int jphi = datavec[vindex].fVecShapeIndex[j].second;
                                int jvec = datavec[vindex].fVecShapeIndex[j].first;
                                
                                TPZFNMatrix<9,STATE> phiVj(fDimension,1),phiVnj(1,1,0.),phiVtj(1,1,0.);
                                
                                for (int e=0; e<fDimension; e++) {
                                    phiVj(e,0)=datavec[vindex].fDeformedDirections(e,jvec)*datavec[vindex].phi(jphi,0);
                                    phiVnj(0,0)+=phiVj(e,0)*n[e];
                                    phiVtj(0,0)+=phiVj(e,0)*t[e];
                                    
                                }
                                
                                ek(i,j) += weight * fBigNumber * (phiVni(0,0)) * (phiVnj(0,0)) ;
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
            }
            break;
            
            
            
        case 10: //Penetração com slip for continuous formulation
            {
                
                if(bc.HasForcingFunctionBC())
                {
                    TPZManVector<STATE> vbc(3);
                    TPZFMatrix<STATE> gradu;
                    bc.ForcingFunctionBC()(datavec[vindex].x,vbc,gradu);
                    v_2[0] = vbc[0];
                    v_2[1] = vbc[1];
                    p_D=vbc[2];
                    
                }
                
                //Componente tangencial -> imposta fracamente:
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
                    TPZFNMatrix<9,STATE> GradVni(fDimension,1,0.),phiVi(fDimension,1),phiVni(1,1,0.),phiVti(1,1,0.);
                    GradVni.Zero();
                    
                    TPZFNMatrix<4,STATE> GradVi(fDimension,fDimension,0.),GradVit(fDimension,fDimension,0.),Dui(fDimension,fDimension,0.),Duni(fDimension,1,0.);
                    
                    for (int e=0; e<fDimension; e++) {
                        
                        for (int f=0; f<fDimension; f++) {
                            GradVi(e,f) = datavec[vindex].fDeformedDirections(e,ivec)*dphiVx(f,iphi);
                            //termo transposto:
                            GradVit(f,e) = datavec[vindex].fDeformedDirections(e,ivec)*dphiVx(f,iphi);
                            
                        }
                    }
                    
                    //Du = 0.5(GradU+GradU^T)
                    for (int e=0; e<fDimension; e++) {
                        for (int f=0; f<fDimension; f++) {
                            Dui(e,f)= (1./2.) * (GradVi(e,f) + GradVit(e,f));
                        }
                    }
                    
                    //Duni
                    for (int e=0; e<fDimension; e++) {
                        for (int f=0; f<fDimension; f++) {
                            Duni(e,0) += Dui(e,f)*normal[f] ;
                        }
                    }
                    
                    //GradVni
                    for (int e=0; e<fDimension; e++) {
                        for (int f=0; f<fDimension; f++) {
                            GradVni(e,0) += GradVi(e,f)*normal[f] ;
                        }
                    }
                    
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*datavec[vindex].phi(iphi,0);
                        phiVni(0,0)+=phiVi(e,0)*normal[e];
                        
                    }
                    
                    TPZManVector<REAL> n = data.normal;
                    TPZManVector<REAL> t(2);
                    t[0]=-n[1];
                    t[1]=n[0];
                    
                    
                    
                    phiVti(0,0)= t[0] * phiVi(0,0) + t[1] * phiVi(1,0);
                    TPZFNMatrix<9,STATE> phiVtit(fDimension,1,0.);
                    phiVtit(0,0)=phiVti(0,0)*t[0];
                    phiVtit(1,0)=phiVti(0,0)*t[1];
                    
                    TPZFNMatrix<9,STATE> phiVnin(fDimension,1,0.);
                    phiVnin(0,0)=phiVni(0,0)*n[0];
                    phiVnin(1,0)=phiVni(0,0)*n[1];
                    
                    
                    REAL vh_t = v_h[1];
                    REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
                    TPZManVector<REAL> v_tt(2);
                    v_tt[0]=v_t*t[0];
                    v_tt[1]=v_t*t[1];
                    TPZManVector<REAL> vh_tt(2);
                    vh_tt[0]=vh_t*t[0];
                    vh_tt[1]=vh_t*t[1];
                    TPZFNMatrix<9,STATE> diffvt(fDimension,1,0.);
                    diffvt(0,0)=v_tt[0];
                    diffvt(1,0)=v_tt[1];
                    
                    REAL vh_n = v_h[0];
                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                    TPZManVector<REAL> v_nn(2);
                    v_nn[0]=v_n*n[0];
                    v_nn[1]=v_n*n[1];
                    TPZManVector<REAL> vh_nn(2);
                    vh_nn[0]=vh_n*n[0];
                    vh_nn[1]=vh_n*n[1];
                    TPZFNMatrix<9,STATE> diffvn(fDimension,1,0.);
                    diffvn(0,0)=v_nn[0];
                    diffvn(1,0)=v_nn[1];
                    
                    STATE factefn = weight * fSigma * v_n * phiVni(0,0);
                    ef(i,0) += factefn;
                    STATE factn= 2. * weight * fViscosity * InnerVec(diffvn, Duni);
                    ef(i,0) += fTheta*factn;
                    
                    STATE facteft = weight * fSigma * v_t * phiVti(0,0);
                    ef(i,0) += facteft;
                    STATE factt= 2. * weight * fViscosity * InnerVec(diffvt, Duni);
                    ef(i,0) += fTheta*factt;
                    
                    for(int j = 0; j < nshapeV; j++){
                        int jphi = datavec[vindex].fVecShapeIndex[j].second;
                        int jvec = datavec[vindex].fVecShapeIndex[j].first;
                        TPZFNMatrix<9,STATE> GradVnj(fDimension,1),phiVtj(1,1,0.),phiVnj(1,1,0.),phiVj(fDimension,1);
                        for (int e=0; e<fDimension; e++) {
                            phiVj(e,0)=datavec[vindex].fDeformedDirections(e,jvec)*datavec[vindex].phi(jphi,0);
                        }
                        
                        phiVnj(0,0)= n[0] * phiVj(0,0) + n[1] * phiVj(1,0);
                        phiVtj(0,0)= t[0] * phiVj(0,0) + t[1] * phiVj(1,0);
                        
                        TPZFNMatrix<4,STATE> GradVj(fDimension,fDimension,0.),GradVjt(fDimension,fDimension,0.),Duj(fDimension,fDimension,0.),Dunj(fDimension,1,0.);
                        for (int e=0; e<fDimension; e++) {
                            
                            for (int f=0; f<fDimension; f++) {
                                GradVj(e,f) = datavec[vindex].fDeformedDirections(e,jvec)*dphiVx(f,jphi);
                                //termo transposto:
                                GradVjt(f,e) = datavec[vindex].fDeformedDirections(e,jvec)*dphiVx(f,jphi);
                            }
                        }
                        
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Duj(e,f)= (1./2.) * (GradVj(e,f) + GradVjt(e,f));
                            }
                        }
                        
                        //Du2nj
                        for (int e=0; e<fDimension; e++) {
                            for (int f=0; f<fDimension; f++) {
                                Dunj(e,0) += Duj(e,f)*normal[f] ;
                            }
                        }
                        
                        STATE factek = weight * fSigma * phiVtj(0,0)* phiVti(0,0);
                        ek(i,j) +=  factek;
                        
                        STATE fact =(-1.) * weight * 2. * fViscosity * InnerVec(phiVtit, Dunj) ;
                        ek(i,j) += fact ;
                        ek(j,i) += -fTheta*fact;
                        
                        
                        STATE factekn = weight * fSigma * phiVnj(0,0)* phiVni(0,0);
                        ek(i,j) +=  factekn;
                        
                        STATE factn =(-1.) * weight * 2. * fViscosity * InnerVec(phiVnin, Dunj) ;
                        ek(i,j) += factn ;
                        ek(j,i) += -fTheta*factn;
                        
                    }
                    
                }
                
                
            }
            break;
            
            
            
        }
            
            
            
    }
    
    
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    if(0)
    {
        std::ofstream fileEK("FileEKContributeBCInterf.txt");
        std::ofstream fileEF("FileEFContributeBCInterf.txt");
        ek.Print("MatrizBCint = ",fileEK,EMathematicaInput);
        ef.Print("ForceBCint = ",fileEF,EMathematicaInput);
    }
    
    
    
}


TPZManVector<REAL,3> TPZMHMNavierStokesMaterial::ComputeNormal(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright){

    int vindex = VIndex();
    int pindex = PIndex();
    
    TPZManVector<REAL,3> xcenterL = datavecleft[vindex].XCenter;
    TPZManVector<REAL,3> xcenterR = datavecright[pindex].XCenter;
    TPZManVector<REAL,3> normalV = data.normal;
    
    if(xcenterR[0]>xcenterL[0]){
        normalV[0]=1.;

    }else if(xcenterR[0]>xcenterL[0]){
        normalV[0]=-1.;
    }
    return normalV;
    
    if(xcenterR[1]>xcenterL[1]){
        
    }
    
    return normalV;
    
}

TPZFMatrix<STATE> TPZMHMNavierStokesMaterial::Transpose(TPZFMatrix<STATE> &MatrixU ){

    int dim = Dimension();
    TPZFMatrix<STATE> MatrixUt(dim,dim);
    
    if((MatrixU.Rows()!=dim)&&(MatrixU.Cols()!=dim)){
        DebugStop();
    }
    
    for (int i = 0; i < dim; i++ ) {
        for (int j = 0; j < dim; j++ ) {
            MatrixUt(i,j) = MatrixU(j,i);
         }
    }
    
    return MatrixUt;
    
}

TPZFNMatrix<9,REAL> TPZMHMNavierStokesMaterial::Tangential(TPZManVector<REAL,3> normal){

    TPZFNMatrix<9,REAL>  MatrixTan(3,3);

    for (int i = 0; i < 3; i++ ) {
        for (int j = 0; j < 3; j++ ) {
            MatrixTan(i,j) = normal[i]*normal[j];
            if(i==j){
                MatrixTan(i,j)=1-MatrixTan(i,j);
            }
        }
    }

    return MatrixTan;
}

////////////////////////////////////////////////////////////////////



