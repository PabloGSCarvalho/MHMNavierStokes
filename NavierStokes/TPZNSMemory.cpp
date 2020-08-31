//
//  TPZNSMemory.cpp
//  PZ
//
//  Created by Pablo Carvalho on 08/31/20.
//
//

#include "TPZNSMemory.h"


/** @brief Default constructor */
TPZNSMemory::TPZNSMemory(){
    
    fu.Resize(3, 0.0);
    fu_last.Resize(3, 0.0);
    
    // Required
    fp              = 0.0;
    fp_last            = 0.0;
    fp_avg          = 0.0;
    fp_avg_last        = 0.0;
    fx.Resize(3,0.0);
    frhs_last = 0.0;
    
}

TPZNSMemory::TPZNSMemory(const TPZNSMemory &copy)
{

    fu = copy.fu;
    fu_last = copy.fu_last;
    fp = copy.fp;
    fp_last = copy.fp_last;
    fp_avg = copy.fp_avg;
    fp_avg_last = copy.fp_avg_last;
    frhs_last = copy.frhs_last;
    fx = copy.fx;
}

TPZNSMemory &TPZNSMemory::operator=(const TPZNSMemory &cp)
{

    fp = cp.fp;
    fp_last = cp.fp_last;
    fp_avg = cp.fp_avg;
    fp_avg_last = cp.fp_avg_last;
    fu = cp.fu;
    fu_last = cp.fu_last;
    frhs_last = cp.frhs_last;
    fx = cp.fx;

    return *this;
}

/** @brief Default destructor */
TPZNSMemory::~TPZNSMemory(){
    
}


const std::string TPZNSMemory::Name() const{
    return "TPZNSMemory";
}

void TPZNSMemory::Write(TPZStream &buf, int withclassid) const{

}

void TPZNSMemory::Read(TPZStream &buf, void *context){

}

void TPZNSMemory::Print(std::ostream &out) const{
    out << "TPZNSMemory " << std::endl;
}


int TPZNSMemory::ClassId() const {
    return Hash("TPZNSMemory");
}