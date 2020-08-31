//
//  TPZNSMemory.h
//  PZ
//
//  Created by Pablo Carvalho on 08/31/20.
//
//

#ifndef __PZ__TPZNSMemory__
#define __PZ__TPZNSMemory__

#include <stdio.h>
#include "pzreal.h"
//#include "pzfilebuffer.h"
#include "pzfmatrix.h"


/*! @brief Store the information required on a integration point.
 *         Brief description continued.
 *
 *  Store the solution at n time step.
 *  Also it can store the nonlinear part of the flux at n step.
 *  Store the xyz of the spatial properties.
 */

class TPZNSMemory {

    /**
     * @defgroup Required Memory quantities
     * @{
     */
    
    /** @brief Total flux */
    TPZManVector<REAL,3> fu;
    
    /** @brief Total flux at the previous timestep */
    TPZManVector<REAL,3> fu_last;
    
    /** @brief Weighted pressure */
    REAL fp;
    
    /** @brief Weighted pressure at the previous timestep */
    REAL fp_last;
    
    /** @brief Average weighted pressure */
    REAL fp_avg;
    
    /** @brief Average weighted pressure at the previous timestep */
    REAL fp_avg_last;
    
    /** @brief Right hand side */
    REAL frhs_last;
    
    /** @brief Spatial coordinate */
    TPZManVector<REAL,3> fx;

    
    //@}
    
public:
    
    /** @brief Default constructor */
    TPZNSMemory();
    
    /** @brief Default destructor */
    ~TPZNSMemory();
    
    TPZNSMemory(const TPZNSMemory &copy);

    TPZNSMemory &operator=(const TPZNSMemory &cp);

    /// Class name
    const std::string Name() const;

    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;

    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);

    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;

    /// Print class attributes
    friend std::ostream & operator<<( std::ostream& out, const TPZNSMemory & memory ){
        memory.Print(out);
        return out;
    }

    virtual int ClassId() const;

    /**
     * @defgroup Set and Get methods
     * @{
     */
    
    /** @brief Set total flux */
    void Set_u(TPZManVector<REAL,3> &u){
        fu = u;
    }
    
    /** @brief Get total flux */
    TPZManVector<REAL,3> u(){
        return fu;
    }

    /** @brief Set total flux at last step */
    void Set_u_last(TPZManVector<REAL,3> &u_last){
        fu_last = u_last;
    }
    
    /** @brief Get total flux at last step */
    TPZManVector<REAL,3> u_last(){
        return fu_last;
    }
    
    /** @brief Set the weighted pressure */
    void Set_p(REAL p){
        fp = p;
    }
    
    /** @brief Get the weighted pressure */
    REAL p(){
        return fp;
    }
    
    /** @brief Set the weighted pressure at the previous timestep */
    void Set_p_n(REAL p_last){
        fp_last = p_last;
    }
    
    /** @brief Get the weighted pressure at the previous timestep */
    REAL p_last(){
        return fp_last;
    }
    
    /** @brief Set the average weighted pressure */
    void Set_p_avg(REAL p_avg){
        fp_avg = p_avg;
    }
    
    /** @brief Get the average weighted pressure */
    REAL p_avg(){
        return fp_avg;
    }
    
    /** @brief Set the average weighted pressure at the previous timestep */
    void Set_p_avg_last(REAL p_avg_last){
        fp_avg_last = p_avg_last;
    }

    
};



#endif /* defined(__PZ__TPZNSMemory__) */

