#ifndef TPZInterfaceMemory_h
#define TPZInterfaceMemory_h

#include "pzvec.h"
#include "pzfmatrix.h"
#include "TPZMonoPhasicMemoryDFN.h"

/**
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief class to store data of fracture propagation simulation
 */
class TPZInterfaceMemory {//}: public TPZMonoPhasicMemoryDFN {
  
public:
  
    /** @brief Default Constructor */
    TPZInterfaceMemory();
  
    /** @brief Destructor */
    ~TPZInterfaceMemory();
  
    /** @brief copy constructor */
    TPZInterfaceMemory(const TPZInterfaceMemory &copy);
  
    /** @brief operator equal */
    TPZInterfaceMemory &operator=(const TPZInterfaceMemory &copy);
  
private:
  
    /** @brief State: Stiffness or Mass Matrix Calculations */
    enum nState { n = 0, nplusone = 1 };
    nState State;
  
    /** @brief Solution left at current state */
    TPZManVector<STATE,3> m_forceFrac_n;
    
    /** @brief Solution right at current state */
    TPZManVector<STATE,3> m_u_n;

//    /** @brief Solution left at last state */
//    TPZManVector<STATE,3> SolL_n;
//
//    /** @brief Solution right at last state */
//    TPZManVector<STATE,3> SolR_n;
//
    
public:
    
    /** @brief Enum for materials Ids */
    enum MadIds {EMatDarcy = 1, EBCBottom = 2, EBCRight = 3, EBCTop = 4, EBCLeft = 5, EMatFrac = 6, EMatInterFrac = 20, EBCInlet = 7, EBCOutlet = 8 , EBCAuxBottom = 10};
    
    /** @brief Sets the postprocess file for fracture vtk */
    void SetPostProcessFileName(std::string &postProcessFileName);
  
    /** @brief Returns the postprocess file for fracture vtk */
    std::string PostProcessFileName();
    
    /** @brief Set left solution at current state. */
    void SetForceFrac_n(TPZManVector<STATE,3> solL){
        m_forceFrac_n = solL;
    }
    
    /** @brief Returns left solution at current state. */
    TPZManVector<STATE,3> GetForceFrac_n() const {
        return m_forceFrac_n;
    }

    /** @brief Set right solution at current state. */
    void Setu_n(TPZManVector<STATE,3> solR){
        m_u_n = solR;
    }
    
    /** @brief Returns right solution at current state. */
    TPZManVector<STATE,3> Getu_n() const {
        return m_u_n;
    }

    /// Class name
    const std::string Name() const;
    
    /// Write class attributes
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /// Read class attributes
    virtual void Read(TPZStream &buf, void *context);
    
    /// Print class attributes
    virtual void Print(std::ostream &out = std::cout) const;
    
    virtual int ClassId() const;
    
    
};

inline std::ostream &operator<<(std::ostream &out, const TPZInterfaceMemory&obj)
{
    obj.Print(out);
    return out;
}


#endif
