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
class TPZInterfaceMemory : public TPZMonoPhasicMemoryDFN {
  
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
  
    /** @brief Solution left */
    TPZManVector<STATE,3> SolL;
    
    /** @brief Solution right */
    TPZManVector<STATE,3> SolR;
  
    
public:
    
    /** @brief Enum for materials Ids */
    enum MadIds {EMatDarcy = 1, EBCBottom = 2, EBCRight = 3, EBCTop = 4, EBCLeft = 5, EMatFrac = 6, EMatInterFrac = 20, EBCInlet = 7, EBCOutlet = 8 , EBCAuxBottom = 10};
    
    /** @brief Sets the postprocess file for fracture vtk */
    void SetPostProcessFileName(std::string &postProcessFileName);
  
    /** @brief Returns the postprocess file for fracture vtk */
    std::string PostProcessFileName();

    /** @brief Set left solution. */
    void SetLeftSol(TPZManVector<STATE,3> solL){this->SolL = solL;}
    
    /** @brief Returns left solution. */
    TPZManVector<STATE,3> LeftSol() const {return this->SolL;}

    /** @brief Set right solution. */
    void SetRightSol(TPZManVector<STATE,3> solL){this->SolR = solL;}
    
    /** @brief Returns right solution. */
    TPZManVector<STATE,3> RightSol() const {return this->SolR;}
    
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

#endif
