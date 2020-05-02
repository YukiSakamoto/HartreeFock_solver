
namespace MOSolver {

inline bool 
numerical_check(REAL result, REAL answer, REAL delta)
{
    if (std::abs(answer - result) <= delta) {  
        return true;    
    } else {
        return false;
    }
}

}


