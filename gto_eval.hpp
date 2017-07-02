#pragma once
#include "common.hpp"
#include "gto.hpp"
#include "system.hpp"

// Prototypes
REAL overlap_PGTO(const PrimitiveGTO &lhs, const PrimitiveGTO &rhs);
REAL overlap_CGTO(const ContractedGTO &lhs, const ContractedGTO &rhs);
REAL kinetic_PGTO(const PrimitiveGTO &lhs, const PrimitiveGTO &rhs);
REAL kinetic_CGTO(const ContractedGTO &lhs, const ContractedGTO &rhs);
REAL nuclear_attraction_CGTO(const ContractedGTO &lhs, const ContractedGTO &rhs, const Atom &atom);
REAL electron_repulsion_CGTO(const ContractedGTO &cgto1, const ContractedGTO &cgto2, const ContractedGTO &cgto3, const ContractedGTO &cgto4);


