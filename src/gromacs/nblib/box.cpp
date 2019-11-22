//
// Created by sebkelle on 19.11.19.
//

#include "gromacs/utility/exceptions.h"

#include "box.h"

namespace nblib {

Box::Box(real l) : box_{0}
{
    if (std::isnan(l) or std::isinf(l))
    {
        GMX_THROW(gmx::InvalidInputError("Cannot have NaN or Inf box length."));
    }

    box_[XX][XX] = l;
    box_[YY][YY] = l;
    box_[ZZ][ZZ] = l;
}

Box::Box(real x, real y, real z) : box_{0}
{
    if (std::isnan(x) or std::isinf(x) or
        std::isnan(y) or std::isinf(y) or
        std::isnan(z) or std::isinf(z))
    {
        GMX_THROW(gmx::InvalidInputError("Cannot have NaN or Inf box length."));
    }

    box_[XX][XX] = x;
    box_[YY][YY] = y;
    box_[ZZ][ZZ] = z;
}

Box::data_type Box::matrix()
{
    return box_;
}

void convertBoxToGMXFormat(Box &box, matrix matrix)
{
    matrix[XX][XX] = box.matrix()[XX][XX];
    matrix[XX][YY] = box.matrix()[XX][YY];
    matrix[XX][ZZ] = box.matrix()[XX][ZZ];
    matrix[YY][XX] = box.matrix()[YY][XX];
    matrix[YY][YY] = box.matrix()[YY][YY];
    matrix[YY][ZZ] = box.matrix()[YY][ZZ];
    matrix[ZZ][XX] = box.matrix()[ZZ][XX];
    matrix[ZZ][YY] = box.matrix()[ZZ][YY];
    matrix[ZZ][ZZ] = box.matrix()[ZZ][ZZ];
}

} // namespace nblib
