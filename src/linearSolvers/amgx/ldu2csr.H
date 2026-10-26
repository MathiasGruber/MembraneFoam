/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | openfoam-amgx
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of openfoam-amgx.

    openfoam-amgx is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    openfoam-amgx is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with openfoam-amgx.  If not, see <https://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef membraneLdu2csr_H
#define membraneLdu2csr_H

#include "lduMatrix.H"
#include <vector>

namespace Foam
{
// Serial LDU-to-CSR conversion. Cache only the sparsity and coefficient map;
// matrix values are refreshed on every solve. Duplicate columns are summed.
class ldu2csr
{
    struct Item
    {
        label column;
        scalar value;
        label kind;
        label index;
    };
    using Row = std::vector<Item>;
    labelList rowPtr_, columns_, kinds_, indices_;
    scalarField values_;
    bool reusable_ = false;

    void build(const lduMatrix&);

public:
    explicit ldu2csr(const lduMatrix& matrix) { build(matrix); }
    void updateValues(const lduMatrix&);
    label nRows() const { return rowPtr_.size() - 1; }
    label nNonZeros() const { return values_.size(); }
    const label* rowPtr() const { return rowPtr_.cbegin(); }
    const label* colInd() const { return columns_.cbegin(); }
    scalar* values() { return values_.begin(); }
    const scalar* values() const { return values_.cbegin(); }
};
}
#endif
