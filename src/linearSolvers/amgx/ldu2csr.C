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

#include "ldu2csr.H"
#include <algorithm>
#include "executionPolicy.H"
#include <limits>

void Foam::ldu2csr::build(const lduMatrix& matrix)
{
    const label cells = matrix.diag().size();
    const label faces = matrix.upper().size();
    if (cells <= 0 || faces > (std::numeric_limits<int>::max()-cells)/2)
    {
        FatalErrorInFunction << "AmgX requires a nonempty matrix with int32 CSR offsets"
            << exit(FatalError);
    }
    const auto& upperAddr = matrix.lduAddr().upperAddr();
    const auto& lowerAddr = matrix.lduAddr().lowerAddr();
    if (upperAddr.size() != faces || lowerAddr.size() != faces)
        FatalErrorInFunction << "Inconsistent LDU addressing" << exit(FatalError);

    // One compact temporary array replaces per-row heap allocations. Values
    // remain independent per matrix; only their origin is stored while sorting.
    struct Entry { label column, kind, index; };
    const label count = cells + 2*faces;
    labelList offsets(cells+1, 1);
    offsets[0] = 0;
    for (label f = 0; f < faces; ++f)
    {
        const label lo = lowerAddr[f], hi = upperAddr[f];
        if (lo < 0 || lo >= cells || hi < 0 || hi >= cells)
            FatalErrorInFunction << "LDU address outside matrix" << exit(FatalError);
        ++offsets[lo+1]; ++offsets[hi+1];
    }
    for (label i = 0; i < cells; ++i) offsets[i+1] += offsets[i];
    labelList cursor(offsets);
    std::vector<Entry> entries(count);
    for (label i = 0; i < cells; ++i) entries[cursor[i]++] = {i, 0, i};
    for (label f = 0; f < faces; ++f)
    {
        entries[cursor[upperAddr[f]]++] = {lowerAddr[f], 1, f};
        entries[cursor[lowerAddr[f]]++] = {upperAddr[f], 2, f};
    }
    reusable_ = true;
    rowPtr_.setSize(cells+1); rowPtr_[0] = 0;
    for (label i = 0; i < cells; ++i)
    {
        auto first = entries.begin()+offsets[i];
        const auto last = entries.begin()+offsets[i+1];
        std::sort(first, last,
            [](const Entry& a, const Entry& b) { return a.column < b.column; });
        label unique = 0, previous = -1;
        for (auto p = first; p != last; ++p)
        {
            if (p->column != previous) { ++unique; previous = p->column; }
            else reusable_ = false;
        }
        rowPtr_[i+1] = rowPtr_[i] + unique;
    }
    const label unique = rowPtr_[cells];
    columns_.setSize(unique); values_.setSize(unique);
    kinds_.setSize(unique); indices_.setSize(unique);
    const scalar* diag = matrix.diag().cbegin();
    const scalar* lower = matrix.lower().cbegin();
    const scalar* upper = matrix.upper().cbegin();
    for (label i = 0; i < cells; ++i)
    {
        label j = rowPtr_[i];
        for (label k = offsets[i]; k < offsets[i+1]; ++k)
        {
            const auto& e = entries[k];
            const scalar value = e.kind == 0 ? diag[e.index]
                : e.kind == 1 ? lower[e.index] : upper[e.index];
            if (j > rowPtr_[i] && columns_[j-1] == e.column) values_[j-1] += value;
            else
            {
                columns_[j] = e.column; values_[j] = value;
                kinds_[j] = e.kind; indices_[j] = e.index;
                ++j;
            }
        }
    }
}

void Foam::ldu2csr::updateValues(const lduMatrix& matrix)
{
    if (!reusable_)
    {
        build(matrix);
        return;
    }
    const scalar* diag = matrix.diag().cbegin();
    const scalar* lower = matrix.lower().cbegin();
    const scalar* upper = matrix.upper().cbegin();
    std::transform
    (
        amgxExecution,
        kinds_.cbegin(), kinds_.cend(), indices_.cbegin(), values_.begin(),
        [diag, lower, upper](label kind, label index)
        { return kind == 0 ? diag[index] : kind == 1 ? lower[index] : upper[index]; }
    );
}
