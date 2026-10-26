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

#ifndef membraneAmgXExecutionPolicy_H
#define membraneAmgXExecutionPolicy_H
#include <execution>
namespace Foam
{
// CPU OpenFOAM builds retain serial conversion without a TBB dependency.
#if defined(__NVCOMPILER_STDPAR_GPU)
inline constexpr auto amgxExecution = std::execution::par_unseq;
#else
inline constexpr auto amgxExecution = std::execution::seq;
#endif
}
#endif
