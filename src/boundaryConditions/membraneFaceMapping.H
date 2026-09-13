// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef membraneFaceMapping_H
#define membraneFaceMapping_H
#include "fvPatch.H"
#include "labelList.H"
#include <array>
#include <map>
#include <vector>
#include <cmath>

namespace Foam
{
// Coincident, oppositely oriented baffle faces must remain on the same rank.
// Spatial bins avoid the legacy quadratic all-pairs search. Check neighbouring
// bins too, so matching does not depend on rounding at a bin boundary.
inline void membraneFaceMapping(const fvPatch& patch, labelList& mapping)
{
    const scalar tol = 1e-9;
    using Key = std::array<long long,3>;
    std::map<Key,std::vector<label>> bins;
    const vectorField& centres = patch.Cf();
    const auto normals = patch.nf();
    mapping.setSize(patch.size());
    mapping = -1;
    auto key = [tol](const vector& p)
    {
        return Key{{static_cast<long long>(std::floor(p.x()/tol)),
                    static_cast<long long>(std::floor(p.y()/tol)),
                    static_cast<long long>(std::floor(p.z()/tol))}};
    };
    forAll(centres,i) bins[key(centres[i])].push_back(i);
    forAll(centres,i)
    {
        const auto k = key(centres[i]);
        label matches = 0;
        for (int x=-1;x<=1;++x) for (int y=-1;y<=1;++y) for (int z=-1;z<=1;++z)
        {
            const auto found = bins.find(Key{{k[0]+x,k[1]+y,k[2]+z}});
            if (found == bins.end()) continue;
            for (const label j : found->second)
                if (i != j && mag(centres[i]-centres[j]) < tol)
                { mapping[i] = j; ++matches; }
        }
        if (matches != 1)
            FatalErrorInFunction << "Membrane patch " << patch.name()
                << " face " << i << " has " << matches << " coincident partners. "
                << "Use a two-sided conformal baffle and keep each pair on the same MPI rank."
                << exit(FatalError);
        const label j = mapping[i];
        if ((normals()[i] & normals()[j]) > -1+1e-8
            || mag(patch.magSf()[i]-patch.magSf()[j]) > 1e-8*patch.magSf()[i])
            FatalErrorInFunction << "Membrane partners need opposite normals and equal areas"
                << exit(FatalError);
    }
}
}
#endif
