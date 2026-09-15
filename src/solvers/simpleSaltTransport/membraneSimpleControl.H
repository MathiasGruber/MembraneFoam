// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef membraneSimpleControl_H
#define membraneSimpleControl_H
#include "simpleControl.H"
#include <functional>
#include <utility>

namespace Foam
{
// Residual normalization uses the bulk salt field; a membrane transfers only a
// small fraction of that salt. Require a separate integral conservation check.
class membraneSimpleControl : public simpleControl
{
    std::function<scalar()> saltCheck_;
protected:
    bool criteriaSatisfied() override
    {
        if (!simpleControl::criteriaSatisfied()) return false;
        const scalar limit = dict().getOrDefault<scalar>("maxRelativeSaltImbalance", 0.01);
        if (!(limit > 0 && limit < 1))
            FatalErrorInFunction << "maxRelativeSaltImbalance must be between 0 and 1"
                << exit(FatalError);
        if (!saltCheck_) return false;
        const scalar error = saltCheck_();
        Info<< "MEMBRANE_CONVERGENCE relative_salt_imbalance=" << error
            << " limit=" << limit << nl;
        return error <= limit;
    }
public:
    explicit membraneSimpleControl(fvMesh& mesh) : simpleControl(mesh) {}
    void setSaltCheck(std::function<scalar()> check) { saltCheck_ = std::move(check); }
};
}
#endif
