#include "fluxModel.H"
#include <cassert>
#include <iostream>
int main()
{
    const double a=0.44/3.6e11, b=0.087/3.6e6, k=0.72e6, phi=80510000;
    for (bool advanced : {false,true})
    {
        assert(membrane::flux(.01,.01,a,b,k,phi,advanced)==0);
        for (double feed : {0.,.001,.01}) for (double draw : {.01,.03,.08})
        {
            const auto j = membrane::flux(feed,draw,a,b,k,phi,advanced);
            assert(j>=0 && j<=a*phi*(draw-feed));
            assert(std::abs(membrane::residual(j,feed,draw,a,b,k,phi,advanced))<1e-12);
            const auto noIcp = membrane::flux(feed,draw,a,b,0,phi,advanced);
            assert(std::abs(noIcp-a*phi*(draw-feed))<1e-13);
            assert(std::abs(membrane::flux(feed,draw,a,0,k,phi,true)
                           -membrane::flux(feed,draw,a,0,k,phi,false))<1e-14);
        }
    }
    // Independent long-double bisection oracle over weak and strong ICP,
    // vanishing salt permeability, and nearly equal concentrations.
    for (bool advanced : {false,true})
    for (double permeability : {1e-13,a,1.61111e-12,1e-10})
    for (double salt : {0.,b,8.33333e-8,1e-5})
    for (double resistance : {0.,1.,1e2,150666.,k,6.64e6,1e9,1e12})
    for (double feed : {0.,1e-12,.00065,.01,.08})
    for (double draw : {.01,.065,.09})
    {
        if (draw < feed) continue;
        const double j = membrane::flux(feed,draw,permeability,salt,resistance,phi,advanced);
        long double lo=0, hi=(long double)permeability*phi*(draw-feed);
        const long double solute=advanced ? salt : 0;
        for (int i=0; i<200; ++i)
        {
            const long double x=(lo+hi)/2;
            const long double f=x+solute+(long double)permeability*phi*feed
                -(solute+(long double)permeability*phi*draw)*std::exp(-x*resistance);
            if (f > 0) hi=x;
            else lo=x;
        }
        const double expected=(lo+hi)/2;
        assert(std::abs(j-expected) <= 1e-16+1e-12*expected);
        assert(j>=0 && j<=permeability*phi*(draw-feed));
    }
    bool rejected=false;
    try { membrane::flux(.1,0,a,b,k,phi,true); }
    catch (const std::domain_error&) { rejected=true; }
    assert(rejected);
    std::cout << "FO residual, bounds, equal-concentration, K=0, B=0 and invalid-domain checks passed\n";
}
