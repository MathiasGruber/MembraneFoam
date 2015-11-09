
// Use vector identity Ei * grad(p) = div( p Ei ) 
// Solve as being equal to mu * laplacian( U ) for getting the real pressure

surfaceScalarField sf = fvc::interpolate( UnityVectorField ) & mesh.Sf(); 

volScalarField RightHand = UnityVectorField & mu_laplacianU;

tmp<fvScalarMatrix> rpEqn 
( 
    fvm::div(sf, rP) - (UnityVectorField & (mu * fvc::laplacian( U )))
); 
rpEqn().relax();
// solve(rpEqn);
