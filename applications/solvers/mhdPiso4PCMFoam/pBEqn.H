for (int Bcorr=0; Bcorr<nBcorr; Bcorr++)
{
    fvVectorMatrix BEqn
    (
        fvm::ddt(B)
      + fvm::div(phi, B)
      - fvc::div(phiB, U)
      - fvm::laplacian(DB, B)
    );

    BEqn.solve();

    volScalarField rBA = 1.0/BEqn.A();

    phiB = (fvc::interpolate(B) & mesh.Sf())
         + fvc::ddtPhiCorr(rBA, B, phiB);

    fvScalarMatrix pBEqn
    (
        fvm::laplacian(rBA, pB) == fvc::div(phiB)
    );
    pBEqn.solve();

    phiB -= pBEqn.flux();

#   include "magneticFieldErr.H"
}
