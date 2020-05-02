namespace MOSolver2
{
struct SlaterDet {
    std::vector<REAL> coeffs;
    std::vector<REAL> eigenvalues;
};
struct WaveFunction {
    std::vector<GaussBasis> bfs;
    std::vector<SlaterDet>  configuration;
    std::vector<Real>       coeffs;
};

}
