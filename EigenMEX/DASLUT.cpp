/*
 * DASLUT - Delay-And-Sum beamforming using lookup tables
 *
 * MEX function using the MATLAB C++ Data API (R2018a+).
 *
 * Performs DAS beamforming on analytic (Hilbert-transformed) RF data
 * using precomputed TX and RX delay lookup tables. Delays are factored
 * into separable lateral (X) and axial (Z) components for TX, and a
 * combined RX delay table with reduced redundancy.
 *
 * Syntax (MATLAB):
 *   Recon = DASLUT(p, RFData, RXDelay, TXDelayX, TXDelayZ)
 *
 * Inputs:
 *   inputs[0] - p:        struct with fields:
 *                            numEl  (int32) - number of elements
 *                            szRF   (int32) - total RF samples per frame
 *                            szAcq  (int32) - samples per acquisition
 *                            szX    (int32) - number of lateral pixels
 *                            szZ    (int32) - number of axial pixels
 *                            na     (int32) - number of TX angles
 *                            nc     (int32) - number of channels
 *                            ConnMap(int32) - element-to-channel mapping
 *   inputs[1] - RFData:   complex single [szRF x nc]
 *                          Analytic RF data (after Hilbert transform).
 *   inputs[2] - RXDelay:  single [(numEl+szX-1) x szZ]
 *                          RX delay LUT with reduced redundancy (samples).
 *   inputs[3] - TXDelayX: single [szX x na]
 *                          TX lateral delay component (samples).
 *   inputs[4] - TXDelayZ: single [szZ x na]
 *                          TX axial delay component (samples).
 *
 * Outputs:
 *   outputs[0] - Recon: complex single [szZ x szX]
 *                        Beamformed image (coherent sum over elements
 *                        and TX angles).
 *
 * Algorithm:
 *   For each TX angle:
 *     For each pixel (ix, iz):
 *       total_delay = TXDelayX(ix) + TXDelayZ(iz) + RXDelay(iD, iz)
 *       Accumulate linearly interpolated RF samples across all elements.
 *   The RX delay table exploits Toeplitz-like structure: the index
 *   iD = element + (szX - ix) maps lateral pixel position into the
 *   reduced-redundancy table.
 *
 * Compilation:
 *   See CompileAndTestLUTDAS.m
 *
 * Dependencies: Eigen 3.4.0, OpenMP
 */
#include <iostream>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <math.h>
#include <omp.h>

// Parameters Class Definition
class Parameters {
public:
    int numEl, szRF, szAcq, szX, szZ, nPoints, na, nc;
    
    Eigen::VectorXi ConnMap;

    // Constructor
    Parameters() = default;
};

class MexFunction : public matlab::mex::Function {
public:
    
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    matlab::data::ArrayFactory factory;
    
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
        // Load Main Parameters
        matlab::data::StructArray inStructArray = inputs[0];
        
        Parameters p;
        initParams(p, inStructArray);
        
        // Assign input
        auto ptr = getDataPtr<std::complex<float>>(inputs[1]);
        Eigen::Map< const Eigen::MatrixXcf > RFData( ptr, p.szRF, p.nc );
        
        // Assign LUTs
        auto RXptrX = getDataPtr<float>(inputs[2]);
        Eigen::Map<const Eigen::MatrixXf> RXDelay(RXptrX, p.numEl+p.szX-1, p.szZ);
        
        auto TXptrX = getDataPtr<float>(inputs[3]);
        Eigen::Map<const Eigen::MatrixXf> TXDelayX(TXptrX, p.szX, p.na);
        
        auto TXptrZ = getDataPtr<float>(inputs[4]);
        Eigen::Map<const Eigen::MatrixXf> TXDelayZ(TXptrZ, p.szZ, p.na);
        
        // Initialize Output
        outputs[0] = factory.createArray<std::complex<float>>({static_cast<size_t>(p.szZ),static_cast<size_t>(p.szX)});
        auto ptrRecon = getOutDataPtr<std::complex<float>>(outputs[0]);
        Eigen::Map<Eigen::MatrixXcf> Recon(ptrRecon,p.szZ,p.szX);
        Recon.setZero();

        // Get num threads
        int numThreads = omp_get_max_threads();
        int nProc = omp_get_num_procs();
        omp_set_num_threads(nProc);


        for(int i = 0; i < p.na; i++) {
            bfmDASCMPLX(p, TXDelayX.col(i), TXDelayZ.col(i), RXDelay, Recon,
                RFData(Eigen::seq(i*p.szAcq,(i+1)*(p.szAcq)-1),p.ConnMap.array()-1));
        }
        
    }
    
    
    void bfmDASCMPLX(const Parameters& p,
            const Eigen::Ref<const Eigen::VectorXf> TXDelayX,
            const Eigen::Ref<const Eigen::VectorXf> TXDelayZ,
            const Eigen::Ref<const Eigen::MatrixXf> RXDelay,
            Eigen::Ref<Eigen::MatrixXcf> Recon,
            const Eigen::Ref<const Eigen::MatrixXcf>& RFData) {

        // Parallelize over lateral pixel positions
        #pragma omp parallel for
        for (int ix = 0; ix < p.szX; ix++) {
            for (int iz = 0; iz < p.szZ; iz++) {
                // Combined TX delay for this pixel and angle
                float tx = TXDelayX(ix) + TXDelayZ(iz);

                // Map lateral pixel to reduced-redundancy RX delay index
                int ixD = p.szX-ix;
                for (int ie = 0; ie < p.numEl; ++ie) {

                    const int iD = ie + ixD - 1;

                    // Total round-trip delay in samples
                    const float s = RXDelay(iD,iz) + tx;

                    // Linear interpolation between adjacent samples
                    int base = static_cast<int>(s);
                    const float w = s - static_cast<float>(base);

                    // Need base and base+1 in bounds
                    if (base >= 0 && (base + 1) < p.szAcq) {
                        const std::complex<float> v0 = RFData(base,     ie);
                        const std::complex<float> v1 = RFData(base + 1, ie);

                        Recon(iz,ix) += (1.0f - w) * v0 + w * v1;
                    }
                }
            }
        }

    }
    
    
    void initParams(Parameters& p, matlab::data::StructArray obj) {

        using namespace matlab::data;
        
        // Load and initialize parameters

        TypedArray<int> tempI0 = obj[0]["numEl"];
        TypedArray<int> tempI1 = obj[0]["szRF"];
        TypedArray<int> tempI2 = obj[0]["szAcq"];
        TypedArray<int> tempI3 = obj[0]["szX"];
        TypedArray<int> tempI4 = obj[0]["szZ"];
        TypedArray<int> tempI5 = obj[0]["na"];
        TypedArray<int> tempI6 = obj[0]["nc"];

        p.numEl = std::move(tempI0[0]);
        p.szRF = std::move(tempI1[0]);
        p.szAcq = std::move(tempI2[0]);
        p.szX = std::move(tempI3[0]);
        p.szZ = std::move(tempI4[0]);
        p.na = std::move(tempI5[0]);
        p.nc = std::move(tempI6[0]);
        
        matlab::data::TypedArray<int> temp = obj[0]["ConnMap"];
        p.ConnMap = Eigen::Map<Eigen::VectorXi>(temp.release().get(), p.numEl);
        
        p.nPoints = p.szX * p.szZ;
    }
    
    template <typename T>
    const T* getDataPtr(matlab::data::Array arr) {
        const matlab::data::TypedArray<T> arr_t = arr;
        matlab::data::TypedIterator<const T> it(arr_t.begin());
        return it.operator->();
    }

    template <typename T>
    T* getOutDataPtr(matlab::data::Array& arr) {
        auto range = matlab::data::getWritableElements<T>(arr);
        return range.begin().operator->();
    }
};
