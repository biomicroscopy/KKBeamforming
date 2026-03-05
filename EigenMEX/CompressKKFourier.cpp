/*
 * CompressKKFourier - Frequency-domain data compression for KK beamforming
 *
 * MEX function using the MATLAB C++ Data API (R2018a+).
 *
 * Compresses RF data in the frequency domain by applying precomputed
 * shift factors (complex exponential phase terms) to the FFT of the
 * channel data. The result is a compressed dataset indexed by TX and
 * RX angle combinations, ready for beamforming via KKLUT.
 *
 * Syntax (MATLAB):
 *   RFDataKK = CompressKKFourier(param, DataFFT, shiftFac)
 *
 * Inputs:
 *   inputs[0] - param:    int32 vector [szAcq, numEl, na, nRX, midpt]
 *                         Sizes and the one-sided FFT midpoint index.
 *   inputs[1] - DataFFT:  complex single [szAcq x (na*numEl)]
 *                         FFT of the descrambled RF data.
 *   inputs[2] - shiftFac: complex single [(numEl*nRX) x szAcq]
 *                         Precomputed shift factors (phase ramps with
 *                         one-sided Hilbert weighting).
 *
 * Outputs:
 *   outputs[0] - RFDataKK: complex single [szAcq x (na*nRX)]
 *                          Compressed data in the frequency domain.
 *
 * Algorithm:
 *   For each frequency bin t (up to midpt):
 *     1. Extract one row of DataFFT as a [na x numEl] matrix (via stride).
 *     2. Extract the corresponding column of shiftFac as [numEl x nRX].
 *     3. Multiply: RFRow = DataRow * sFacCol  =>  [na x nRX].
 *     4. Write result into the output via strided mapping.
 *   Only the first half of the spectrum is computed (one-sided due to
 *   Hilbert weighting embedded in shiftFac), so the IFFT of the output
 *   produces an analytic signal.
 *
 * Compilation:
 *   See CompileScriptCompressKK.m
 *
 * Dependencies: Eigen 3.4.0, OpenMP
 */
#include <iostream>
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <Eigen/Dense>
#include <MatlabDataArray.hpp>
#include <cmath>
#include <math.h>
#include <omp.h>

class MexFunction : public matlab::mex::Function {
public:
    
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    matlab::data::ArrayFactory factory;
    
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
        // Initialize input parameters
        int tSize = inputs[0][0];
        int xSize = inputs[0][1];
        int TXSize = inputs[0][2];
        int RXSize = inputs[0][3];
        int midpt = inputs[0][4];
        
        auto ptr = getDataPtr<std::complex<float>>(inputs[1]);
        Eigen::Map< const Eigen::MatrixXcf > DataFFT( ptr, tSize, TXSize*xSize );

        
        auto ptr2 = getDataPtr<std::complex<float>>(inputs[2]);
        Eigen::Map< const Eigen::MatrixXcf > shiftFac( ptr2, xSize*RXSize, tSize );

        // Allocate Output
        outputs[0] = factory.createArray<std::complex<float>>({static_cast<size_t>(tSize),static_cast<size_t>(TXSize*RXSize)});
        auto ptrRecon = getOutDataPtr<std::complex<float>>(outputs[0]);
        Eigen::Map<Eigen::MatrixXcf> RFDataKK(ptrRecon,tSize,TXSize*RXSize);
        RFDataKK.setZero();
        
        // Get num threads
        int numThreads = omp_get_max_threads();
        int nProc = omp_get_num_procs();
        omp_set_num_threads(nProc);
        
        
        // Define strided map types for non-contiguous access into column-major data
        using Stride2 = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
        using CMatStrided = Eigen::Map<const Eigen::MatrixXcf, 0, Stride2>;
        using MatStrided  = Eigen::Map<Eigen::MatrixXcf, 0, Stride2>;
        
        // Process each frequency bin independently across threads
        #pragma omp parallel for
        for (int t = 0; t < midpt; ++t) {
            
            // Map to shiftFac col
            Eigen::Map<const Eigen::MatrixXcf> sFacCol(shiftFac.col(t).data(), xSize, RXSize);
            
            // Map row t of DataFFT as [na x numEl] using column-major stride
            CMatStrided DataRow(ptr + t, TXSize, xSize, Stride2(TXSize * tSize, tSize));
            
            // Map row t of output as [na x nRX] using column-major stride
            MatStrided RFRow(ptrRecon + t, TXSize, RXSize, Stride2(TXSize * tSize, tSize));
            
            // Matrix multiply: compress numEl channels into nRX virtual RX angles
            RFRow.noalias() = DataRow * sFacCol;
        }

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