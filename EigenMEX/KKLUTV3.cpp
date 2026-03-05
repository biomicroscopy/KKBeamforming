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
    int numEl, szAcq, szX, szZ, nPoints, na, nRX;

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
        Eigen::Map< const Eigen::MatrixXcf > RFDataKK( ptr, p.szAcq, p.na*p.nRX );
        
        // Assign LUTs
        auto RXptrX = getDataPtr<float>(inputs[2]);
        Eigen::Map<const Eigen::MatrixXf> RXDelayX(RXptrX, p.szX, p.nRX);
        
        auto RXptrZ = getDataPtr<float>(inputs[3]);
        Eigen::Map<const Eigen::MatrixXf> RXDelayZ(RXptrZ, p.szZ, p.nRX);
        
        auto TXptrX = getDataPtr<float>(inputs[4]);
        Eigen::Map<const Eigen::MatrixXf> TXDelayX(TXptrX, p.szX, p.na);
        
        auto TXptrZ = getDataPtr<float>(inputs[5]);
        Eigen::Map<const Eigen::MatrixXf> TXDelayZ(TXptrZ, p.szZ, p.na);
        
        // Initialize Output
        outputs[0] = factory.createArray<std::complex<float>>({static_cast<size_t>(p.szZ),static_cast<size_t>(p.szX)});
        auto ptrRecon = getOutDataPtr<std::complex<float>>(outputs[0]);
        Eigen::Map<Eigen::MatrixXcf> Recon(ptrRecon,p.szZ,p.szX);

        // Get num threads
        int numThreads = omp_get_max_threads();
        int nProc = omp_get_num_procs();
        omp_set_num_threads(nProc);
        
        int maxIDX = p.szAcq - 2;
        
        #pragma omp parallel for
        for (int n = 0; n < p.nRX; n++) {
            
            for (int j = 0; j < p.na; j++) {
                
                int RFCol = n*p.na + j;
                calcLinKK(p, maxIDX, TXDelayX.col(j), TXDelayZ.col(j), RXDelayX.col(n), 
                        RXDelayZ.col(n), Recon, RFDataKK.col(RFCol));
            }
        }
    }
    
    
    void calcLinKK(const Parameters& p, int& maxIDX,
            const Eigen::Ref<const Eigen::VectorXf> TXDelayX,
            const Eigen::Ref<const Eigen::VectorXf> TXDelayZ,
            const Eigen::Ref<const Eigen::VectorXf> RXDelayX,
            const Eigen::Ref<const Eigen::VectorXf> RXDelayZ,
            Eigen::Ref<Eigen::MatrixXcf> Recon,
            const Eigen::Ref<const Eigen::VectorXcf>& RFDataKK) {
        
        for (int ix = 0; ix < p.szX; ix++) {
            for (int iz = 0; iz < p.szZ; iz++) {
                float idxt = TXDelayX(ix) + TXDelayZ(iz) + RXDelayX(ix) + RXDelayZ(iz);
                int vIDX = static_cast<int>(idxt);
                if (vIDX >= 0 && vIDX < maxIDX) {
                    float interp2 = idxt - vIDX;
                    float interp1 = 1 - interp2;
                    Recon(iz,ix) += (RFDataKK(vIDX)*interp1 + RFDataKK(vIDX+1)*interp2);
                }
            }
        }

    }
    
    void initParams(Parameters& p, matlab::data::StructArray obj) {

        using namespace matlab::data;
        
        // Load and initialize parameters

        TypedArray<int> tempI0 = obj[0]["numEl"];
        TypedArray<int> tempI2 = obj[0]["szAcq"];
        TypedArray<int> tempI3 = obj[0]["szX"];
        TypedArray<int> tempI4 = obj[0]["szZ"];
        TypedArray<int> tempI5 = obj[0]["na"];
        TypedArray<int> tempI7 = obj[0]["nRX"];

        p.numEl = std::move(tempI0[0]);
        p.szAcq = std::move(tempI2[0]);
        p.szX = std::move(tempI3[0]);
        p.szZ = std::move(tempI4[0]);
        p.na = std::move(tempI5[0]);
        p.nRX = std::move(tempI7[0]);

        // Dependent Parameters
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
