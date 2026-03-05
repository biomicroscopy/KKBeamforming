# KKBeamforming

Source code for KK beamforming — an efficient ultrasound beamforming method that operates in the frequency domain using precomputed shift factors and lookup-table-based delay computation. This repository accompanies the associated publication and provides implementations of both KK and conventional Delay-And-Sum (DAS) beamforming for comparison.

## Repository Structure

```
KKBeamforming/
├── SubFunctions/                     Core MATLAB functions
│   ├── initParams.m                  Load RF data and build parameter struct
│   ├── initParamsLUTV2.m             Compute DAS delay lookup tables
│   ├── initKKFreqDomain.m            Compute KK shift factors and delay LUTs
│   ├── VSonicsInit.m                 Extract parameters from Verasonics data
│   ├── BfmDAS.m                      DAS beamforming wrapper
│   ├── BfmKKFreqSum.m                KK beamforming wrapper
│   └── bfmAndProcessFreq.m           Run DAS + multiple KK configurations
│
├── EigenMEX/                         C++ MEX functions and compilation scripts
│   ├── CompressKKFourier.cpp         Frequency-domain KK data compression
│   ├── DASLUT.cpp                    DAS beamforming with delay LUTs
│   ├── KKLUT.cpp                     KK beamforming with delay LUTs
│   ├── CompileAndTestLUTDAS.m        Compile and validate DASLUT
│   ├── CompileScriptBeamformKK.m     Compile and validate KKLUT
│   └── CompileScriptCompressKK.m     Compile and validate CompressKKFourier
│
├── KKTimingExperiment.m              Benchmark KK vs DAS execution time
├── KKfiguresAcquisition.m            Acquisition comparison figures
├── KKfiguresBiologicalData.m         Biological data comparison figures
├── KKfiguresContrastTargets.m        Contrast target analysis figures
├── KKfiguresResolutionTargets.m      Resolution target analysis figures
├── LICENSE                           GPLv3
└── README.md
```

## Dependencies

| Dependency | Version | Notes |
|------------|---------|-------|
| MATLAB | R2018a+ | Required for C++ Data API MEX support |
| MinGW-w64 / GCC | — | C++ compiler with OpenMP support |
| Eigen | 3.4.0 | Header-only linear algebra library |
| Verasonics | — | Ultrasound platform (for data acquisition and format) |

Additional custom MATLAB functions required (not included in this repository):

- `convertParamsToSingle`, `convertParamsToInteger` — type conversion utilities
- `computeNewGrid` — pixel grid resampling
- `computeContrastMatch` — gamma-scale contrast matching
- `plotGammaScaleImage` — B-mode visualization with gamma compression
- `RegionSelector` — ROI definition for contrast metrics
- `computeCRF` — channel RF data reorganization
- `UnscrambleRF` — RF data descrambling

## Compilation

MEX functions must be compiled before use. Each compilation script in `EigenMEX/` handles one MEX function:

1. Open the compilation script in the MATLAB editor
2. Update the Eigen include path (`ipath`) to point to your Eigen 3.4.0 installation
3. Run the script — it will compile the MEX function and run a basic validation test

```matlab
% Example: compile the DAS beamformer
cd EigenMEX
CompileAndTestLUTDAS      % compiles DASLUT.cpp -> DASLUT.mexw64

% Compile KK functions
CompileScriptCompressKK   % compiles CompressKKFourier.cpp
CompileScriptBeamformKK   % compiles KKLUT.cpp
```

The compilation scripts assume a MinGW-w64 compiler on Windows. Compiler flags include `-O3`, `-fopenmp`, `-ffast-math`, and `-march=native`.

## Usage

After compiling the MEX functions, add `SubFunctions/` and `EigenMEX/` to the MATLAB path:

```matlab
addpath('SubFunctions', 'EigenMEX');
```

### Quick start

```matlab
% Load data
dataFile = {'path/to/your/VerasonicsData.mat'};
[p, RFData] = initParams(dataFile, 0);    % filetype 0 = standard VSX export
p.szAcq = int32(p.szRFframe + 1);

% Run DAS + KK beamforming with M=7 RX angles
M = 7;
images = bfmAndProcessFreq(p, RFData, M);

% Display results
figure; imagesc(abs(images(1).data));  title('DAS');
figure; imagesc(abs(images(6).data));  title('KK confocal');
```

### Data paths

The figure and timing scripts locate data relative to a parent `Ultrasound/Datasets/` directory. Update `dataFilePath` at the top of each script to match your local directory structure. Input data must be in Verasonics format (`.mat` files containing `RcvData`, `TX`, `Receive`, `Trans`, `Resource`, etc.).

## Citation

If you use this code in your research, please cite:

```
[CITATION]
```

## License

This project is licensed under the GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.
