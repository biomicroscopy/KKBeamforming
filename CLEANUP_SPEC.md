# KKBeamforming Repository Cleanup Specification

**Purpose:** Prepare the KKBeamforming repository for public release as a paper companion.
**Audience:** Researchers reproducing results from the associated publication.

---

## 1. README.md Overhaul

**Current state:** Single line — "Source code for KK Beamforming."

**Target:** A practical README covering:
- Project title and one-paragraph description (what this repo provides; refer to the method only as "KK beamforming" — do NOT expand the acronym, as the method is unpublished)
- Repository structure (table listing each folder and key files)
- Dependencies (MATLAB version, MinGW/GCC compiler, Eigen 3.4.0, Verasonics platform)
- Compilation instructions (how to compile the MEX files)
- Usage instructions (how to run the figure/timing scripts)
- Citation placeholder: `[CITATION]` block for the associated paper
- License reference (GPLv3)

---

## 2. .gitignore

**Current state:** None.

**Target:** Add a `.gitignore` covering:
- MATLAB artifacts: `*.asv`, `*.mex*`, `*.mlappinstall`, `slprj/`
- C++ build artifacts: `*.o`, `*.obj`, `*.exe`, `*.dll`
- OS files: `Thumbs.db`, `.DS_Store`
- Editor files: `*.swp`, `*.swo`, `*~`

---

## 3. Code Commenting — Full Headers + Inline

Add documentation to every file following this structure:

### 3a. MATLAB Functions (SubFunctions/)

Each `.m` file gets a header block:
```matlab
%% FUNCTIONNAME - Brief description
%
% Purpose:   What this function does and why
% Algorithm: Brief description of the approach
%
% Inputs:
%   param1 - Description (type, dimensions)
%   param2 - Description (type, dimensions)
%
% Outputs:
%   out1   - Description (type, dimensions)
%
% Dependencies: List of called functions/MEX
%
% See also: RELATEDFUNCTION1, RELATEDFUNCTION2
```

Add inline comments for non-obvious logic (signal processing operations, indexing tricks, coordinate transformations).

**Files:**
- `initParams.m` — Document filetype 0 and 13 formats, parameter struct fields, coordinate system
- `initParamsLUTV2.m` — Document LUT computation, delay model, Hilbert prefactor
- `initKKFreqDomain.m` — Document shift factor computation, TX/RX delay model
- `VSonicsInit.m` — Document Verasonics parameter extraction
- `BfmDAS.m` — Document pipeline stages, fix DASLUTV5 → DASLUTV4 call
- `BfmKKFreqSum.m` — Document pipeline stages, compression + beamforming flow
- `bfmAndProcessFreq.m` — Document the 6 image variants, coherent/incoherent combination

### 3b. C++ MEX Functions (EigenMEX/)

These files use the **MATLAB C++ Data API** (`mex.hpp`, `mexAdapter.hpp`, `MexFunction` class),
**not** the legacy C MEX API (`mex.h`, `mexFunction`, `prhs`/`plhs`).

Each `.cpp` file gets a header block:
```cpp
/*
 * FUNCTIONNAME - Brief description
 *
 * MEX function using the MATLAB C++ Data API (R2018a+).
 *
 * Syntax (MATLAB): output = FunctionName(input1, input2, ...)
 *
 * Inputs:
 *   inputs[0] - Description (MATLAB type, dimensions)
 *   inputs[1] - Description (MATLAB type, dimensions)
 *
 * Outputs:
 *   outputs[0] - Description (MATLAB type, dimensions)
 *
 * Algorithm:
 *   Brief description of the computational approach.
 *
 * Compilation:
 *   See CompileScript*.m
 *
 * Dependencies: Eigen 3.4.0, OpenMP
 */
```

Add inline comments for: memory layout decisions, parallelization strategy, interpolation logic, boundary checking.

**Files:**
- `CompressKKFourierV4.cpp` — Document Eigen stride mapping, shift factor application
- `DASLUTV4.cpp` — Document delay model, interpolation, OpenMP strategy
- `KKLUTV3.cpp` — Document KK vs DAS differences, RX angle loop

### 3c. Root-Level Scripts

Each script gets a top-level comment block:
```matlab
%% ScriptName - Brief description
%
% This script generates [figures/timing data] for [paper section].
%
% Required data: Description of expected .mat files
% Required functions: List of SubFunctions/ and MEX dependencies
% Outputs: What figures/data this produces
```

**Files:**
- `KKTimingExperiment.m`
- `KKfiguresAcquisition.m`
- `KKfiguresBiologicalData.m`
- `KKfiguresContrastTargets.m`
- `KKfiguresResolutionTargets.m`

### 3d. Compilation Scripts (EigenMEX/)

Add header explaining purpose, then document which functions are compiled and what flags are used.

**Files:**
- `CompileAndTestLUTDAS.m`
- `CompileScriptBeamformKKV2.m`
- `CompileScriptCompressKKV3.m`

---

## 4. Version Cleanup — Keep Latest Only

### 4a. Remove Old Version References

Strip out all commented-out code referencing old versions (V1, V2, V3 where V4 is current, etc.) from:
- `CompileAndTestLUTDAS.m` — Remove references to DASLUTV2, V3, V5; keep only V4 compile + test
- `CompileScriptBeamformKKV2.m` — Remove references to KKLUT V1, V2; keep only V3 compile + test
- `CompileScriptCompressKKV3.m` — Remove old compression approaches (V1, V2, V3 inline code); keep only V4 compile + test

### 4b. Rename to Clean Names

Remove version suffixes from filenames and function calls:
- `CompressKKFourierV4.cpp` → `CompressKKFourier.cpp`
- `DASLUTV4.cpp` → `DASLUT.cpp`
- `KKLUTV3.cpp` → `KKLUT.cpp`
- `CompileAndTestLUTDAS.m` — Update MEX call target
- `CompileScriptBeamformKKV2.m` → `CompileScriptBeamformKK.m`; update MEX call target
- `CompileScriptCompressKKV3.m` → `CompileScriptCompressKK.m`; update MEX call target
- `BfmDAS.m` — Fix `DASLUTV5` → `DASLUT`
- `BfmKKFreqSum.m` — Update `CompressKKFourierV4` → `CompressKKFourier`, `KKLUTV3` → `KKLUT`
- `bfmAndProcessFreq.m` — Update any version-suffixed calls
- All other files referencing the old names

---

## 5. Code Cleanup

### 5a. Remove Commented-Out Code
- Strip dead/commented-out code blocks throughout all files
- Keep only comments that explain logic (not old implementations)

### 5b. Remove export_fig Calls
- Remove all `export_fig(...)` calls from figure scripts
- Do not replace with alternatives

### 5c. Fix initParams.m
- Remove the TODO comment ("Nikunj will insert documentation here")
- Fix filetype 2: change to filetype 0 (duplicate path, merge/correct)
- Keep only filetype 0 and filetype 13 code paths
- Add documentation for both formats

### 5d. Clean Up Compilation Scripts
- Minimize testing/benchmarking code to only test the current (kept) versions
- Remove timing comparisons of old versions
- Keep a basic validation test for each compiled MEX function

### 5e. Remove Magic Number Opacity
- Where feasible, add brief inline comments explaining hardcoded values in the figure scripts
- Do not restructure the figure scripts themselves

---

## 6. Path Documentation

Hardcoded paths will remain as-is, but each script/function that uses them will have a comment at the top indicating:
- What paths need to be modified by the user
- What the expected directory structure looks like

The README will also contain a section explaining the expected data/dependency layout.

---

## 7. File List Summary

| Action | File |
|--------|------|
| Rewrite | `README.md` |
| Create | `.gitignore` |
| Comment + clean | `SubFunctions/initParams.m` |
| Comment + clean | `SubFunctions/initParamsLUTV2.m` |
| Comment + clean | `SubFunctions/initKKFreqDomain.m` |
| Comment + clean | `SubFunctions/VSonicsInit.m` |
| Comment + fix + rename refs | `SubFunctions/BfmDAS.m` |
| Comment + rename refs | `SubFunctions/BfmKKFreqSum.m` |
| Comment + rename refs | `SubFunctions/bfmAndProcessFreq.m` |
| Comment + rename | `EigenMEX/CompressKKFourierV4.cpp` → `CompressKKFourier.cpp` |
| Comment + rename | `EigenMEX/DASLUTV4.cpp` → `DASLUT.cpp` |
| Comment + rename | `EigenMEX/KKLUTV3.cpp` → `KKLUT.cpp` |
| Clean + rename | `EigenMEX/CompileAndTestLUTDAS.m` |
| Clean + rename | `EigenMEX/CompileScriptBeamformKKV2.m` → `CompileScriptBeamformKK.m` |
| Clean + rename | `EigenMEX/CompileScriptCompressKKV3.m` → `CompileScriptCompressKK.m` |
| Minimal cleanup | `KKTimingExperiment.m` |
| Minimal cleanup | `KKfiguresAcquisition.m` |
| Minimal cleanup | `KKfiguresBiologicalData.m` |
| Minimal cleanup | `KKfiguresContrastTargets.m` |
| Minimal cleanup | `KKfiguresResolutionTargets.m` |
| Delete after rename | Old-named files (handled by git mv) |

---

## 8. Missing MATLAB Dependencies (User Action Required)

The following custom MATLAB functions are called throughout the codebase but are **not included in the repository**.
These must be added by the repository owner — they will NOT be created as part of this cleanup.

| Function | Called in |
|----------|----------|
| `convertParamsToSingle` | `initParamsLUTV2.m`, `initKKFreqDomain.m`, `CompileAndTestLUTDAS.m`, `CompileScriptCompressKKV3.m` |
| `convertParamsToInteger` | `initParamsLUTV2.m`, `initKKFreqDomain.m`, `CompileAndTestLUTDAS.m`, `CompileScriptCompressKKV3.m` |
| `computeNewGrid` | All 4 `KKfigures*.m` scripts |
| `computeContrastMatch` | All 4 `KKfigures*.m` scripts |
| `plotGammaScaleImage` | All 4 `KKfigures*.m`, all 3 `Compile*.m` scripts |
| `RegionSelector` | `KKfiguresAcquisition.m`, `KKfiguresContrastTargets.m` |
| `computeCRF` | `CompileAndTestLUTDAS.m`, `CompileScriptBeamformKKV2.m`, `CompileScriptCompressKKV3.m` |
| `UnscrambleRF` | `CompileAndTestLUTDAS.m` |
| `BfmKKPixelCPP` | `CompileScriptBeamformKKV2.m` |
| `BfmKKFreqMTX` | `CompileScriptBeamformKKV2.m` |
| `reconraw.DASBModeOffline` | `CompileAndTestLUTDAS.m` (external class/package) |

---

## 9. Out of Scope

The following are explicitly **not** part of this cleanup:
- Adding unit tests or a test framework
- Restructuring the figure scripts' layout/positioning logic
- Adding error handling or input validation
- Creating build scripts (CMake, Makefile)
- Including sample data
- Changing the algorithm or computational logic
- Modifying the license
