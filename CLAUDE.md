# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains research work on microtremor array analysis (微動アレイ観測解析) for subsurface structure estimation. The project focuses on implementing and analyzing seismic data using methods like SPAC (Spatial Autocorrelation) and FK (Frequency-Wavenumber) analysis.

## Repository Structure

```
.
├── microtremor_analysis_notes.md  # 詳細な技術ドキュメント
├── microtremor_tutorial.ipynb     # 学生向け実習ノートブック
├── utils/
│   └── microtremor_utils.py      # 解析用ユーティリティ関数
└── data/
    └── sample_microtremor.csv     # サンプルデータ
```

## Key Commands

### Running Jupyter Notebooks
```bash
jupyter notebook microtremor_tutorial.ipynb
```

### Using the utility functions
```python
from utils.microtremor_utils import MicrotremorArray, spac_coefficient, fk_analysis
```

## Code Architecture

### Core Components

1. **MicrotremorArray Class** (`utils/microtremor_utils.py`): Main data container and preprocessing
   - Handles array geometry and coordinate systems
   - Provides preprocessing methods (detrending, filtering, tapering)
   - Computes inter-station distances for SPAC analysis

2. **SPAC Analysis Pipeline**:
   - `spac_coefficient()`: Computes spatial autocorrelation between station pairs
   - `estimate_phase_velocity_spac()`: Extracts phase velocities via Bessel function fitting
   - Quality control through coherence analysis

3. **FK Analysis Pipeline**:
   - `fk_analysis()`: Implements frequency-wavenumber beamforming
   - Returns phase velocities from slowness estimation

4. **Dispersion Curve Processing**:
   - `combine_dispersion_curves()`: Merges results from multiple array sizes
   - `theoretical_dispersion_curve()`: Computes forward model for inversion

### Key Technical Parameters

- **Array radius selection**: Controls analyzable frequency range (f = 100/r to 500/r Hz)
- **Window length**: 20-40 seconds for stable SPAC coefficients
- **Overlap**: 50% for spectral estimation
- **Bessel function fitting range**: c = 50-2000 m/s

## Important Implementation Notes

### SPAC Method Specifics
- Uses J₀ (zeroth-order Bessel function) for circular arrays: ρ(r,f) = J₀(2πrf/c(f))
- Requires sufficient array size relative to wavelength (r < λ/3)
- Sensitive to noise at low frequencies when array is too small

### FK Method Considerations
- Requires regular or well-distributed array geometry
- Computational cost scales with O(n_stations × n_slowness²)
- Resolution depends on array aperture and station distribution

### Data Quality Requirements
- Time synchronization accuracy: < 0.01 seconds (GPS recommended)
- Minimum recording duration: 20-30 minutes for stable statistics
- Coherence threshold: typically > 0.8 for reliable results

## Language and Documentation Conventions

- Code: English variable names and functions
- Documentation: 日本語での説明を優先
- Mathematical notation: LaTeX in markdown (例: $c(f) = f \times \lambda$)
- Jupyter notebooks: 日本語での解説と英語コード
- Markdown and Jupyter notebooks: 日本語の記述言語として使用可能

## Workflow Guidelines

1. データ解析時は必ず品質チェックを実施
2. 新しい解析手法追加時は`microtremor_utils.py`に実装
3. 実装完了時はREADME.mdに更新履歴を記載
4. タスク完了時にgit commitする