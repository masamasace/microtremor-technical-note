# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains research work on microtremor array analysis (微動アレイ観測解析) for subsurface structure estimation. The project focuses on implementing and analyzing seismic data using methods like SPAC (Spatial Autocorrelation) and FK (Frequency-Wavenumber) analysis.

## Technical Context

### Key Analysis Methods
- **SPAC Method**: Spatial autocorrelation analysis for phase velocity estimation
- **FK Method**: Frequency-wavenumber analysis for wave propagation
- **Dispersion Curve Analysis**: Surface wave dispersion for subsurface structure inversion

### Important Technical Considerations
- Microtremor frequency range: 0.1-20 Hz
- Exploration depth estimation: λ/2 to λ/3 (where λ is wavelength)
- Data processing requires careful windowing (20-40 second segments)
- Quality control through SNR and coherence analysis

## Development Guidelines

### When implementing analysis code:
1. Follow the theoretical framework outlined in `microtremor_analysis_notes.md`
2. Implement robust error handling for signal processing operations
3. Use appropriate windowing functions (e.g., Hanning) for spectral analysis
4. Include data quality checks (SNR, coherence thresholds)

### Common data formats to support:
- SEG-Y, SAC, MiniSEED for seismic data
- GPS time synchronization is critical

### Typical analysis workflow:
1. Data preprocessing (detrending, filtering)
2. Spectral analysis (FFT with proper windowing)
3. Spatial correlation or FK analysis
4. Dispersion curve extraction
5. Inversion for subsurface structure

## Language Considerations

This is a Japanese research project. When working with:
- Variable names: Use English for code, Japanese comments are acceptable
- Documentation: Maintain Japanese technical terms where appropriate
- File naming: English preferred for code files