# HAUDI 1.0.0

## Breaking changes

- Removed direct support for VCF genotype and FLARE local ancestry input
- Renamed several arguments in `haudi` and `gaudi`

## New features

- Create FBM more efficiently using plink2 pgen and Admix-kit .lanc input
- Added efficient c++ functions for converting FLARE and RFMix input to .lanc
- `haudi` now takes a vector of gamma values as input
