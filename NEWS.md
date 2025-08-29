## HAUDI 1.0.4

- Fixed bug where sample IDs may be read as numeric

## HAUDI 1.0.3

- Fixed bug where no splits argument is available GAUDI
- Removed support for CMSA in GAUDI

## HAUDI 1.0.2

- Fixed bug in add_to_fbm when no columns meet ac threshold

## HAUDI 1.0.1

- Fixed bug with incorrect switch points when converting to .lanc

## HAUDI 1.0.0

### Breaking changes

- Removed direct support for VCF genotype and FLARE local ancestry input
- Renamed several arguments in `haudi` and `gaudi`

### New features

- Create FBM more efficiently using plink2 pgen and Admix-kit .lanc input
- Added efficient c++ functions for converting FLARE and RFMix input to .lanc
- `haudi` now takes a vector of gamma values as input
