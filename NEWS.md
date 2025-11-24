## HAUDI 1.0.9

- Fix bugs in ancestry queries (wrong variant indices and incorrect sample subsetting)

## HAUDI 1.0.8

- Now use flattened structure to represent ancestry tracts

## HAUDI 1.0.7

- In 'convert_to_lanc', now extend final tract index to the end of the chromosome

## HAUDI 1.0.6

- Now does not throw error in make_fbm if a pvar file does not contain a matching variant

## HAUDI 1.0.5

- Now read psam files as character to avoid big integer issues with sample IDs
- Previously read with automatic colClasses then converted

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
