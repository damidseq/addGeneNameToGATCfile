# addGeneNameToGATCfile
Add gene name to a RDamID GATC site file 

This R script can be used to add gene name to a RDamID GATC site file.


From:
- "chr" "start.position" "gfpDamS1" "gfpDamS2" "TS1" "TS2"
- "1" "chrI" 2492 0 0 0 0
- "2" "chrI" 2675 0 0 0 4
- "3" "chrI" 3681 7 13 19 13

To:
- "chr" "start.position"  "gfpDamS1"  "gfpDamS2"  "TS1" "TS2" "inGene"
- "1" "chrI"  2492  0 0 0 0 "0"
- "2" "chrI"  2675  0 0 0 4 "0"
- "3" "chrI"  3681  7 13 19 13 "0"
- "4" "chrI"  3750  2 12  2 2 "WBGene00023193"
- "5" "chrI"  3975  1  9  5  0  "0"
- "6" "chrI"  5215  0 8 7 0 "WBGene00022277"
