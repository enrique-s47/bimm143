Class 11: Structural Bioinformatics
================

The PDB Database
----------------

The \[PDB\] (<http://www.rcsb.org/>) is the main repository for biomolecular structure data.

Here we examine the contents of the PDB:

``` r
db<- read.csv("Data Export Summary.csv", row.names=1)
head(db)
```

    ##                     Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## X-Ray                 126880          2012               6547     8 135447
    ## NMR                    11062          1279                259     8  12608
    ## Electron Microscopy     2277            31                800     0   3108
    ## Other                    256             4                  6    13    279
    ## Multi Method             129             5                  2     1    137

How many are X-ray

``` r
(db$Total/sum(db$Total)) *100
```

    ## [1] 89.35736481  8.31777489  2.05041595  0.18406244  0.09038191

What percent are protein

``` r
sum(db$Proteins/sum(db$Total)) *100
```

    ## [1] 92.75955

We could also try "Addins"&gt; "Paste as data.frame"

``` r
library(datapasta)

tmp <-data.frame(stringsAsFactors=FALSE,
   Experimental.Method = c("X-Ray", "Other", "NMR", "Multi Method",
                           "Electron Microscopy", "Total"),
              Proteins = c(126880, 256, 11062, 129, 2277, 140604),
         Nucleic.Acids = c(2012, 4, 1279, 5, 31, 3331),
    ProteinComplex = c(6547, 6, 259, 2, 800, 7614),
                 Other = c(8, 13, 8, 1, 0, 30),
                 Total = c(135447, 279, 12608, 137, 3108, 151579)
)

tmp
```

    ##   Experimental.Method Proteins Nucleic.Acids ProteinComplex Other  Total
    ## 1               X-Ray   126880          2012           6547     8 135447
    ## 2               Other      256             4              6    13    279
    ## 3                 NMR    11062          1279            259     8  12608
    ## 4        Multi Method      129             5              2     1    137
    ## 5 Electron Microscopy     2277            31            800     0   3108
    ## 6               Total   140604          3331           7614    30 151579

There are 1157 as of 2019-05-07 See: <http://www.rcsb.org/pdb/results/results.do?tabtoshow=Current&qrid=81599ADF>

Section 3 Using Bio3d
---------------------

``` r
library(bio3d)

pdb<- read.pdb("1hsg") 
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

Q6. How many amino acid residues are there in this pdb object and what are the two nonprotein
---------------------------------------------------------------------------------------------

residues? Theres 198 protein residues The two nonprotein residues are HOH (127) and MK1(1)

Q7. What type of R object is pdb$atom? HINT: You can always use the str() function to get a useful summery of any R object.
---------------------------------------------------------------------------------------------------------------------------

``` r
attributes(pdb)
```

    ## $names
    ## [1] "atom"   "xyz"    "seqres" "helix"  "sheet"  "calpha" "remark" "call"  
    ## 
    ## $class
    ## [1] "pdb" "sse"

``` r
str(pdb$atom)
```

    ## 'data.frame':    1686 obs. of  16 variables:
    ##  $ type  : chr  "ATOM" "ATOM" "ATOM" "ATOM" ...
    ##  $ eleno : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ elety : chr  "N" "CA" "C" "O" ...
    ##  $ alt   : chr  NA NA NA NA ...
    ##  $ resid : chr  "PRO" "PRO" "PRO" "PRO" ...
    ##  $ chain : chr  "A" "A" "A" "A" ...
    ##  $ resno : int  1 1 1 1 1 1 1 2 2 2 ...
    ##  $ insert: chr  NA NA NA NA ...
    ##  $ x     : num  29.4 30.3 29.8 28.6 30.5 ...
    ##  $ y     : num  39.7 38.7 38.1 38.3 37.5 ...
    ##  $ z     : num  5.86 5.32 4.02 3.68 6.34 ...
    ##  $ o     : num  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ b     : num  38.1 40.6 42.6 43.4 37.9 ...
    ##  $ segid : chr  NA NA NA NA ...
    ##  $ elesy : chr  "N" "C" "C" "O" ...
    ##  $ charge: chr  NA NA NA NA ...

``` r
pdb$atom$resid
```

    ##    [1] "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "GLN" "GLN" "GLN" "GLN"
    ##   [12] "GLN" "GLN" "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##   [23] "ILE" "ILE" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "LEU" "LEU"
    ##   [34] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "TRP" "TRP" "TRP" "TRP" "TRP"
    ##   [45] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "GLN" "GLN"
    ##   [56] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "ARG" "ARG" "ARG" "ARG"
    ##   [67] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "PRO" "PRO" "PRO" "PRO"
    ##   [78] "PRO" "PRO" "PRO" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##   [89] "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "THR" "THR" "THR" "THR"
    ##  [100] "THR" "THR" "THR" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [111] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ILE" "ILE"
    ##  [122] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY"
    ##  [133] "GLY" "GLY" "GLY" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [144] "GLN" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LYS" "LYS"
    ##  [155] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "GLU" "GLU" "GLU" "GLU"
    ##  [166] "GLU" "GLU" "GLU" "GLU" "GLU" "ALA" "ALA" "ALA" "ALA" "ALA" "LEU"
    ##  [177] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [188] "LEU" "LEU" "LEU" "LEU" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [199] "ASP" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "GLY" "GLY" "GLY"
    ##  [210] "GLY" "ALA" "ALA" "ALA" "ALA" "ALA" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [221] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [232] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "VAL" "VAL" "VAL" "VAL"
    ##  [243] "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [254] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ##  [265] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "MET" "MET" "MET" "MET"
    ##  [276] "MET" "MET" "MET" "MET" "SER" "SER" "SER" "SER" "SER" "SER" "LEU"
    ##  [287] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "PRO" "PRO" "PRO" "PRO"
    ##  [298] "PRO" "PRO" "PRO" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG"
    ##  [309] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "TRP" "TRP" "TRP" "TRP"
    ##  [320] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "LYS"
    ##  [331] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "PRO" "PRO" "PRO"
    ##  [342] "PRO" "PRO" "PRO" "PRO" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ##  [353] "LYS" "LYS" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "ILE"
    ##  [364] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY"
    ##  [375] "GLY" "GLY" "GLY" "GLY" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [386] "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "PHE" "PHE"
    ##  [397] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "ILE" "ILE"
    ##  [408] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "LYS" "LYS" "LYS" "LYS" "LYS"
    ##  [419] "LYS" "LYS" "LYS" "LYS" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ##  [430] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [441] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "TYR" "TYR"
    ##  [452] "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "ASP"
    ##  [463] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "GLN" "GLN" "GLN" "GLN"
    ##  [474] "GLN" "GLN" "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [485] "ILE" "ILE" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ILE"
    ##  [496] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLU" "GLU" "GLU" "GLU"
    ##  [507] "GLU" "GLU" "GLU" "GLU" "GLU" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [518] "ILE" "ILE" "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "GLY" "GLY" "GLY"
    ##  [529] "GLY" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS"
    ##  [540] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ALA" "ALA"
    ##  [551] "ALA" "ALA" "ALA" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [562] "GLY" "GLY" "GLY" "GLY" "THR" "THR" "THR" "THR" "THR" "THR" "THR"
    ##  [573] "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU"
    ##  [584] "LEU" "LEU" "LEU" "LEU" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ##  [595] "GLY" "GLY" "GLY" "GLY" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ##  [606] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "PRO" "PRO" "PRO" "PRO"
    ##  [617] "PRO" "PRO" "PRO" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ASN"
    ##  [628] "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ILE" "ILE" "ILE" "ILE"
    ##  [639] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [650] "ILE" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [661] "ARG" "ARG" "ARG" "ARG" "ARG" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN"
    ##  [672] "ASN" "ASN" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [683] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "THR" "THR" "THR" "THR"
    ##  [694] "THR" "THR" "THR" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [705] "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY"
    ##  [716] "GLY" "GLY" "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "THR" "THR" "THR"
    ##  [727] "THR" "THR" "THR" "THR" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [738] "LEU" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "PHE" "PHE"
    ##  [749] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PRO" "PRO"
    ##  [760] "PRO" "PRO" "PRO" "PRO" "PRO" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ##  [771] "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ##  [782] "THR" "THR" "THR" "THR" "THR" "THR" "THR" "LEU" "LEU" "LEU" "LEU"
    ##  [793] "LEU" "LEU" "LEU" "LEU" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP"
    ##  [804] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "GLN" "GLN" "GLN" "GLN"
    ##  [815] "GLN" "GLN" "GLN" "GLN" "GLN" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ##  [826] "ARG" "ARG" "ARG" "ARG" "ARG" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ##  [837] "PRO" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "VAL" "VAL"
    ##  [848] "VAL" "VAL" "VAL" "VAL" "VAL" "THR" "THR" "THR" "THR" "THR" "THR"
    ##  [859] "THR" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "LYS" "LYS"
    ##  [870] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ILE" "ILE" "ILE" "ILE"
    ##  [881] "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY"
    ##  [892] "GLY" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "LEU"
    ##  [903] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LYS" "LYS" "LYS" "LYS"
    ##  [914] "LYS" "LYS" "LYS" "LYS" "LYS" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ##  [925] "GLU" "GLU" "GLU" "ALA" "ALA" "ALA" "ALA" "ALA" "LEU" "LEU" "LEU"
    ##  [936] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ##  [947] "LEU" "LEU" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "THR"
    ##  [958] "THR" "THR" "THR" "THR" "THR" "THR" "GLY" "GLY" "GLY" "GLY" "ALA"
    ##  [969] "ALA" "ALA" "ALA" "ALA" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP"
    ##  [980] "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "ASP" "THR" "THR"
    ##  [991] "THR" "THR" "THR" "THR" "THR" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL"
    ## [1002] "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "GLU" "GLU"
    ## [1013] "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ## [1024] "GLU" "GLU" "GLU" "GLU" "GLU" "MET" "MET" "MET" "MET" "MET" "MET"
    ## [1035] "MET" "MET" "SER" "SER" "SER" "SER" "SER" "SER" "LEU" "LEU" "LEU"
    ## [1046] "LEU" "LEU" "LEU" "LEU" "LEU" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1057] "PRO" "GLY" "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ## [1068] "ARG" "ARG" "ARG" "ARG" "ARG" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP"
    ## [1079] "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "TRP" "LYS" "LYS" "LYS"
    ## [1090] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1101] "PRO" "PRO" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ## [1112] "MET" "MET" "MET" "MET" "MET" "MET" "MET" "MET" "ILE" "ILE" "ILE"
    ## [1123] "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY"
    ## [1134] "GLY" "GLY" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY"
    ## [1145] "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "GLY" "PHE" "PHE" "PHE" "PHE"
    ## [1156] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "ILE" "ILE" "ILE" "ILE"
    ## [1167] "ILE" "ILE" "ILE" "ILE" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS"
    ## [1178] "LYS" "LYS" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ARG" "ARG"
    ## [1189] "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "GLN" "GLN"
    ## [1200] "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "TYR" "TYR" "TYR" "TYR"
    ## [1211] "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "TYR" "ASP" "ASP" "ASP"
    ## [1222] "ASP" "ASP" "ASP" "ASP" "ASP" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN"
    ## [1233] "GLN" "GLN" "GLN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1244] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ILE" "ILE" "ILE"
    ## [1255] "ILE" "ILE" "ILE" "ILE" "ILE" "GLU" "GLU" "GLU" "GLU" "GLU" "GLU"
    ## [1266] "GLU" "GLU" "GLU" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1277] "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "GLY" "GLY" "GLY" "GLY" "HIS"
    ## [1288] "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "HIS" "LYS" "LYS"
    ## [1299] "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "LYS" "ALA" "ALA" "ALA" "ALA"
    ## [1310] "ALA" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY"
    ## [1321] "GLY" "GLY" "THR" "THR" "THR" "THR" "THR" "THR" "THR" "VAL" "VAL"
    ## [1332] "VAL" "VAL" "VAL" "VAL" "VAL" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ## [1343] "LEU" "LEU" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "GLY" "GLY"
    ## [1354] "GLY" "GLY" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO" "THR" "THR"
    ## [1365] "THR" "THR" "THR" "THR" "THR" "PRO" "PRO" "PRO" "PRO" "PRO" "PRO"
    ## [1376] "PRO" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "VAL" "ASN" "ASN" "ASN"
    ## [1387] "ASN" "ASN" "ASN" "ASN" "ASN" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE"
    ## [1398] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY"
    ## [1409] "GLY" "GLY" "GLY" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG" "ARG"
    ## [1420] "ARG" "ARG" "ARG" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN"
    ## [1431] "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU"
    ## [1442] "LEU" "LEU" "LEU" "LEU" "LEU" "THR" "THR" "THR" "THR" "THR" "THR"
    ## [1453] "THR" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "GLN" "ILE"
    ## [1464] "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "ILE" "GLY" "GLY" "GLY" "GLY"
    ## [1475] "CYS" "CYS" "CYS" "CYS" "CYS" "CYS" "THR" "THR" "THR" "THR" "THR"
    ## [1486] "THR" "THR" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "LEU" "ASN"
    ## [1497] "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "ASN" "PHE" "PHE" "PHE" "PHE"
    ## [1508] "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "PHE" "MK1" "MK1" "MK1" "MK1"
    ## [1519] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1530] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1541] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1"
    ## [1552] "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "MK1" "HOH" "HOH" "HOH"
    ## [1563] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1574] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1585] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1596] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1607] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1618] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1629] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1640] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1651] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1662] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1673] "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH" "HOH"
    ## [1684] "HOH" "HOH" "HOH"

Atom selection is done via the function **atom.select()**

``` r
prot.pdb<- atom.select(pdb, "protein", value= TRUE)
write.pdb(prot.pdb, file="1hsg_protein.pdb")
```

``` r
lig.pdb<- atom.select(pdb, "ligand", value= TRUE)
write.pdb(lig.pdb, file = "1hsg_protein.pdb")
```

Section 5

``` r
aa <- get.seq("1ake_A")
```

    ## Warning in get.seq("1ake_A"): Removing existing file: seqs.fasta

``` r
b <- blast.pdb(aa)
```

    ##  Searching ... please wait (updates every 5 seconds) RID = D3R22SR8015 
    ##  .
    ##  Reporting 97 hits
