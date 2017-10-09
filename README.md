# MRes
Supplementary material and R source code for my Master of Research thesis, Explaining Imperfect Ant Mimicry.

## Contents
### Supplementary material

* [Table S1.csv](Supplementary%20material/Table%20S1.csv)
  List of references from literature analysis in the introduction. 

* [Table S2.csv](Supplementary%20material/Table%20S2.csv)
  Specimens used in chapter 1 for morphological analysis.

* [Table S3.csv](Supplementary%20material/Table%20S3.csv)
  Species used in chapter 1 for morphological analysis.

* [Table S4.csv](Supplementary%20material/Table%20S4.csv)
  Specimens used in chapter 2 for behavioural analysis.

* [Table S5.csv](Supplementary%20material/Table%20S5.csv)
  Species used in chapter 2 for behavioural analysis.

### Statistical analysis

Contains R source code for analyses from chapters 1 and 2.

[morphometrics.r](Statistical%20analysis/morphometrics.r), function doAnalysis() is the entry point for morphometric analysis, chapter 1. It is automatically run when the file is sourced from the command line.

[motion-analysis.r](Statistical%20analysis/motion-analysis.r), function runMotionAnalysis() is the entry point for behavioural analysis, chapter 2. It is automatically run when the file is sourced from the command line.

Remaining files contain functions invoked from the two above.
