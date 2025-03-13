# Functional-Genomics
RNAseq codes
Experimental Design: 


Paragraph Description of Sample Selection Rules:
For this study, we aim to explore the role of the IIS molecular network in body size evolution by analyzing gene expression differences in liver tissue between large and small dog breeds using RNA sequencing (RNAseq) data from the NCBI Sequence Read Archive (SRA). To ensure a robust and reproducible analysis, we established the following sample selection criteria: (1) Tissue specificity: Only liver tissue samples were included to maintain a consistent biological context relevant to metabolic regulation by IIS. (2) Health status: We restricted our selection to healthy control samples, excluding those from diseased dogs or experimental treatments to avoid confounding gene expression changes unrelated to size. (3) Size classification: Breeds were classified as large (>50 lbs or >20 inches at the shoulder) or small (<20 lbs or <15 inches at the shoulder) based on standard adult metrics from veterinary and kennel club data. (4) Diversity: To prevent bias from over-representation, we limited samples to no more than two per breed per size group. (5) Sample size: A minimum of four samples per group (large and small) was targeted to provide sufficient statistical power for differential expression analysis using DESeq2. These rules aim to isolate size-related differences in IIS gene expression while minimizing extraneous variability.
Technical Noise
To control technical noise, we prioritized RNAseq samples sequenced on Illumina platforms to ensure consistency in sequencing technology. Where possible, samples were selected from the same BioProject to reduce batch effects; however, when multiple BioProjects were necessary, we verified that library preparation methods (e.g., poly-A enrichment) and read lengths (e.g., paired-end 100 bp) were comparable. All samples will be aligned to the Tasha the Boxer reference genome (GCF_000002285.5) to standardize mapping and further minimize technical variability. Batch effects will be assessed and corrected during DESeq2 analysis if significant.
Biological Noise: Age, Sex, Location/Region
Biological noise was minimized by controlling for key variables:
Age: Samples were restricted to Adult dogs to reduce gene expression variability due to developmental or aging effects.
Sex: We aimed for only males to reduce sex effect.
Location/Region: While geographic origin was not always specified in SRA metadata, we assumed samples from different BioProjects might reflect regional variation. This was mitigated by focusing on breed and size as the primary variables, with any residual variation modeled as a covariate in DESeq2 if needed.




Parameters for Large and Small Dogs
Large Dogs: Breeds with an average adult weight >50 lbs or height >20 inches at the shoulder. Examples: Newfoundland, Belgian Malanois, Labrador Retriever, Tibetan Mastiff.
Small Dogs: Breeds with an average adult weight <20 lbs or height <15 inches at the shoulder. Examples: Yorkshire Terrier, Toy Poodle, Shih Tzu, West highland white terrier. 
Breeds with ambiguous size classifications (e.g., Beagle) were excluded to maintain a clear size distinction.
Summary Table of RNAseq Samples
Below is a table of RNAseq samples selected for this study, sourced from BioProjects in the NCBI SRA database. All are liver tissue samples from healthy control dogs.


| SRR | Project | Breed | Size | Class | Sex | Age | Group | Tissue |
|-----|---------|-------|------|-------|-----|-----|-------|--------|
|SRR8996966 | PRJNA396033 | Newfoundland | Large | Male | 11 | Liver|






Project Management Plan
Task Assignments
All group members will contribute to coding, method development, and result interpretation, but primary responsibilities are assigned as follows:


Run the Standard Bioinformatic Pipeline:
Primary: ?
Support: All members will review pipeline outputs (e.g., FASTQC, alignment rates).


Run Differential Gene Expression Using DESeq2:
Primary: ?
Support: All members will validate differentially expressed genes and statistical significance.


Run GSEA for KEGG Pathways:
Primary: Ava
Support: All members will interpret enriched pathways biologically.


Run GSEA for IIS Molecular Pathway:
Primary: Ava
Support: All members will confirm IIS gene set accuracy and relevance.


Evaluate Protein-Coding Sequence Variation:
Primary: ?
Support: All members will analyze variants in top IIS regulators (e.g., using VCF files or Ensembl).


Manuscript Organization:
Lead: ?
Support: All members will draft their sections and review the final manuscript.


Presentation Organization:
Lead: ?
Support: All members will create slides and rehearse the presentation.


GitHub and Readme Management:
Lead: Anne Marie (GitHub: https://github.com/amk0121/LIVER_BIOL6850)
Support: All members will upload annotated code and update the README.


Timeline for Project Completion
March 10-15: Finalize sample selection, index the Tasha genome, and run the bioinformatic pipeline (e.g., trimming, alignment, counting).
March 16-19: Generate count matrix and perform DESeq2 differential expression analysis.
March 20-25: Run GSEA for KEGG pathways and IIS pathway.
March 26-31: Assess protein-coding sequence variation in top IIS genes.
April 1: Group check-in to review progress and troubleshoot.
April 2-15: Draft manuscript sections (introduction, methods, results, discussion).
April 16-17: Finalize manuscript draft for peer review.
April 18-21: Conduct peer reviews and revise manuscript.
April 22-30: Prepare and rehearse presentation.
May 1: Deliver group presentation.
May 2-5: Finalize manuscript and GitHub repository.
May 5, 5 PM: Submit final manuscript.
May 7: Complete group reflection survey.
