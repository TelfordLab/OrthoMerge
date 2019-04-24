# OrthoMerge
Python 3 pipeline to compare and merge orthology inference results of OrthoMCL, OrthoFinder and OrthoInspector

## OrthologyMerge.py


Requirements:
- [BioPython](https://biopython.org/)

Use:
`python3 OrthologyMerge.py -p <phylum> --orthomcl <OrthoMCL_groups> --orthofinder <OrthoFinder_groups> --orthoinspector <OrthoInspector_groups> -fd <FASTA_directory>`

Input:
- orthologous groups output from OrthoMCL, OrthoFinder and OrthoInspector (all in OrthoMCL output format)
- phyla definitions in phyla.txt file in the following format:
`  <phylum>:<species1>,<species2>[...]`
` <superphylum>:<phylum>|<phylum>`
- phylum of interest (must be defined in phyla.txt)

Output:
- log file containing numbers of groups inferred congruently between methods (default name: OrthologousGroups_identical+split.log)
- directory containing merged results of all orthology inference methods (default name: OrthologousGroups_identical+split)

## OG_monophyly_test.py

Requirements:
- [BioPython](https://biopython.org/)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Clustal Omega](http://www.clustal.org/omega/)
- [RAxML-NG](https://github.com/amkozlov/raxml-ng)
- [DendroPy](https://dendropy.org/)

Use:
`python3 OG_monophyly_test.py -od <OrthoGroups_directory> -db <sequence_database>`

Input:
- directory containing orthogroups inferred as congruently inferred by `OrthologyMerge.py`
- database containing non-clade sequences to search for potential homologous sequences outside clade of interest
- *optional*: e-value for BLAST search (`-e <e-value>`, default: 10)

Output:
- text file containing the results of the monophyly test (default name: OrthologousGroups_monophyly.txt)
- directory containing BLAST results for each orthogroup
- directory containing all orthogroups without BLAST hits
- directory containing all orthogroups with BLAST hits including found hits and reconstructed gene trees

## OG_rBBH_test.py

Requirements:
- [BioPython](https://biopython.org/)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [DendroPy](https://dendropy.org/)

Use:
`python3 OG_rBBH_test.py -od <OrthoGroups_directory> -hd <BLAST_hits_directory> -m <Monophyly_test_results> -dbd <sequence_database_directory>`

Input:
- directory containing orthogroups inferred as congruently inferred by `OrthologyMerge.py`
- directory containing all orthogroups with BLAST hits including found hits and reconstructed gene trees
- text file containing the results of the monophyly test
- directory containing BLAST databases for each individual species present in orthogroups

Output:
- text file containing the results of the rBBH test (default name: OrthologousGroups_rBBH.txt)
- directory containing BLAST results for each orthogroup (default name: OrthologousGroups_rBBHits)