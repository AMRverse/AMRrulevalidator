# Current version of the AMR rules specification
SPEC_VERSION = "v1.1"

# for the current spec version, this are the columns that are expected in the rules file
CANONICAL_COLUMNS = ["ruleID", "txid", "organism", "gene", "nodeID", "protein accession", "HMM accession", "nucleotide accession", "ARO accession", "gene family", "mutation", "variation type", "gene context", "drug", "drug class", "resistance mechanism", "phenotype", "clinical category", "breakpoint", "breakpoint standard", "breakpoint condition", "PMID", "evidence code", "evidence grade", "evidence limitations", "rule curation note", "version of last update"]

# the following variables are the allowed values for columns in the rules file, where these are not drawn from outside sources like NCBI or CARD
# NCBI/CARD values are included in the ResourceManager class, inside resources.py
PHENOTYPE = ['wildtype', 'nonwildtype']

GENE_CONTEXT = ['core', 'acquired']

CLINICAL_CAT = ["S", "I", "R"]

BREAKPOINT_CONDITIONS = [
            "-", "Endocarditis", "Endocarditis with combination treatment",
            "Intravenous", "Meningitis", "Meningitis, Endocarditis",
            "Non-endocarditis", "Non-meningitis", "Non-meningitis, Non-endocarditis",
            "Non-pneumonia", "Oral", "Oral, Infections originating from the urinary tract",
            "Oral, Other indications", "Oral, Uncomplicated urinary tract infection",
            "Pneumonia", "Prophylaxis", "Respiratory", "Screen", "Skin",
            "Uncomplicated urinary tract infection"
        ]

EVIDENCE_CODES = ["ECO:0001091 knockout phenotypic evidence", "ECO:0000012 functional complementation evidence", "ECO:0001113 point mutation phenotypic evidence", "ECO:0000024 protein-binding evidence", "ECO:0001034 crystallography evidence", "ECO:0000005 enzymatic activity assay evidence", "ECO:0000042 gain-of-function mutant phenotypic evidence", "ECO:0007000 high throughput mutant phenotypic evidence", "ECO:0001103 natural variation mutant evidence", "ECO:0005027 genetic transformation evidence", "ECO:0000020 protein inhibition evidence", "ECO:0006404 experimentally evolved mutant phenotypic evidence", "ECO:0000054 double mutant phenotype evidence", "ECO:0000154 heterologous protein expression evidence", "ECO:0000006 experimental evidence", "ECO:0001583 small interfering RNA knockdown evidence", "ECO:0005528 site-directed mutagenesis phenotypic evidence", "ECO:0006264 competitive binding evidence", "ECO:0006053 pharmacological assay evidence"]

EVIDENCE_GRADES = ["high", "moderate", "low", "very low", "insufficient evidence"]

EVIDENCE_LIMITATIONS = ["lacks evidence for this species", "lacks evidence for this genus", "lacks evidence for this allele", "lacks evidence of the degree to which MIC is affected", "low clinical relevance", "unknown clinical relevance", "statistical geno/pheno evidence but no experimental evidence", "conflicting evidence", "lacks formal breakpoints", "lacks evidence for this drug", "lacks evidence across diverse allelic backgrounds", "lacks evidence across diverse geographies", "limited sources of geno/pheno data"]

# This table is used to determine the call between two rule-based markers for the same drug or drug class. The first key in the dictionary is the call for rule 1. The second key for the dictionary is the call for rule 2, and the value in this second dictionary is the resulting call after comparing the two rules. If the value is IMPOSSIBLE, then the two calls cannot be combined and an error is raised. 
PAIRWISE_TABLE = {
    ('nonwildtype', 'R'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): ('nonwildtype', 'R'),
        ('nonwildtype', 'S'): ('nonwildtype', 'R'), ('nonwildtype', '-'): ('nonwildtype', 'R'),
        ('wildtype', 'R'): 'IMPOSSIBLE',               ('wildtype', 'I'): ('nonwildtype', 'R'),
        ('wildtype', 'S'): ('nonwildtype', 'R'), ('wildtype', '-'): ('nonwildtype', 'R'),
        ('-', '-'): ('nonwildtype', 'R'),
    },
    ('nonwildtype', 'I'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): ('nonwildtype', 'I'),
        ('nonwildtype', 'S'): ('nonwildtype', 'I'), ('nonwildtype', '-'): ('nonwildtype', '-'),
        ('wildtype', 'R'): 'IMPOSSIBLE',               ('wildtype', 'I'): 'IMPOSSIBLE',
        ('wildtype', 'S'): ('nonwildtype', 'I'),  ('wildtype', '-'): ('nonwildtype', 'I'),
        ('-', '-'): ('nonwildtype', 'I'),
    },
    ('nonwildtype', 'S'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): ('nonwildtype', 'I'),
        ('nonwildtype', 'S'): ('nonwildtype', 'S'), ('nonwildtype', '-'): ('nonwildtype', '-'),
        ('wildtype', 'R'): 'IMPOSSIBLE',               ('wildtype', 'I'): 'IMPOSSIBLE',
        ('wildtype', 'S'): ('nonwildtype', 'S'), ('wildtype', '-'): ('nonwildtype', 'S'),
        ('-', '-'): ('nonwildtype', 'S'),
    },
    ('nonwildtype', '-'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): ('nonwildtype', '-'),
        ('nonwildtype', 'S'): ('nonwildtype', '-'), ('nonwildtype', '-'): ('nonwildtype', '-'),
        ('wildtype', 'R'): ('wildtype', 'R'),        ('wildtype', 'I'): ('nonwildtype', '-'),
        ('wildtype', 'S'): ('nonwildtype', '-'), ('wildtype', '-'): ('nonwildtype', '-'),
        ('-', '-'): ('nonwildtype', '-'),
    },
    ('wildtype', 'R'): {
        ('nonwildtype', 'R'): 'IMPOSSIBLE', ('nonwildtype', 'I'): 'IMPOSSIBLE',
        ('nonwildtype', 'S'): 'IMPOSSIBLE', ('nonwildtype', '-'): ('wildtype', 'R'),
        ('wildtype', 'R'): ('wildtype', 'R'), ('wildtype', 'I'): ('wildtype', 'R'),
        ('wildtype', 'S'): ('wildtype', 'R'), ('wildtype', '-'): ('wildtype', 'R'),
        ('-', '-'): ('wildtype', 'R'),
    },
    ('wildtype', 'I'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): 'IMPOSSIBLE',
        ('nonwildtype', 'S'): 'IMPOSSIBLE',             ('nonwildtype', '-'): ('nonwildtype', '-'),
        ('wildtype', 'R'): ('wildtype', 'R'),          ('wildtype', 'I'): ('wildtype', 'I'),
        ('wildtype', 'S'): ('wildtype', 'I'), ('wildtype', '-'): ('wildtype', '-'), ('-', '-'): ('wildtype', 'I'),
    },
    ('wildtype', 'S'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): ('nonwildtype', 'I'),
        ('nonwildtype', 'S'): ('nonwildtype', 'S'), ('nonwildtype', '-'): ('nonwildtype', '-'),
        ('wildtype', 'R'): ('wildtype', 'R'),          ('wildtype', 'I'): ('wildtype', 'I'),
        ('wildtype', 'S'): ('wildtype', 'S'), ('wildtype', '-'): ('wildtype', '-'), ('-', '-'): ('wildtype', 'S'),
    },
    ('wildtype', '-'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): ('nonwildtype', 'I'),
        ('nonwildtype', 'S'): ('nonwildtype', 'S'), ('nonwildtype', '-'): ('nonwildtype', '-'),
        ('wildtype', 'R'): ('wildtype', 'R'),          ('wildtype', 'I'): ('wildtype', '-'),
        ('wildtype', 'S'): ('wildtype', '-'), ('wildtype', '-'): ('wildtype', '-'), ('-', '-'): ('wildtype', '-'),
    },
    ('-', '-'): {
        ('nonwildtype', 'R'): ('nonwildtype', 'R'), ('nonwildtype', 'I'): ('nonwildtype', 'I'),
        ('nonwildtype', 'S'): ('nonwildtype', 'S'), ('nonwildtype', '-'): ('nonwildtype', '-'),
        ('wildtype', 'R'): ('wildtype', 'R'),          ('wildtype', 'I'): ('wildtype', 'I'),
        ('wildtype', 'S'): ('wildtype', 'S'), ('wildtype', '-'): ('wildtype', '-'), ('-', '-'): ('-', '-'),
    },
}
