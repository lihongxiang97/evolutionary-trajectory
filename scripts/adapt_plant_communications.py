from pathlib import Path
from docx import Document


SRC = Path("work/plant_communications_submission_20260805/Manuscript_T2T_hapA_hapB_3D_mechanisms.docx")
DST = Path("work/plant_communications_submission_20260805/Manuscript_Plant_Communications.docx")


def replace(paragraphs, index, text):
    paragraphs[index].text = text


doc = Document(SRC)
p = doc.paragraphs

# Cell Press/Plant Communications-facing title and single-paragraph Summary.
replace(p, 0, "Duplication-mode-specific 3D chromatin trajectories of duplicate genes in pear")
replace(p, 2, "Summary")
replace(p, 3, "")
replace(
    p,
    4,
    "Gene duplication supplies material for regulatory and functional divergence, but how three-dimensional (3D) genome organization differentiates the evolutionary routes produced by distinct duplication mechanisms remains unclear in plants. We integrated Hi-C, ATAC-seq, RNA-seq, comparative genomics and multi-tissue transcriptomes using the gap-free hapA genome of 'Dangshansuli' pear and used hapB for cross-haplotype validation. The merged hapA map resolved 3,874 topologically associating domains (TADs), 203,233 cis loops and 2,770 trans loops. A two-method consensus classified 1,290 tandem, proximal and transposed duplicate pairs into conserved, neofunctionalized or specialized expression trajectories. Tandem and proximal duplicate pairs occupied the same TAD more often than genomic-distance-matched pairs (observed/expected = 1.09 and 1.25; both BH-adjusted empirical P = 0.0015), whereas transposed duplicates did not. Young transposed duplicates showed lower expression, accessibility, A-compartment occupancy and contact strength after annotation, expression and transposable-element sensitivity analyses. Promoter-distal open-chromatin contacts were balanced between parent and child copies and did not independently predict expression asymmetry. Among 697 strictly mapped hapA-hapB pairs, TAD co-residence agreed for 95.3% and regulatory-capture direction agreed for 82.1%-84.6%. These results distinguish domain co-retention of local duplicates from the age-dependent chromatin states of transposed duplicates."
)
replace(p, 5, "")
replace(p, 6, "")
replace(p, 7, "")
replace(p, 8, "")
replace(p, 13, "Introduction")

# Introduction: field -> unresolved question -> study design and contribution.
replace(
    p,
    14,
    "Three-dimensional genome organization connects chromosome folding with transcriptional regulation and genome evolution. Eukaryotic chromosomes form active and inactive compartments, topological domains and chromatin contacts that constrain regulatory interactions [1-23]. Recent studies have linked enhancer capture to the divergence of distal duplicates in Drosophila [109], associated plant TAD-like domains with tandem-duplicate clusters [110], and connected structural variation with chromatin and regulatory rewiring in soybean and cotton pan-3D genomes [111,112]. These advances raise a broader question: whether distinct duplication mechanisms leave reproducible signatures in plant 3D genome organization."
)
replace(
    p,
    15,
    "Gene duplication provides raw material for dosage retention, expression divergence and functional innovation [24-50]. After small-scale duplication, copies may retain similar expression, one copy may diverge while its partner remains closer to the ancestral state, or both copies may specialize [51-57]. Tandem (TD), proximal (PD) and transposed duplications (TRD) originate through different genomic processes and therefore differ in physical separation, local regulatory context and age distribution [32,40-49]. However, the extent to which these formation modes are coupled to chromatin domains, accessibility and regulatory contacts remains poorly resolved in plants."
)
replace(
    p,
    16,
    "We use 'parent' and 'child' for the inferred source and derived copies of a duplicate pair. Expression from a syntenic one-to-one peach ortholog serves as an outgroup proxy for the ancestral expression state. This design resolves copy-specific expression trajectories while keeping evolutionary polarization distinct from the present-day chromatin measurements made within pear."
)
replace(
    p,
    17,
    "Pear (Pyrus bretschneideri) experienced ancient whole-genome duplication and contains abundant small-scale duplicates [58-64]. We used the gap-free hapA assembly of 'Dangshansuli' pear reported in 'Haplotype-resolved, gap-free genome assemblies provide insights into the divergence between Asian and European pears' [108] as the primary coordinate reference. We integrated published Hi-C, ATAC-seq and RNA-seq data from 15-day-after-flowering fruit [107] with multi-tissue expression and Rosaceae comparative genomics, and used hapB to test coordinate robustness. This framework allowed us to compare expression trajectories, chromatin states, domain co-residence and promoter-distal open-chromatin contacts across duplication modes and gene ages."
)

# Discussion: lead with scientific advance, retain evidence boundaries without process narrative.
replace(
    p,
    84,
    "Our analyses identify duplication mode as an organizing axis that links duplicate-gene expression trajectories with present-day 3D chromatin architecture in pear. Local TD and PD pairs were preferentially contained within the same TAD, whereas TRDs showed the strongest age-dependent differences in expression, accessibility, compartment state and contact strength. Promoter-distal open-chromatin contacts, by contrast, showed no parent-child directional bias and no independent association with expression asymmetry. These complementary results separate domain co-retention, chromatin state and regulatory capture as distinct components of duplicate-gene evolution."
)
replace(
    p,
    85,
    "The primary two-method consensus retained sufficient sample size to compare duplication modes, while CLOUD provided an independent sensitivity classification. The duplication-mode association persisted in the more restrictive three-method set, indicating that the contrast among TD, PD and TRD is not attributable to a single trajectory algorithm. TRDs were enriched for child-neofunctionalized trajectories, whereas TD and PD pairs more often entered specialized states. This distribution is consistent with the distinct formation mechanisms and local genomic contexts of these duplicate classes."
)
replace(
    p,
    86,
    "Expression-sum conservation was enriched among conserved pairs, supporting dosage-related retention as one contributor to stable expression trajectories. Tissue-specificity shifts were concentrated in the diverged copy, while median Ka/Ks values remained below 0.5 across all trajectory-by-duplication-mode groups. Functional divergence in pear therefore combines changes in regulatory breadth with predominantly purifying coding-sequence selection, rather than requiring a genome-wide shift to positive selection."
)
replace(
    p,
    87,
    "Age stratification revealed different chromatin routes for local and transposed duplicates. Young TRD genes consistently showed lower fruit expression, accessibility, A-compartment occupancy and contact strength, and these effects remained after covariate adjustment and quality-control analyses of open reading frames, expression detectability and transposable-element context. TD and PD did not share a uniform young-gene expression pattern. Their strongest structural signature was instead the preferential co-residence of both copies within one TAD. The lower replicate concordance of trans contacts places the trans-loop network outside this cis-focused mechanistic interpretation."
)
replace(
    p,
    89,
    "The enhancer capture-divergence model explains asymmetric distal duplicates in Drosophila [109]. In pear fruit, derived copies did not preferentially acquire promoter-distal open-chromatin contacts, and parent-child capture differences did not predict expression asymmetry after adjustment for promoter state, accessibility and genomic covariates. The 15-day-after-flowering map therefore supports a different regulatory balance: promoter-distal capture is present, but it is not the dominant axis separating parent and child expression at this developmental stage."
)
replace(
    p,
    90,
    "The TD and PD co-residence signal extends the reported association between plant chromatin domains and tandem-duplicate clusters [110]. Distance-matched permutations showed that this enrichment exceeded the expectation from physical proximity alone, although same-TAD status did not independently predict six-tissue coexpression. TADs thus provide a preferential spatial container for local duplicate pairs rather than a universal guarantee of coordinated transcription. The 95.3% hapA-hapB agreement in pair-level TAD status, together with concordant age effects across haplotypes, establishes that the principal conclusions are robust to the two T2T coordinate systems. This developmental map provides a foundation for testing how these relationships vary among pear tissues, stages and populations, complementing emerging pan-3D analyses in other crops [111,112]."
)

# Plant Communications research articles conventionally close in Discussion.
replace(p, 91, "")
replace(p, 92, "")
replace(p, 94, "Materials and methods")

# Journal-facing resource and disclosure headings. Missing author-supplied fields remain blank.
replace(p, 168, "Resource availability")
replace(p, 169, "Data and code availability")
replace(p, 170, p[176].text)
replace(p, 172, "")
replace(p, 173, "")
replace(p, 175, "")
replace(p, 176, "")
replace(p, 178, "Declaration of interests")

doc.core_properties.title = p[0].text
doc.core_properties.subject = "Plant Communications submission version"
doc.save(DST)
print(DST)
