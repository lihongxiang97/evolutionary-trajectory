from docx import Document


path = "work/plant_communications_submission_20260805/Manuscript_Plant_Communications.docx"
doc = Document(path)
p = doc.paragraphs
p[14].text = (
    "Three-dimensional genome organization connects chromosome folding with transcriptional regulation and genome evolution. "
    "Eukaryotic chromosomes form active and inactive compartments, topological domains and chromatin contacts that constrain regulatory interactions [1-23]. "
    "Recent studies linked enhancer capture to the divergence of distal duplicates in Drosophila [109] and associated plant TAD-like domains with tandem-duplicate clusters [110]. "
    "Pan-3D analyses in soybean and cotton further connected structural variation with chromatin and regulatory rewiring [111,112]. "
    "These advances raise a broader question: whether distinct duplication mechanisms leave reproducible signatures in plant 3D genome organization."
)
p[17].text = (
    "Pear (Pyrus bretschneideri) experienced ancient whole-genome duplication and contains abundant small-scale duplicates [58-64]. "
    "Our primary coordinate reference was the gap-free hapA assembly of 'Dangshansuli' pear [108]. "
    "This assembly was reported in 'Haplotype-resolved, gap-free genome assemblies provide insights into the divergence between Asian and European pears.' "
    "We integrated published Hi-C, ATAC-seq and RNA-seq data from 15-day-after-flowering fruit [107] with multi-tissue expression and Rosaceae comparative genomics. "
    "HapB provided an independent coordinate system for robustness tests. "
    "This framework allowed us to compare expression trajectories, chromatin states, domain co-residence and promoter-distal open-chromatin contacts across duplication modes and gene ages."
)
p[87].text = (
    "Age stratification revealed different chromatin routes for local and transposed duplicates. "
    "Young TRD genes consistently showed lower fruit expression, accessibility, A-compartment occupancy and contact strength. "
    "These effects remained after covariate adjustment and quality-control analyses of open reading frames, expression detectability and transposable-element context. "
    "TD and PD did not share a uniform young-gene expression pattern. "
    "Their strongest structural signature was instead the preferential co-residence of both copies within one TAD. "
    "The lower replicate concordance of trans contacts places the trans-loop network outside this cis-focused mechanistic interpretation."
)
doc.save(path)
