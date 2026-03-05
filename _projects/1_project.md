---
layout: page
title: CAFA6 Challenge
description: Predicting a protein's function from its amino acid sequence
img: assets/img/cafa_image.jpeg
importance: 1
category: work
related_publications: false
---

Can we predict a protein's function from its amino acid sequence? The CAFA (Critical Assessment of protein Function Annotation) 6th challenge aims to predict protein function prediction using Gene Ontology (GO) annotations. GO annotations are the function annotations associated with the protein, the GO terms are structured as a Directed Acyclic Graph (DAG), and belong to three non-overlapping ontologies: Molecular Function (MF), Biological Process (BP) and Cellular Component (CC). Each protein can have multiple GO terms, and annotations must respect the hierarchical structure of the ontology (true path rule).

**This task is a multi-label classification problem:**
	
	- Input: Protein amino acid sequence
	- Output: A probability score for each GO term


**Challenges:**

    - Sparce labelling / Class imbalance: most GO terms are rare
    - Hierarchical Label Structure: the labels are associated with each other. If a protein is associated with a GO term, it must also be annotated with that term's ancestor terms.
    - Varying protein sequence length: We must find a reliable way to represent proteins that are 30 to 30,000 amino acids long
    - Large dataset, and limited computational power on my laptop


**Pipeline Outline:**

    - Compute protein embeddings from ESM2 (6 layers and embedding dimension 320). Frozen embeddings used.
	- Design a lightweight MLP to predict the associated GO terms
	- Use a loss function which includes
        1. assymetric BCE (to weigh false positives and false negatives according to GO term)
        2. coverage loss
        3. GO term heirarchy enforcing loss
        4. Negative GO term Loss

---

# Pipeline:

### GO term Processing

We first split the GO terms into the three sub-ontologies (Molecular Function, Biological Process and Cellular Component). We then elimiate any GO terms not present in the training data. Associate each protein id to a set of positive GO terms.

### Information Content (IC) Computation
Two complementary information-content measures were computed:

     1. Structural IC
     Derived purely from GO DAG topology. For each ontology separately, IC was defined as:
    $$
    IC_{\text{struct}}(t) = 1 - \frac{\log_2 \left( \left| \text{descendants}(t) \right| \right)}
    {\log_2 \left( \left| \text{ontology} \right| \right)}
    $$
    This assigns higher IC to more specific (leaf-like) terms and lower IC to generic ancestors.
	
    2. Abundance-Based IC
    Using training annotations, I computed IC based on empirical frequency:
    $$
    IC_{\text{abundance}}(t) = -\log_2 \, p(t)
    $$
     where p(t)accounts for both direct annotations and descendant propagation, using goatools.TermCounts.
	 
Both IC measures were retained and later used for weighting false positives vs false negatives, and selecting informative negative GO terms.


### Negative GO Term Sampling

Standard CAFA training treats all unannotated GO terms as negatives, which is problematic: many “negatives” are actually unknown positives, the loss is dominated by trivial easy negatives, models learn to predict everything as zero. To address this, I explicitly constructed informative negative GO sets per protein. This analysis was done as described in the [NegGOA study](https://academic.oup.com/bioinformatics/article/32/19/2996/2196619).
    
1. Conditional GO Co-Occurrence Matrix
        - Constructing the GO co-occurrence graph/matrix
    $$
    P(B \mid A) = \frac{\#(A, B)}{\#(A)}
    $$
        This captures how often GO term B appears given A. For computational ease, I sparsified it to the top-100 conditional neighbors per GO term, and then row-normalized it to preserve probability semantics.
        - Conditional Random Walk 
        We perform a random walk to the neighbours in the co-occurance graph, in accordance to the weight of co-occurance. We take 4 steps, with a 50% probability of restart. 
        $$
        Rc = alpha * (Rc @ Wc) + (1 - alpha) * I
        $$
        where Wc is the GO Co-occurance matrix, alpha is the restart probability (0.5) and Rc is the matrix of the conditional random walk, which determines the weight of possible co-occurances. The random walk is performed over 4 iteration. 

2. Heirarchial GO Matrix
        - Contruct the heirarchial graph/matrix: map the ancestor node to all its children. The weight of the edges (values in the matrix) is determined by
        $$
        wi,j​= 2IC(i) / IC(i)+IC(j)
        $$
        where IC is the information content determined by the structure of the GO graph. Note that i is the current node and j is the child node. Transition is favoured to child nodes of higher informational content.
        - Conditional Random walk
        Perform the same random walk as for the co-occurance, and obtain the Rh matrix that determines the possibility that if a parent node was annotated, the child node will later be annotated. 
        
3. Compute the negative GO terme:
        - Combine Rh and Wh into a single sparce matrix
        $$
        R = βR_h+(1-β)R_c
        $$
        We want the negative GO terms, so we must take 1-R. This is dense and huge, so computed only in the end (for the specific GO terms required, the full matrix remains sparse). 
        - For each  protein id, we take the list of positive GO terms, find the potential postitve GO terms for each of the protein's associated GO terms (R matrix). Find the negative GO terms by L = 1-R. Subset those GO terms that have high informational content and take the ones with the top 50 scores. This is the list of NegGO terms per protein​
    

### Protein epresentation with ESM2
Protein sequences were embedded using ESM2-T6 (8M parameters):
	- Encoder frozen for stability and efficiency
	- CLS token embeddings extracted from the final layer
For long sequences:
	- sequences were chunked (max length ~1000, stride ~970),
	- CLS embeddings were averaged across chunks to form a single protein representation.
All embeddings were precomputed and cached to disk, enabling fast experimentation without recomputing ESM features.

### Model Architecture
The prediction model is a lightweight MLP head:
ESM2 embedding (CLS)
 → Linear → ReLU → Dropout
 → Linear → ReLU → Dropout
 → Linear → logits over GO terms

### Loss Function Design
The final objective is a weighted combination of four losses, each addressing a specific failure mode.

1. Asymmetric BCE with IC Weighting:
False positives and false negatives are weighted per GO term:
	high-IC (specific) terms → stronger FP penalty,
	low-IC (general) terms → stronger FN penalty.
Losses are computed separately for BP, MF, and CC and averaged to prevent BP dominance.
We use the Information Content scores computed from the abundance (of GO terms across training proteins). We first nornalise it (to range from 0 to 1), then the weights for the false positive terms are 1 - ic_norm, while weights for false negative terms are 1 + ic_norm. While we are weighing each GO term based on tis IC, we must also note that we have 25k GO terms, at most 50 will be annotated per protein. This will mean that our loss will be min even when all go are predicted negative. We need to force some positives to be predicted, so we penalise false negatives harder than false positives.

2. Coverage Loss (Early Training Stabilization)
To prevent collapse to all-zero predictions, I introduced a coverage loss: enforces a minimum average predicted probability over true positive GO terms, applied in early epochs.
This ensures the model learns to “activate” meaningful outputs before fine-grained discrimination.

3. Hierarchy Consistency Loss
Enforces: P(child)≤P(parent). To keep computation tractable, I only applied to top-K predicted GO terms per batch, hierarchy pairs are precomputed once from the GO DAG.

4. Negative GO Margin Loss
Explicitly penalizes confident predictions on precomputed negative GO terms. This forces separation between true positives and hard negatives. This was activated only after a few epochs of training.
 
Final Loss
$$
L_{asym} + \lambda_{hier} L_{hier} + \lambda_{neg} L_{neg} + \lambda_{cov} L_{cov}
$$

### Performance
Given the computational power of my laptop, I could run the code only for 5 epochs (took 24 hrs), during which I observed a big decrease in training loss, but only a mild decrease in validation loss. The same was trend observed in the F1 score and AUROC scores (micro). 

The test set showed F1 micro score : XXXX and AUROC micro score: XXXX.

There are many other models that are computationally advance that could yield better results. 
1. Extract embeddings from multiple protein language models and input into an dense neural network, as done in [PROTGOAT](https://www.biorxiv.org/content/10.1101/2024.04.01.587572v1.full)
2. Use embeddings for GO terms (from the GO term description) rather than the one hot encoded like system here. Although I did account for relationships between the GO terms, using embeddings might increase accuracy. 
3. Incorporate more protein information, rather than just its sequence, like taxonomy (provided) and structural information (eg. propensity to disorder, alpha helix)

### My Reflections
My goal was to learn how to construct data structures and loss functions for a complex classification problem like this: I didn't have the means to test the above mentioned possibilites. In this process, alongside understanding the computer science concepts of deep learning, I also learnt to design a model that is RAM efficient and structure aware. 
1.	Dataset construction was the cornerstone of ensuring RAM efficiency: this enabled me to take in sequences between 30-35k amino acids long, and map them to 25k GO-terms, while ensuring limits such as the negative GO-terms and hierarchical relationships. 
2.	Incorporate as much prior knowledge as possible into the loss functions: my first model didn’t account for hierarchical IC increase, co-occurrence probabilities and non-overlapping nature of sub-ontologies, as I wanted the model to learn them implicitly through training. But by defining the loss function to incorporate these limits, I saw a much better performance (however, there could be other cases where it could lead to bias). 
3.	I also saw that scaling down can lead to massive improvements in time-memory efficiency, while retaining model accuracy: taking only top hits, defining cut-offs and making the data sparce were the difference between models that ran vs crashed. 
