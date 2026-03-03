---
layout: page
title: CAFA6 Challenge
description: Predicting a protein's function from its amino acid sequence
img: assets/img/12.jpg
importance: 1
category: work
related_publications: true
---

Can we predict a protein's function from its amino acid sequence? The CAFA (Critical Assessment of protein Function Annotation) 6th challenge aims to predict protein function prediction using Gene Ontology (GO) annotations. GO annotations are the function annotations associated with the protein, the GO terms are structured as a Directed Acyclic Graph (DAG), and belong to three non-overlapping ontologies: Molecular Function (MF), Biological Process (BP) and Cellular Component (CC). Each protein can have multiple GO terms, and annotations must respect the hierarchical structure of the ontology (true path rule).

This task is a multi-label classification problem:
    - Input: Protein amino acid sequence
    - Output: A probability score for each GO term

Challenges:
    - Sparce labelling / Class imbalance: most GO terms are rare
    - Hierarchical Label Structure: the labels are associated with each other. If a protein is associated with a GO term, it must also be annotated with that term's ancestor terms.
    - Varying protein sequence length: We must find a reliablw way to represent proteins that are 30 to 30,000 amino acids long
    - Large dataset, and limited computational power on my laptop
    
GO Term Standardization and Indexing: 
    - All GO terms from the CAFA training set were parsed using go-basic.obo. 
    - Only GO terms appearing in the training annotations were retained, yielding ~25k terms.
    - The GO terms were split into the 3 subontologies to process separately


<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/1.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/3.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/5.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    Caption photos easily. On the left, a road goes through a tunnel. Middle, leaves artistically fall in a hipster photoshoot. Right, in another hipster photoshoot, a lumberjack grasps a handful of pine needles.
</div>
<div class="row">
    <div class="col-sm mt-3 mt-md-0">
        {% include figure.liquid loading="eager" path="assets/img/5.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    This image can also have a caption. It's like magic.
</div>

You can also put regular text between your rows of images, even citations {% cite einstein1950meaning %}.
Say you wanted to write a bit about your project before you posted the rest of the images.
You describe how you toiled, sweated, _bled_ for your project, and then... you reveal its glory in the next row of images.

<div class="row justify-content-sm-center">
    <div class="col-sm-8 mt-3 mt-md-0">
        {% include figure.liquid path="assets/img/6.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm-4 mt-3 mt-md-0">
        {% include figure.liquid path="assets/img/11.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
    </div>
</div>
<div class="caption">
    You can also have artistically styled 2/3 + 1/3 images, like these.
</div>

The code is simple.
Just wrap your images with `<div class="col-sm">` and place them inside `<div class="row">` (read more about the <a href="https://getbootstrap.com/docs/4.4/layout/grid/">Bootstrap Grid</a> system).
To make images responsive, add `img-fluid` class to each; for rounded corners and shadows use `rounded` and `z-depth-1` classes.
Here's the code for the last row of images above:

{% raw %}

```html
<div class="row justify-content-sm-center">
  <div class="col-sm-8 mt-3 mt-md-0">
    {% include figure.liquid path="assets/img/6.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
  </div>
  <div class="col-sm-4 mt-3 mt-md-0">
    {% include figure.liquid path="assets/img/11.jpg" title="example image" class="img-fluid rounded z-depth-1" %}
  </div>
</div>
```

{% endraw %}
