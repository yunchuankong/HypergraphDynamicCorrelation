# HypergraphDynamicCorrelation

Relevant code for the project Hypergraph-based Dynamic Correlation Detection in High-throughput Data.

## Prerequisites

The following packages are required for executing the main code file: igraph, coin, ks, dynamicTreeCut, GOstats. For the sample dataset in the below example, org.Sc.sgd.db is also needed.

## Usage

### Data formats

* A gene expression matrix (called "array"), seen in spellman_73_filled.bin for this example
* A list of GO gene module selection (called "GO.select"), seen in GO_select_yeast_ 0.7 0.4 1000 20 .bin for this example

### Preliminary process

Suposse all relevant files are in the same folder, and the work directory has been set properly.

```
source("main.r")
load('spellman_73_filled.bin')
load("GO_select_yeast_ 0.7 0.4 1000 20 .bin")
array = array[sample(nrow(array), 1000),]
```

### Unsupervised hypergraph construction

To construct module-level hyergraph, use
```
hp_us <- hypergraph_unsup(array, fdr=0.2, save_glist=T)
```
and visualize the full hypergraph and top-connected sub-hypergraph by
```
full_hypergraph_us <- plot_entrie_unsup(hp_us$net, hp_us$label, folds=2)
plot_top_unsup(full_hypergraph_us$g, full_hypergraph_us$elist, hp_us$label, full_hypergraph_us$folds)
```
.
