# PGR: A Silver Standard Corpus of Human Phenotype-Gene Relations

The PGR corpus is a silver standard corpus of human phenotype and gene annotations and their relations. This corpus is available in the **corpora/10_12_2018_corpus/** directory (in *.tsv* and *.xml* formats). Later, a new corpus was created using a different query, available at the **corpora/11_03_2019_corpus/** directory (in *.tsv* and *.xml* formats).
If you intend to create a new corpus you can follow the bellow guidelines.

Our academic paper which describes PGR in detail can be found [here](https://aclweb.org/anthology/papers/N/N19/N19-1152/).

## Dependencies

* Python >= 3.5

* Pre-processing:
    * [Genia Sentence Splitter](http://www.nactem.ac.uk/y-matsu/geniass/)
    
* Term Recognition:
    * [MER (Minimal Named-Entity Recognizer)](https://github.com/lasigeBioTM/MER) (Gene Entities)
    * [IHP (Identifying Human Phenotype Entities)](https://github.com/lasigeBioTM/IHP) (Human Phenotype Entities)
        
* Relation Extraction:
    * [Human Phenotype Ontology Gold Standard Relations](https://hpo.jax.org/app/download/annotation) (Knowledge Base)
    * [Gene2Go Correspondence File](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz) (To facilitate the use of the [BO-LSTM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2584-5#Abs1) application.) [OPTIONAL]
    
## Getting Started

````
 cd bin/
 git clone git@github.com:lasigeBioTM/MER.git
 git clone -b IHP_Python3.6 --single-branch git@github.com:lasigeBioTM/IHP.git
````
Use the Dockerfile to setup the rest of the experimental environment or the [PGR Image](https://hub.docker.com/r/dpavot/pgr) available at Docker Hub.
   
## Usage

Run Stanford CoreNLP for the IHP to be able to annotate the human phenotype entities.

````
 cd bin/IHP/bin/stanford-corenlp-full-2015-12-09/
 java -mx4g -cp "*" edu.stanford.nlp.pipeline.StanfordCoreNLPServer -timeout 500000 &
````

### Retrieving Abstracts

````
 python3 src/pubmed_corpus.py [NUMBER]
````

where [NUMBER] (integer) corresponds to the intended number of abstracts per gene that participates in human phenotype-gene relations.

* Creates: 
    * **corpora/pubmed_corpus/**

### Annotating Genes, Human Phenotypes and Relations

````
 python3 src/annotations.py
````

* Creates: 
    * **corpora/gene_phenotype_annotations/** 
    * __corpora/relations.tsv__
    
* Changes:
    * **corpora/pubmed_corpus/** (removes abstracts that do not have entities from both types)

### Creating a XML Format Corpus

````
 python3 src/pgr_corpus.py [ENTITY TYPE]
````

where [ENTITY TYPE] (*gene* or *go*) corresponds to the intended pair of entities (human phenotype-gene pair or human phenotype-go pair) to generate an XML format corpus with. The GO (Gene Ontology) term corresponds to the most representative term for the gene that establishes the relation with that human phenotype.

* Creates: 
    * **corpora/pgr_gene/** (with [ENTITY TYPE] = *gene*)
    * **corpora/go_phenotype_annotations/** (with [ENTITY TYPE] = *go*)
    * **corpora/pgr_go/** (with [ENTITY TYPE] = *go*)
    
### General Statistics

````
 python3 src/statistics.py
````

* Creates: 
    * __report.txt__

## Configuration

* ### bin/
    * **MER/**
        * **data/**
            * __genes.txt__
            * __genes_links.tsv__
    * **IHP/**
    * **geniass/**
    
* ### corpora/
    * **10_12_2018_corpus/**
        * **pgr_test/**
            * **pgr_gene/**
            * **pgr_go/**
        * **pgr_train/**
            * **pgr_gene/**
            * **pgr_go/**
        * __test.tsv__
        * __train.tsv__
    * **11_03_2019_corpus/**
        * **pgr_train/**
            * **pgr_gene/**
            * **pgr_go/**
        * __train.tsv__    
          

* ### data/
    * __ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt__
    * __ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt__
    * __gene2go.gz__
    
* ### src/
    * **annotations.py**
    * **pgr_corpus.py**
    * **pubmed_corpus.py**
    * **relations.py**
    * **statistics.py**
    
## Reference

- Diana Sousa, Andre Lamurias, and Francisco M. Couto. 2019. A Silver Standard Corpus of Human Phenotype-Gene Relations. In Proceedings of the 2019 Conference of the North American Chapter of the Association for Computational Linguistics: Human Language Technologies, Volume 1 (Long and Short Papers), pages 1487–1492.
