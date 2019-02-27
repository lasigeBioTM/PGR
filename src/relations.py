import collections


def dict_g2p(file_g2p, dict_type = None):
    """Creates a dictionary of type {gene1:[phenotype1, phenotype2, ...], }

    :param file_g2p: file with relations gene to phenotype
    :param dict_type: dict_type (optional) to create a dict with the names
                      (default) or ids of the genes and respective phenotypes
    :return: dict with the names or ids of type {gene1:[phenotype1, phenotype2, ...], }
    """

    gene2phenotype = open(file_g2p, 'r')

    gene2phenotype.readline()

    relations_g2p = gene2phenotype.readlines()
    
    relations_g2p.pop()

    gene2phenotype.close()

    if not dict_type:

        dict_gene_phenotype = {}

        for line in relations_g2p:

            line = line.split('\t')

            gene_name = line[1].lower()
            phenotype_name = line[2].lower()

            if gene_name not in dict_gene_phenotype:
                dict_gene_phenotype[gene_name] = []
                dict_gene_phenotype[gene_name].append(phenotype_name)

            else:
                dict_gene_phenotype[gene_name].append(phenotype_name)

        return dict_gene_phenotype

    else:
        dict_geneID_phenotypeID = {}

        for line in relations_g2p:

            line = line.split('\t')

            geneID = line[0]
            phenotypeID = line[3][:-1].replace(':', '_')

            if geneID not in dict_geneID_phenotypeID:
                dict_geneID_phenotypeID[geneID] = []
                dict_geneID_phenotypeID[geneID].append(phenotypeID)

            else:
                dict_geneID_phenotypeID[geneID].append(phenotypeID)

        return dict_geneID_phenotypeID


def dict_p2g(file_p2g, dict_type = None):
    """Creates a dictionary of type {gene1:[phenotype1, phenotype2, ...], }

    :param file_p2g: file with relations phenotype to gene
    :param dict_type: dict_type (optional) to create a dict with the names
                      (default) or ids of the genes and respective phenotypes
    :return: dict with the names or ids of type {gene1:[phenotype1, phenotype2, ...], }
    """

    phenotype2gene = open(file_p2g, 'r')

    phenotype2gene.readline()

    relations_p2g = phenotype2gene.readlines()

    relations_p2g.pop()

    phenotype2gene.close()

    if not dict_type:

        dict_gene_phenotype = {}

        for line in relations_p2g:

            line = line.split('\t')

            phenotype_name = line[1].lower()
            gene_name = line[3][:-1].lower()

            if gene_name not in dict_gene_phenotype:
                dict_gene_phenotype[gene_name] = []
                dict_gene_phenotype[gene_name].append(phenotype_name)
        
            else:
                
                dict_gene_phenotype[gene_name].append(phenotype_name)

        return dict_gene_phenotype

    else:
        
        dict_geneID_phenotypeID = {}

        for line in relations_p2g:

            line = line.split('\t')

            phenotypeID = line[0].replace(':', '_')
            geneID = line[2]

            if geneID not in dict_geneID_phenotypeID:
                dict_geneID_phenotypeID[geneID] = []
                dict_geneID_phenotypeID[geneID].append(phenotypeID)

            else:
                dict_geneID_phenotypeID[geneID].append(phenotypeID)

        return dict_geneID_phenotypeID


def join_dicts(file_g2p, file_p2g, dict_type = None):
    """Creates a joined merged dictionary of type {gene1:[phenotype1, phenotype2, ...], }

    :param file_g2p: file with relations gene to phenotype
    :param file_p2g: file with relations phenotype to gene
    :param dict_type: dict_type (optional) to create a dict with the names
                      (default) or ids of the genes and respective phenotypes
    :return: dict with the names or ids of type {gene1:[phenotype1, phenotype2, ...], }
    """

    if not dict_type:

        dict_phenotype_gene = dict_p2g(file_p2g)
        dict_gene_phenotype = dict_g2p(file_g2p)
        
        dict_relations = collections.defaultdict(list, dict_phenotype_gene)

        for key, value in dict_gene_phenotype.items():

            dict_relations[key].extend(value)

        for key, value in dict_relations.items():  # unique values only

            dict_relations[key] = list(set(value))

        return dict(dict_relations)
    
    else:

        dict_phenotypeID_geneID = dict_p2g(file_p2g, 1)
        dict_geneID_phenotypeID = dict_g2p(file_g2p, 1)

        dict_relationsID = collections.defaultdict(list, dict_phenotypeID_geneID)

        for key, value in dict_geneID_phenotypeID.items():

            dict_relationsID[key].extend(value)

        for key, value in dict_relationsID.items():  # unique values only

            dict_relationsID[key] = list(set(value))

        return dict(dict_relationsID)
