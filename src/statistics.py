import os
import datetime


#### GENERAL STATISTICS FOR THE CORPUS CREATION ####

def general_statistics(abstracts_path, annotations_path):
    """Reports the final numbers for each processing step, from retrieving the abstracts to the relations extraction

    :param abstracts_path: abstracts path
    :param annotations_path: annotations path
    :return: file .txt that reports the count of the abstracts, their annotations
             (phenotype and gene) and the relations between the annotations
    """

    report = open('report.txt', 'w', encoding = 'utf-8')

    current_time = datetime.datetime.now()

    report.write('--------------------------- STATISTICS REPORT (' + str(current_time) + ') ---------------------------')
    report.write('\n\n\n')

    report.write('------------------- ABSTRACTS -------------------')
    report.write('\n\n')

    final_number_abstracts = os.popen('ls -l ' + abstracts_path + '* | egrep -c \'^-\'').read()

    report.write('FINAL NUMBER OF ABSTRACTS ------> ' + str(final_number_abstracts))

    report.write('\n\n')
    report.write('------------------ ANNOTATIONS ------------------')
    report.write('\n\n')

    total_final_annotations = 0
    total_final_phenotype_annotations = 0
    total_final_gene_annotations = 0

    for (dir_path, dir_names, file_names) in os.walk(annotations_path + 'final_gene_annotations/'):

        for filename in file_names:

            annotation_file = open(annotations_path + 'final_gene_phenotype_annotations/' + filename, 'r', encoding = 'utf-8')
            contents = annotation_file.readlines()
            annotation_file.close()

            total_final_annotations += len(contents)

            for annotation in contents:

                if 'HP_' in annotation:

                    total_final_phenotype_annotations += 1

                else:

                    total_final_gene_annotations += 1

    report.write('FINAL TOTAL NUMBER OF ANNOTATIONS ---------> ' + str(total_final_annotations))
    report.write('\n')
    report.write('FINAL NUMBER OF PHENOTYPE ANNOTATIONS -----> ' + str(total_final_phenotype_annotations))
    report.write('\n')
    report.write('FINAL NUMBER OF GENE ANNOTATIONS ----------> ' + str(total_final_gene_annotations))
    report.write('\n')

    report.write('\n\n')
    report.write('------------------- RELATIONS -------------------')
    report.write('\n\n')

    relations = open(annotations_path + 'relations.tsv', 'r', encoding = 'utf-8')
    relations_contents = relations.readlines()
    relations.close()

    total_relations = len(relations_contents) - 1
    total_true_relations = 0
    total_false_relations = 0

    for relation in relations_contents:

        if 'True' == relation.split('\t')[-1][:-1]:

            total_true_relations += 1

        elif 'False' == relation.split('\t')[-1][:-1]:

            total_false_relations += 1

    report.write('FINAL TOTAL NUMBER OF RELATIONS -------> ' + str(total_relations))
    report.write('\n')
    report.write('FINAL NUMBER OF TRUE RELATIONS --------> ' + str(total_true_relations))
    report.write('\n')
    report.write('FINAL NUMBER OF FALSE RELATIONS -------> ' + str(total_false_relations))
    report.write('\n')

    return

#### RUN ####

general_statistics('corpora/pubmed_corpus/', 'corpora/')
