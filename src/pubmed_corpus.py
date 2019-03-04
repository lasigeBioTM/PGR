import os
import subprocess
import sys
from polyglot.detect import Detector

import relations


#### ORIGINAL ABSTRACTS ####

def get_pubmed_ids_list(file_g2p, file_p2g):
    """Creates a list of lists of type [[pubmed_ID1, pubmed_ID2], ...]

    :param file_g2p: file with relations gene to phenotype
    :param file_p2g: file with relations phenotype to gene
    :return: list of lists with the ids of PubMed articles with genes in phenotype-gene
             relations of type [[pubmed_ID1, pubmed_ID2], ...]
    """

    pubmed_list = []

    dict_pg = relations.join_dicts(file_g2p, file_p2g)

    # QUERY ONE

    # for gene_name  in dict_pg:
    #
    #     os.system('curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=' + gene_name + '+homo+sapiens&retmax=50&retmode=xml" > articles.xml')
    #
    #     exit_file = open('articles.xml', 'r', encoding = 'utf-8')
    #     list_ids = exit_file.read().split('<IdList>')[-1].split('</IdList>')[0].split('\n')[1:-1]
    #     list_ids = [x.strip('</Id>') for x in list_ids]
    #
    #     pubmed_list.append(list_ids)
    #
    #     exit_file.close()

    # QUERY TWO

    for gene_name  in dict_pg:

        os.system('curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=' + gene_name + '+homo+sapiens+disease&retmax=50&retmode=xml" > articles.xml')

        exit_file = open('articles.xml', 'r', encoding = 'utf-8')
        list_ids = exit_file.read().split('<IdList>')[-1].split('</IdList>')[0].split('\n')[1:-1]
        list_ids = [x.strip('</Id>') for x in list_ids]

        pubmed_list.append(list_ids)

        exit_file.close()


    os.system('rm articles.xml')

    return pubmed_list


def write_text(file_g2p, file_p2g, number_abstracts_per_gene, destination_path):
    """Creates a file for each retrieved abstract

    :param file_g2p: file with relations gene to phenotype
    :param file_p2g: file with relations phenotype to gene
    :param number_abstracts_per_gene: int that indicates the number of abstracts intended for each gene
    :param destination_path: destination path
    :return: file for each retrieved abstract
    """

    pubmed_list = get_pubmed_ids_list(file_g2p, file_p2g)

    for list_ids in pubmed_list:

        number_requests = 0
        counter = 0

        while number_requests < number_abstracts_per_gene and counter < len(list_ids):

            os.system('curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=' + list_ids[counter] + '&retmode=xml" > abstract.xml')

            presence = subprocess.Popen("grep '<AbstractText>' 'abstract.xml'", shell = True)

            return_code = presence.wait()

            if return_code != 1:

                try:

                    exit_file = open('abstract.xml', 'r', encoding = 'utf-8')
                    abstract = exit_file.read().split('<AbstractText>', 1)[-1].split('</AbstractText>', 1)[0]
                    abstract = ''.join(x for x in abstract if x in abstract.printable)

                    save_language = ''

                    for language in Detector(abstract).languages:
                        save_language = str(language).split()[1]
                        break

                    if save_language == 'English':
                        output = open(destination_path + '/' + list_ids[counter], 'w', encoding = 'utf-8')
                        output.write(abstract)
                        output.close()

                        number_requests += 1

                    else:
                        print(list_ids[counter], 'was discarded:', save_language)

                    exit_file.close()

                except UnicodeDecodeError:

                    pass

                except FileNotFoundError:

                    pass

            counter += 1

    os.system('rm abstract.xml')

    return


#### DIVIDED BY SENTENCES ABSTRACTS ####

def divided_by_sentences(corpus_path, geniass_path, destination_path):  # needs to be run from directory bin/geniass/
    """Creates a file for each divided by sentences abstract

    :param corpus_path: edited corpus path
    :param geniass_path: GENIA Sentence Splitter path
    :param destination_path: destination path
    :return: file for each divided by sentences abstract
    """

    os.system('rm -rf ' + destination_path + '* || true')

    for (dir_path, dir_names, file_names) in os.walk(corpus_path):

        for filename in file_names:

            os.system('./' + geniass_path + ' ' + corpus_path + filename + ' ' + destination_path + filename)

    return


# #### OBSOLETE (FIRST TESTS) - EDITED ABSTRACTS ####
#
# def clean_abstract(corpus_path, destination_path):
#     """Creates a file for each "cleaned" abstract
#
#     :param corpus_path: original corpus path
#     :param destination_path: destination path
#     :return: file for each "cleaned" abstract
#     """
#
#     os.system('rm -rf ' + destination_path + '* || true')
#
#     for (dir_path, dir_names, file_names) in os.walk(corpus_path):
#
#         for filename in file_names:
#
#             misfit_list = ['29171490', '29472449', '29277180', '29868510', '29311298', '28687356',
#                            '29589160', '29518827', '29996362', '29578123', '29301908', '28062395',
#                            '30076191', '29078376', '29786049', '28673981', '27591164', '28235828',
#                            '29171469', '28183797', '29286619', '29431644', '29511098', '29336266',
#                            '29693365', '23897157', '28593008', '29171615']  # different languages, new lines in the middle of text, empty files and bad xml format
#
#             original = open(corpus_path + filename, 'r', encoding = 'utf-8')
#             original_contents = original.read()
#             original.close()
#
#             if filename not in misfit_list:
#
#                 edited = open(destination_path + filename, 'w', encoding = 'utf-8')
#
#                 # Important \xb1 is a +|- character; \u2122 is a TM character; \u2206 is a triangle (increment);
#                 # \xb7 is a middle dot; \u2192 is a rightwards arrow; \u0394 is a big delta; \xa9 is a copyright sign;
#                 # \xae is a registered sign; \u03b4 is a small delta; \u2265 is greater-than or equal to;
#                 # \u223c is a tilde operator; \u226b is much greater-than; \u2248 almost equal to; \u03bb is lambda;
#                 # \u03c1 is small rho; \u2264 is less-than or equal to; \u2267 is greater-than over equal to;
#                 # \u2061 is FA, function application; \u03a8 is small psi; \u2211 is n-ary summations;
#                 # \u2266 is less-than over equal to; \u25b3 is a triangle; \ufb01 is fi, latin small ligature;
#                 # \u03c0 is pi; \u03c6 is phi, greek small letter phi
#
#                 edited.write(original_contents.replace('<i>', '').replace('</i>', '').replace('<sup>', '').replace('</sup>', '')\
#                              .replace('<b>', '').replace('</b>', '').replace('<sub>', '').replace('</sub>', '').replace('[i]', '')\
#                              .replace('[/i]', '').replace("'", '*').replace('\u03b1', 'a').replace('\xa0', ' ')\
#                              .replace('\u03bc', 'u').replace('\xb0', 'o').replace('\u03b2', 'b').replace('\u202f', ' ')\
#                              .replace('\u2009', ' ').replace('\u03ba', 'k').replace('\xd7', 'x').replace('\u3000', ' ')\
#                              .replace('\xb1', '+').replace('\u2122', 'T').replace('\u2206', 't').replace('\u2217', '*')\
#                              .replace('\xe7', 'c').replace('\xb7', '-').replace('\u2011', '-').replace('\u03b5', 'e')\
#                              .replace('\u200a', ' ').replace('\xb5', 'u').replace('\u2192', '-').replace('\u0394', 'D')\
#                              .replace('\u03b3', 'y').replace('\xa9', 'c').replace('\u2113', 'l').replace('\xae', 'R')\
#                              .replace('\u03b4', 'd').replace('\u207a', '+').replace('\u2265', 'G').replace('\u223c', 'a')\
#                              .replace('\xe9', 'e').replace('\u207b','-').replace('\xeb', 'e').replace('\u2082', '2')\
#                              .replace('\xef', 'i').replace('\u2161', '2').replace('\xfc', 'u').replace('\u226b', 'M')\
#                              .replace('\u2248', 'a').replace('\u03bb', 'y').replace('\uff0c', ',').replace('\u2022', '.')\
#                              .replace('\u2005', ' ').replace('\xc5', 'A').replace('\xe8', 'e').replace('\u03c1', 'r')\
#                              .replace('\xf6', 'o').replace('\u0105', 'a').replace('\u2264', 'L').replace('\xe4', 'a')\
#                              .replace('\u03c7', 'x').replace('\u2012', '-').replace('\xed', 'i').replace('\u2267', 'G')\
#                              .replace('\xdf', 'b').replace('\xe0', 'a').replace('\u2061','F').replace('\u03a8', 'p')\
#                              .replace('\u3001', ',').replace('\u03c9', 'w').replace('\u0131', 'i').replace('\u015f', 's')\
#                              .replace('\u011f', 'g').replace('\xf3', 'o').replace('\xe1', 'a').replace('\u2013', '-')\
#                              .replace('\u0440', 'p').replace('\u2003', ' ').replace('\u2211', 'S').replace('\u039c', 'M')\
#                              .replace('\u2029', ' ').replace('\xf4', 'o').replace('\xb2', '2').replace('\u2160', 'I')\
#                              .replace('\u2162', '3').replace('\u2081', '1').replace('\u0430', 'a').replace('\u03c3', 'o')\
#                              .replace('\u03b6', 'z').replace('\u03b8', 'o').replace('\u2163', '4').replace('\u2236', ':')\
#                              .replace('\u201c', '"').replace('\u201d', '"').replace('\u2266', 'L').replace('\u03a7', 'X')\
#                              .replace('\u2019', '"').replace('\u2076', '6').replace('\xb9', '1').replace('\u25b3', 'T')\
#                              .replace('\xec', 'i').replace('\u03b7', 'n').replace('\ufb01', 'f').replace('\u0125', 'h')\
#                              .replace('\u03c0', 'p').replace('\u025b', 'e').replace('\u03c6', 'f').replace('\u02c2', '<'))
#
#                 edited.close()
#
#     return


# #### OBSOLETE - A SENTENCE PER FILE (ORIGINAL SENTENCES) ####
#
# def original_sentences(corpus_path, destination_path):
#     '''Creates a file for each sentence in each abstract
#
#     :param corpus_path: divided by sentences corpus path
#     :param destination_path: destination path
#     :return: file for each sentence in each abstract
#     '''
#
#     for (dir_path, dir_names, file_names) in os.walk(corpus_path):
#
#         for filename in file_names:
#
#             abstract_file = open(corpus_path + filename, 'r', encoding = 'utf-8')
#             abstract = abstract_file.readlines()
#             abstract_file.close()
#
#             sentence_id = 0
#
#             for line in abstract:
#
#                 sentences_file = open(destination_path + filename + 's' + str(sentence_id), 'w', encoding = 'utf-8')
#                 sentences_file.write(line)
#                 sentences_file.close()
#
#                 sentence_id += 1
#
#     return


#### RUN ####

def main():
    """Creates a directory with a file for each retrieved abstract divided by sentences

    :return: directory with a file for each retrieved abstract divided by sentences
    """

    number_of_abstracts_per_gene = int(sys.argv[1])

    os.system('mkdir -p corpora/pubmed_corpus/ || true')
    write_text('data/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt',
               'data/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt', number_of_abstracts_per_gene, 'corpora/pubmed_corpus/')
    os.system('mkdir -p corpora/edited_corpus/ || true')
    os.system('cp corpora/pubmed_corpus/* corpora/edited_corpus/')
    os.chdir('bin/geniass/')
    divided_by_sentences('../../corpora/edited_corpus/', 'geniass', '../../corpora/pubmed_corpus/')
    os.chdir('../..')
    os.system('rm -rf corpora/edited_corpus/')

    return


if __name__ == "__main__":
    main()
