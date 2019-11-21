import os
import re
import collections
import random

import relations


#### GENES PREPARATION FOR MER ####

def get_genes_ids(mer_data_path):
    """Creates a lexicon file .tsv with a gene name and gene ID per line

    :param mer_data_path: source and destination path
    :return: lexicon file .tsv with a gene name and gene ID per line
    """

    lexicon = open(mer_data_path + 'genes.txt', 'r', encoding = 'utf-8')  # https://github.com/macarthur-lab/gene_lists/blob/master/lists/universe.tsv
    lexicon_names = lexicon.readlines()
    lexicon.close()

    lexicon_ids = open(mer_data_path + 'genes_links.tsv', 'w', encoding = 'utf-8')  # links are NCBI identifiers

    empty_list = []

    for gene_name in lexicon_names:

        gene_id = os.popen('curl http://mygene.info/v3/query?q=' +\
                           gene_name[:-1] + '&species=human').readlines()  # API that matches gene names to gene NCBI IDs

        if gene_id[1][-6:-2] == 'null':

            empty_list.append(gene_name[:-1])

        else:

            index = 0

            for line in gene_id:

                if 'entrezgene' in line:
                    gene_id = line.split('"entrezgene": "')[-1][:-3]

                    index = 1

                    break

            if index == 0:

                empty_list.append(gene_name[:-1])

            else:

                lexicon_ids.write(gene_name[:-1] + '\t' + gene_id + '\n')

    lexicon_ids.close()

    return empty_list


#### GENES DICTIONARY: GENE ANNOTATION - ID ####

def dict_genes(mer_genes_links_file):
    """Creates a dictionary of type {gene_name:gene_id, }

    :param mer_genes_links_file: .tsv file with one gene name and respective gene ID per line
    :return: dict of type {gene_name:gene_id, }
    """

    lexicon = open(mer_genes_links_file, 'r', encoding = 'utf-8')
    lexicon_names_id = lexicon.readlines()
    lexicon.close()

    dict_gene_name_id = {}

    for gene_name_id in lexicon_names_id:

        dict_gene_name_id[gene_name_id.split('\t')[0]] = gene_name_id.split('\t')[1][:-1]

    return dict_gene_name_id


#### GENES AND GENES LINKS FILE UPDATED WITH SYNONYMS ####

def get_genes_synonyms(data_path, mer_data_path):
    """Updates genes.txt and genes_links.tsv files with the synonyms of the genes + their ID number

    :param data_path: path to data
    :param mer_data_path: path to mer data
    :return: updated genes.txt and genes_links.tsv files with the synonyms of the genes + their ID number
    """

    synonyms_file = open(data_path + 'gene_symbol_thesaurus.txt', 'r', encoding = 'utf-8')  # https://github.com/macarthur-lab/gene_lists/blob/master/other_data/gene_symbol_thesaurus.txt
    synonyms_file.readline()  # skip header
    synonyms = synonyms_file.readlines()

    try:

        dict_gene_name_id = dict_genes(data_path + 'no_synonyms_genes_links.tsv')  # GENE ANNOTATION - ID

    except FileNotFoundError:

        dict_gene_name_id = dict_genes(mer_data_path + 'genes_links.tsv')  # GENE ANNOTATION - ID
        os.system('mv ' + mer_data_path + 'genes.txt ' + data_path + 'no_synonyms_genes.txt || true')
        os.system('mv ' + mer_data_path + 'genes_links.tsv ' + data_path + 'no_synonyms_genes_links.tsv || true')

    synonym_id_dict = {}

    for line in synonyms:

        synonym = line.split('\t')[1][:-1]

        if synonym == 'NA':
            pass

        else:

            synonyms_list = synonym.split(',')

            for element in synonyms_list:

                # Skips synonyms of genes not present in the original list of genes (because this list is older)
                # and synonyms with only one character (ratio TP/FP doesn't pay off)
                if line.split('\t')[0] in dict_gene_name_id and element not in synonym_id_dict and len(element) > 1:

                    synonym_id_dict[element] = dict_gene_name_id[line.split('\t')[0]]

    dict_join_all_genes = {**dict_gene_name_id, **synonym_id_dict}  # only works for python versions > 3.5
    dict_join_all_genes = dict(sorted(dict_join_all_genes.items(), key=lambda t : tuple(t[0].lower())))  # sort alphabetically

    new_lexicon_file = open(mer_data_path + 'genes.txt', 'w', encoding = 'utf-8')
    new_lexicon_id_file = open(mer_data_path + 'genes_links.tsv', 'w', encoding = 'utf-8')

    for gene_name, id in dict_join_all_genes.items():
        new_lexicon_file.write(gene_name + '\n')
        new_lexicon_id_file.write(gene_name.lower() + '\t' + id + '\n')

    new_lexicon_file.close()
    new_lexicon_id_file.close()

    return


#### HP DICTIONARY: HP ANNOTATION - ID ####

def dict_hp(mer_hp_links_file):
    """Creates a dictionary of type {hp_name:hp_id, }

    :param mer_hp_links_file: .tsv file with one hp name and respective hp ID per line
    :return: dict of type {hp_name:hp_id, }
    """

    hp_links_file = open(mer_hp_links_file, 'r', encoding = 'utf-8')
    hp_links = hp_links_file.readlines()
    hp_links_file.close()
    print(hp_links)
    dict_hp_name_id = {}

    for line in hp_links:

        line = line.split('\t')

        hp_name = line[0]
        hp_id = line[-1].split('/obo/')[-1][:-1]

        dict_hp_name_id[hp_name] = hp_id

    return dict_hp_name_id


#### IHP NEEDS FEWER FILES TO BE ABLE TO ANNOTATE HP ####

def divide_directory(corpus_path, directory_path):
    """Divides files from the corpus directory into seven different directories

    :param corpus_path: corpus path
    :param directory_path: per directory corpus path
    :return: directories with divided by sentences corpus (each with 500 abstracts and one with the leftovers)
    """

    for (dir_path, dir_names, file_names) in os.walk(corpus_path):

        number_directories = int(len(file_names) / 500) + 1

        for u in range(number_directories):
            os.system('mkdir -p ' + directory_path + 'directory_' + str(u) + ' || true')

    for (dir_path, dir_names, file_names) in os.walk(corpus_path):

        for i in range(number_directories):

            counter = 0
            os.system('rm -rf ' + directory_path + '/directory_' + str(i) + '/* || true')

            for filename in file_names:

                counter += 1

                if counter <= 500:
                    os.system('cp ' + corpus_path + filename + ' ' + directory_path + '/directory_' + str(i) + '/')

                else:
                    file_names = file_names[500:]
                    break

    return number_directories


#### JOIN REPORT IHP FILES IN ONE FINAL REPORT ####

def join_report_files(annotations_path):
    """Joins the seven report files generated by the IHP tool (see docstring from divided_directory) in one final report file

    :param annotations_path: per directory annotations path
    :return: file with joined report files
    """

    final_report = open(annotations_path + 'final_report', 'w', encoding = 'utf-8')

    for (dir_path, dir_names, file_names) in os.walk(annotations_path):

        for filename in file_names:

            if filename.startswith('directory'):

                directory_file = open(annotations_path + filename, 'r', encoding = 'utf-8')
                directory_content = directory_file.read().split('>')[-1][2:-1]
                directory_file.close()

                final_report.write(directory_content + '\n')

    final_report.close()

    return


#### HP AND GENES ANNOTATIONS ####

def annotations(corpus_path, mer_path, ihp_report_file, mer_hp_links_file, destination_path):  # needs to be run from directory bin/MER/
    """Creates an annotation file for each abstract in the corpus and an annotations to check file for bad (half/invalid) IHP annotations

    :param corpus_path: divided by sentences corpus path
    :param mer_path: mer data path
    :param ihp_report_file: IHP final report file
    :param mer_hp_links_file: MER HP links file
    :param destination_path: destination path
    :return: annotation file for each abstract in the corpus and a annotations to check file for bad (half/invalid) IHP annotations
             annotation file example:

            55  59  mTOR    2445
            396 401 IF2a    83939
            1917    1923    cancer  HP_0002664

    """

    dict_hp_name_id = dict_hp(mer_hp_links_file)

    for (dir_path, dir_names, file_names) in os.walk(corpus_path):

        check_by_hand = open(destination_path + 'annotations_to_check.tsv', 'w',
                             encoding = 'utf-8')  # half/invalid annotations to verify by hand using hp_links.tsv and
        # final_report (and add to respective annotation file if appropriate)

        for filename in file_names:

            annotations = []

            abstract_file = open(corpus_path + filename, 'r', encoding = 'utf-8')
            abstract = str(abstract_file.read().encode('utf-8'))
            abstract_file.close()

            annotations_genes = os.popen('./' + mer_path + ' ' + abstract + ' genes').readlines()  # genes annotations | abstract needs to be between quotes

            for annotations_gene in annotations_genes:

                gene_index_1 = annotations_gene.split('\t')[0]
                gene_index_2 = annotations_gene.split('\t')[1]
                gene_name = annotations_gene.split('\t')[2]
                gene_id =  annotations_gene.split('\t')[3][:-1]

                annotations.append((gene_index_1, gene_index_2, gene_name, gene_id))

            hp_all_annotations_file = open(ihp_report_file, 'r', encoding = 'utf-8')
            hp_all_annotations_lines = hp_all_annotations_file.readlines()
            hp_all_annotations_file.close()

            for line in hp_all_annotations_lines:

                if filename in line:

                    if line.startswith(filename):  # pass title
                        pass

                    else:
                        if line.split('\t')[-1].split(':')[-1].split('\t')[-1][:-1].lower() in dict_hp_name_id:

                            annotations.append((line.split('\t', 1)[-1].split(':')[0],
                                               line.split('\t', 1)[-1].split(':')[-1].split('\t')[0],
                                               line.split('\t', 1)[-1].split(':')[-1].split('\t')[-1][:-1].lower(),
                                               dict_hp_name_id[line.split('\t')[-1].split(':')[-1].split('\t')[-1][:-1].lower()]))

                        else:  # IHP annotations not present in hp_links.tsv dictionary

                            check_by_hand.write(filename + '\t' + line.split('\t')[-1].split(':')[-1].split('\t')[-1][:-1].lower() + '\n')

            annotations = sorted(annotations, key=lambda position: int(position[0]))  # sort by position in text

            annotation_file = open(destination_path + 'divided_by_sentences_annotations/' + filename, 'w', encoding = 'utf-8')

            for annotation in annotations:
                annotation_file.write(annotation[0] + '\t' + annotation[1] + '\t' + annotation[2] + '\t' + annotation[3] + '\n')

            annotation_file.close()

        check_by_hand.close()

    return


#### CHECK (HALF/INVALID) IHP ANNOTATIONS NOT PRESENT IN HP DICTIONARY: HP ANNOTATION - ID  ####

def check_not_in_dict_annotations(corpus_path, annotations_path):
    """Prints abstract + annotations to compare with annotations_to_check file and manually correct (half/invalid) IHP annotations

    :param corpus_path: divided by sentences corpus path
    :param annotations_path: divided by sentences annotations path
    :return: printed abstract + annotations to compare with annotations_to_check file
             and manually correct (half/invalid) IHP annotations
             usage example:

             Do you want to continue? (y/jump/n) y

             ------------------------------ 28222718 ------------------------------

             Frank-ter Haar syndrome pFTHSp is an autosomal recessive disorder characterized by abnormalities that
             affect the development of bone, heart, and eyes. We report a sibling pair with FTHS caused by a homozygous,
             novel mutation pLys133Glnfs*13 in the SH3PXD2B gene: one sibling had bilateral ocular hypertension and
             unilateral colobomas of iris, choroid and retina; the other, unilateral myelinated nerve fiber layer of
             the optic disk and papilledema due to idiopathic intracranial hypertension.
             Both children had refractive amblyopia and megalocornea.


             248	256	SH3PXD2B	285590
             279	288	bilateral	HP_0012832
             313	323	unilateral	HP_0012833
             436	447	papilledema	HP_0001085
             466	491	intracranial hypertension	HP_0002516

             Do you want to continue? (y/jump/n) jump

             Insert file number: 28223545

             ------------------------------ 28223545 ------------------------------

             (...)

             Do you want to continue? (y/jump/n) n

    """

    dict_files_order = {}

    for (dir_path, dir_names, file_names) in os.walk(annotations_path + 'divided_by_sentences_annotations/'):

        position = 0

        for filename in file_names:

            dict_files_order[position] = filename
            position += 1

        position_grabber = 0

        for iteration in range(len(file_names)):

            abstract = open(corpus_path + dict_files_order[position_grabber], 'r', encoding = 'utf-8')
            abstract_content = abstract.read()
            abstract.close()

            annotations = open(annotations_path + 'divided_by_sentences_annotations/' + dict_files_order[position_grabber], 'r', encoding = 'utf-8')
            annotations_content = annotations.read()
            annotations.close()

            print('------------------------------', dict_files_order[position_grabber], '------------------------------')
            print('\n')
            print(abstract_content)
            print('\n')
            print(annotations_content)

            next_file = input('Do you want to continue? (y/jump/n) ')  # y, jump or something else to stop the program
            print('\n')

            if next_file == 'y':

                position_grabber += 1

            elif next_file == 'jump':  # to jump to specific file

                desired_file = input('Insert file number: ')

                for position, filename in dict_files_order.items():

                    if filename == desired_file:
                        position_grabber = position

            else:
                break

    return


#### DICTIONARY WITH MANUALLY CORRECTED/JOINED ANNOTATIONS FROM ANNOTATIONS TO CHECK FILE ####

def manual_annotations(manual_annotations_file):
    """Creates a dictionary with manually corrected/joined annotations from annotations to check file of type {filename:[annotation1, annotation2, ...], }

    :param manual_annotations_file: .tsv divided by sentences annotations to check file
    :return: dict of type {filename:[annotation1, annotation2, ...], }
    """

    os.system('mv ' + manual_annotations_file + ' data/manual_annotations.tsv || true')  # save

    manual_annotations = open('data/manual_annotations.tsv', 'r', encoding = 'utf-8')
    manual_annotations_content = manual_annotations.readlines()
    manual_annotations.close()

    dict_manual_annotations = {}

    for manual_annotation in manual_annotations_content:

        if len(manual_annotations_content[0].split('\t')) > 2:  # check if manual annotations were added

            if manual_annotation.split('\t')[2] == 'Y':  # true annotations

                if manual_annotation.split('\t')[0] not in dict_manual_annotations:

                    dict_manual_annotations[manual_annotation.split('\t')[0]] = []
                    dict_manual_annotations[manual_annotation.split('\t')[0]].append(manual_annotation.split('Y')[-1][1:])

                else:

                    dict_manual_annotations[manual_annotation.split('\t')[0]].append(manual_annotation.split('Y')[-1][1:])

    return dict_manual_annotations


#### DICTIONARY WITH PREVIOUSLY UNDETECTED ANNOTATIONS ####

def caught_annotations(corpus_path, mer_data_path):
    """Creates a dictionary with previously undetected annotations from missed patterns (by MER and IHP)
       of type {filename:[annotation1, annotation2, ...], }

    :param corpus_path: divided by sentences corpus path
    :param mer_data_path: mer data path
    :return: dict of type {filename:[annotation1, annotation2, ...], }
    """

    dict_hp_name_id = dict_hp(mer_data_path + 'hp_links.tsv')  # HP ANNOTATION - ID

    dict_gene_name_id = dict_genes(mer_data_path + 'genes_links.tsv')  # GENE ANNOTATION - ID

    dict_caught_annotations = {}

    for (dir_path, dir_names, file_names) in os.walk(corpus_path):

        for filename in file_names:

            abstract = open(corpus_path + filename, 'r', encoding = 'utf-8')
            abstract_content = abstract.read()
            abstract.close()

            ### GENES
            for gene_name, gene_id in dict_gene_name_id.items():

                missed_patterns = ['-' + gene_name + '[ ,;]', '[ /]' + gene_name + '-']  # detected missed patterns by MER

                # The rejected list has acronyms or words that match genes if searched by their .lower() form
                # Also includes the genes on super reject list, in .lower() form, that are common acronyms or words
                rejected_list = ['-type[ ,;]', '[ /]up-', '-mb[ ,;]', '-ter[ ,;]', '-an[ ,;]', '[ /]co-',
                                 '[ /]light-', '-net[ ,;]', '-set[ ,;]', '[ /]spart-', '-all[ ,;]',
                                 '-rank[ ,;]', '-ray[ ,;]', '-arm[ ,;]', '[ /]pcd-', '-cal[ ,;]', '[ /]hip-',
                                 '[ /]mul-', '[ /]er-', '-wave[ ,;]', '-flash[ ,;]', '[ /]on-', '[ /]rod-',
                                 '-step[ ,;]', '[ /]fh-', '-ii[ ,;]', '-hip[ ,;]', '[ /]mis-', '-out[ ,;]',
                                 '[ /]ct-', '-stop[ ,;]', '[ /]fat-', '[ /]ph-', '-rod[ ,;]', '[ /]di-', '-ps[ ,;]',
                                 '-sex[ ,;]', '[ /]sex-', '-fat[ ,;]', '[ /]in-', '-end[ ,;]', '-flap[ ,;]',
                                 '-in[ ,;]', '[ /]ende-', '[ /]fast-', '-up[ ,;]', '[ /]end-', '-as[ ,;]',
                                 '[ /]set-', '-apr[ ,;]', '-jun[ ,;]', '[ /]mar-', '[ /]all-', '[ /]type-',
                                 '[ /]oct-', '[ /]cut-', '-lobe[ ,;]', '[ /]gap-', '-mass[ ,;]', '-spatial[ ,;]',
                                 '[ /]task-', '[ /]ii-']

                # -type | up- | -Mb | -ter | -an | co- | light- | -net | -set | SPART- | -ALL | -rank | -ray | -arm
                # PCD- | -cal | hip- | mul- | ER- | -wave | -flash | ON- | rod- | -step | FH- | -II, | -hip | mis- | -out
                # CT- | -stop | fat- | PH- | -rod | di- | -PS | -sex | sex- | -fat | in- | -end | -flap | -in | Ende- | fast-
                # -up | end- | -as | set- | -Apr; | -Jun; | Mar- | all- | type- | Oct- | cut- | -lobe | gap-
                # -mass | -spatial | task- | II-

                # The super rejected list has (capitalized) genes that are common acronyms or words
                super_rejected_list = ['[ /]SPART-', '-ALL[ ,;]', '[ /]PCD-', '[ /]ER-', '[ /]ON-', '[ /]FH-',
                                       '-II[ ,;]', '[ /]CT-', '[ /]PH-', '-PS[ ,;]', '[ /]Ende-', '[ /]II-']

                # SPART- | -ALL | PCD- | ER- | ON- | FH- | -II, | CT- | PH- | -PS | Ende- | II-

                for pattern in missed_patterns:

                    # We want to be able to catch genes in .lower() form but not "false genes"
                    if pattern.lower() not in rejected_list:  # remove some false positives

                        matched_pattern = re.compile(pattern, re.IGNORECASE)

                        for group in matched_pattern.finditer(abstract_content):

                            if filename not in dict_caught_annotations:

                                dict_caught_annotations[filename] = []
                                dict_caught_annotations[filename].append(str(group.start() + 1) + '\t' + \
                                                                    str(group.start() + len(group.group()) - 1) + \
                                                                    '\t' + group.group()[1:-1] + '\t' + gene_id + '\n')

                            else:

                                dict_caught_annotations[filename].append(str(group.start() + 1) + '\t' + \
                                                                         str(group.start() + len(group.group()) - 1) + \
                                                                         '\t' + group.group()[1:-1] + '\t' + gene_id + '\n')

                    # Catches genes in capitalized form that match acronyms in .lower() form, but are able to be truly detected in capitalized form
                    elif pattern not in super_rejected_list:  # filter possible true positives discarded as false positives

                        matched_pattern = re.compile(pattern)

                        for group in matched_pattern.finditer(abstract_content):

                            if filename not in dict_caught_annotations:

                                dict_caught_annotations[filename] = []
                                dict_caught_annotations[filename].append(str(group.start() + 1) + '\t' + \
                                                                    str(group.start() + len(group.group()) - 1) + \
                                                                    '\t' + group.group()[1:-1] + '\t' + gene_id + '\n')

                            else:

                                dict_caught_annotations[filename].append(str(group.start()) + '\t' + \
                                                                         str(group.start() + len(group.group()) - 1) + \
                                                                         '\t' + group.group()[1:-1] + '\t' + gene_id + '\n')

            ### PHENOTYPES
            # for hp_name, hp_id in dict_hp_name_id.items():
            #
            #     missed_phenotypes = [' type 2 diabetes ']
            #
            #     # 30266947 - type 2 diabetes
            #
            #     for phenotype in missed_phenotypes:
            #
            #         matched_phenotype = re.compile(phenotype, re.IGNORECASE)
            #
            #         for group in matched_phenoype.finditer(abstract_content):
            #
            #             if filename not in dict_caught_annotations:
            #
            #                 dict_manual_annotations[filename] = []
            #                 dict_caught_annotations[filename].append(str(group.start() + 1) + '\t' + \
            #                                                     str(group.start() + len(group.group()) - 1) + \
            #                                                     '\t' + group.group()[1:-1] + '\t' + hp_id + '\n')
            #
            #             else:
            #
            #                 dict_manual_annotations[filename].append(str(group.start() + 1) + '\t' + \
            #                                                          str(group.start() + len(group.group()) - 1) + \
            #                                                          '\t' + group.group()[1:-1] + '\t' + hp_id + '\n')

    return dict_caught_annotations


#### UPDATE ANNOTATIONS TO DIRECTORY ADDED_ANNOTATIONS/ WITH MANUALLY CORRECTED/JOINED ANNOTATIONS + NEW ANNOTATIONS ####

def update_annotations(corpus_path, data_path, annotations_path, manual_annotations_file, mer_data_path, destination_path):
    """Creates an added annotation file for each abstract in the corpus with the original annotations + manual and/or caught annotations

    :param corpus_path: divided by sentences corpus path
    :param data_path: path to data
    :param annotations_path: divided by sentences annotations path
    :param manual_annotations_file: .tsv divided by sentences annotations to check file
    :param mer_data_path: mer data path
    :param destination_path: destination path
    :return: annotation file for each abstract in the corpus with the original annotations + manual and/or caught annotations
    """

    dict_manual_annotations = manual_annotations(manual_annotations_file)

    try:

        dict_caught_annotations_file = open(data_path + 'dict_caught_annotations.txt', 'r', encoding = 'utf-8')
        dict_caught_annotations = eval(dict_caught_annotations_file.read())
        dict_caught_annotations_file.close()

    except FileNotFoundError:

        dict_caught_annotations = caught_annotations(corpus_path, mer_data_path)
        dict_caught_annotations_file = open(data_path + 'dict_caught_annotations.txt', 'w', encoding = 'utf-8')
        dict_caught_annotations_file.write(str(dict_caught_annotations))
        dict_caught_annotations_file.close()

    for (dir_path, dir_names, file_names) in os.walk(annotations_path):

        for filename in file_names:

            original_annotations = open(annotations_path + filename, 'r', encoding = 'utf-8')
            list_annotations = original_annotations.readlines()
            original_annotations.close()

            added_annotations = open(destination_path + filename, 'w', encoding = 'utf-8')

            if filename in dict_manual_annotations:

                list_annotations = list_annotations + dict_manual_annotations[filename]

            if filename in dict_caught_annotations:

                list_annotations = list_annotations + dict_caught_annotations[filename]

            list_annotations.sort(key=lambda position: int(position.split('\t', 1)[0]))

            next_position = 1

            for position in range(len(list_annotations)):

                last_position = int(list_annotations[position].split('\t')[1])

                if next_position < len(list_annotations):

                    first_next_position = int(list_annotations[next_position].split('\t')[0])

                    # Some annotations overlap, we chose to keep the larger ones because that usually means more specification
                    # Practical example: between the genes NAF and NAF-1 we opted by NAF-1, a specification of the gene NAF

                    if first_next_position <= last_position:

                        if len(list_annotations[position].split('\t')[2]) > len(list_annotations[next_position].split('\t')[2]):
                            added_annotations.write(list_annotations[position])

                    else:
                        added_annotations.write(list_annotations[position])

                else:
                    added_annotations.write(list_annotations[position])

                next_position += 1

            added_annotations.close()

    return


#### REMOVE ABSTRACTS + ANNOTATIONS WITH ONLY ANNOTATED GENES/PHENOTYPES OR WITHOUT ANNOTATIONS ####

def final_annotations(annotations_path, corpus_path, destination_corpus_path, destination_path):
    """Creates a final annotation file with annotations of genes and phenotypes
       and a final text file, for each abstract in the corpus

    :param annotations_path: added annotations path
    :param corpus_path: divided by sentences corpus path
    :param destination_corpus_path: destination corpus path
    :param destination_path: destination path
    :return: annotation file with annotations of genes and phenotypes and a
             text file, for each abstract in the corpus
    """

    file_to_save_list = []  # final_text directory

    for (dir_path, dir_names, file_names) in os.walk(annotations_path):

        for filename in file_names:

            added_annotations = open(annotations_path + filename, 'r', encoding = 'utf-8')
            list_annotations = added_annotations.readlines()
            added_annotations.close()

            if list_annotations != []:

                check_phenotype = False
                check_gene = False

                for annotation in list_annotations:

                    if 'HP_' in annotation:
                        check_phenotype = True

                    elif annotation.split('\t')[-1][:2] != 'HP':
                        check_gene = True

                if check_phenotype and check_gene:

                    final_annotations = open(destination_path + filename, 'w', encoding = 'utf-8')

                    for annotation in list_annotations:
                        final_annotations.write(annotation)

                    final_annotations.close()

                    file_to_save_list.append(filename)

    for file_to_save in file_to_save_list:  # just to save the abstracts that are going to be "used"

        os.system('cp ' + corpus_path + file_to_save + ' ' + destination_corpus_path)

    return


#### CREATES A FILE WITH THE RELATIONS BETWEEN GENE AND HP ANNOTATIONS FOR EACH SENTENCE IN EACH ABSTRACT ####

def relations_annotations(corpus_path, file_g2p, file_p2g, annotations_path, destination_file_path):
    """Creates a file with the relations between gene and phenotype annotations for each sentence in each abstract of the corpus
       and counts the number of true and false relations between gene and phenotype annotations

    :param corpus_path: divided by sentences corpus path
    :param file_g2p: file with relations gene to phenotype
    :param file_p2g: file with relations phenotype to gene
    :param annotations_path: final annotations path
    :param destination_file_path: destination file path
    :return: file with the relations between gene and phenotype annotations for each sentence in each abstract of the corpus
             and a string with the number of true and false relations between gene and phenotype annotations
    """

    dict_pg = relations.join_dicts(file_g2p, file_p2g, dict_type = 1)

    relations_file = open(destination_file_path, 'w', encoding = 'utf-8')
    relations_file.write('FILE_ID\tSENTENCE\tGENE\tPHENOTYPE\tGENE_ID\tPHENOTYPE_ID\tGENE_START_POSITION\tGENE_END_POSITION\tPHENOTYPE_START_POSITION\tPHENOTYPE_END_POSITION\tRELATION\n')

    for (dir_path, dir_names, file_names) in os.walk(annotations_path):

        for filename in file_names:

            abstract = open(corpus_path + filename, 'r', encoding = 'utf-8')
            abstract_content = abstract.readlines()
            abstract.close()

            final_annotations = open(annotations_path + filename, 'r', encoding = 'utf-8')
            list_annotations = final_annotations.readlines()
            final_annotations.close()

            for sentence in abstract_content:

                sentence_position_count = sentence

                list_genes_in_sentence = []
                list_phenotypes_in_sentence = []

                lenght = len(sentence[:-1])

                for annotation in list_annotations:

                    first_position = str(len(sentence_position_count.split(annotation.split('\t')[2], 1)[0]))
                    last_position = str(len(sentence_position_count.split(annotation.split('\t')[2], 1)[0]) + len(annotation.split('\t')[2]))
                    sentence_position_count = sentence_position_count.replace(annotation.split('\t')[2], len(annotation.split('\t')[2]) * 'A')

                    if int(last_position) < lenght:

                        if 'HP_' in annotation:

                            list_phenotypes_in_sentence.append((annotation.split('\t')[2], annotation.split('\t')[3][:-1],
                                                                first_position, last_position))

                        else:

                            list_genes_in_sentence.append((annotation.split('\t')[2], annotation.split('\t')[3][:-1],
                                                           first_position, last_position))

                if list_phenotypes_in_sentence != [] and list_genes_in_sentence != []:

                    for phenotype_name, phenotype_id, phenotype_first, phenotype_last in list_phenotypes_in_sentence:

                        for gene_name, gene_id, gene_first, gene_last in list_genes_in_sentence:

                            if gene_id in dict_pg and phenotype_id in dict_pg[gene_id]:
                                relation = 'True'

                            else:
                                relation = 'False'

                            relations_file.write(filename + '\t' + sentence[:-1] + '\t' + gene_name + '\t' + phenotype_name + \
                                                 '\t' + gene_id + '\t' + phenotype_id + '\t' + gene_first + '\t' + gene_last + '\t' \
                                                 + phenotype_first + '\t' + phenotype_last + '\t' + relation + '\n')

    relations_file.close()

    return


#### SEPARATES THE RELATIONS IN NOT_VERIFIED (ONE FILE) AND VERIFIED (TEN FILES WITH 50 RELATIONS EACH WITH AN OVERLAP OF 20 RELATIONS (320 RELATIONS TOTAL)) ####

def verify_relations_annotations(annotations_path):
    """Separates the relations in not_verified (one file) and verified (eight files with 50 relations each
       with an overlap of 20 relations (260 relations total))

    :param annotations_path: relations annotations path
    :return: .tsv files, not_verified (one file) and verified (eight files with 50 relations each with an
             overlap of 20 relations (260 relations total))
    """

    relations_annotations = open(annotations_path + 'relations.tsv', 'r', encoding = 'utf-8')
    header = relations_annotations.readline()
    list_annotations = relations_annotations.readlines()
    relations_annotations.close()

    lucky_twenty = random.sample(list_annotations, 20)

    list_annotations = [annotation for annotation in list_annotations if annotation not in lucky_twenty]

    for i in range(8):
        verify_file = open(annotations_path + 'verify_' + str(i) + '.tsv', 'w', encoding = 'utf-8')
        verify_file.write(header[:-1] + '\t' + 'CONFIRMATION (CORRECT(C) | INCORRECT(I) | UNCERTAIN(U))' + '\n')

        for chosen_annotation in lucky_twenty:

            verify_file.write(chosen_annotation)

        random_thirty = random.sample(list_annotations, 30)

        for random_annotation in random_thirty:

            verify_file.write(random_annotation)

        list_annotations = [annotation for annotation in list_annotations if annotation not in random_thirty]

        verify_file.close()

    not_verify_file = open(annotations_path + 'not_verify.tsv', 'w', encoding = 'utf-8')
    not_verify_file.write(header)

    for excluded_annotation in list_annotations:

        not_verify_file.write(excluded_annotation)

    not_verify_file.close()

    return


#### RUN ####

def main():
    """Creates a directory with a file for each retrieved abstract with their respective
       gene and human phenotype annotations and removes from the pubmed directory
       abstracts without annotations from both types (human phenotypes and genes)

    :return: directory with a file for each retrieved abstract with their respective
             gene and human phenotype annotations and removing, from the pubmed
             directory, of abstracts without annotations from both types (human phenotypes and genes)
    """

    os.system('mkdir -p corpora/per_directory_text || true')
    os.system('mkdir -p corpora/per_directory_annotations || true')
    os.system('mkdir -p corpora/divided_by_sentences_annotations || true')
    os.system('mkdir -p corpora/added_annotations || true')
    os.system('mkdir -p corpora/gene_phenotype_annotations || true')
    os.system('mkdir -p corpora/new_corpus || true')
    number_of_directories = divide_directory('corpora/pubmed_corpus/', 'corpora/per_directory_text/')

    for i in range(number_of_directories):

        os.system('cp corpora/per_directory_text/directory_' + str(i) + '/* bin/IHP/corpora/hpo/test_corpus/')
        os.chdir('bin/IHP/')
        os.system('python3 src/main.py load_corpus --goldstd hpo_test --log DEBUG')
        os.system('python3 src/main.py test --goldstd hpo_test -o pickle data/results_hpo_train --models models/hpo_train --log DEBUG')
        os.system('python3 src/evaluate.py evaluate hpo_test --results data/results_hpo_train --models models/hpo_train --log DEBUG')
        os.system('mv data/results_hpo_train_report.txt ../../corpora/per_directory_annotations/directory_' + str(i) + '.txt')
        os.chdir('../..')
        os.system('rm bin/IHP/corpora/hpo/test_corpus/* || true')

    join_report_files('corpora/per_directory_annotations/')
    os.chdir('bin/MER/')
    os.system('cd data; ../produce_data_files.sh genes.txt')
    annotations('../../corpora/pubmed_corpus/', 'get_entities.sh','../../corpora/per_directory_annotations/final_report', 'data/hp_links.tsv', '../../corpora/')
    os.chdir('../..')
    update_annotations('corpora/pubmed_corpus/', 'data/', 'corpora/divided_by_sentences_annotations/','corpora/annotations_to_check.tsv', 'bin/MER/data/', 'corpora/added_annotations/')
    final_annotations('corpora/added_annotations/', 'corpora/pubmed_corpus/', 'corpora/new_corpus/', 'corpora/gene_phenotype_annotations/')
    os.system('rm corpora/pubmed_corpus/* || true')
    os.system('mv corpora/new_corpus/* corpora/pubmed_corpus/ || true')
    relations_annotations('corpora/pubmed_corpus/', 'data/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt', 'data/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt', 'corpora/gene_phenotype_annotations/', 'corpora/relations.tsv')
    #verify_relations_annotations('corpora/')  # if curator confirmation for a test corpus
    os.system('rm -rf corpora/per_directory_text corpora/per_directory_annotations corpora/divided_by_sentences_annotations corpora/added_annotations corpora/new_corpus || true')

    return


if __name__ == "__main__":
    main()
