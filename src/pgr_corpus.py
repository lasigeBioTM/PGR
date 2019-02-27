import os
import collections
import sys
import re

import relations


#### GENE - GO DICTIONARY: GENE ID - (GO ID, EVIDENCE CODE, GO NAME, CATEGORY), (...) ####

def dict_g2go(file_g2go):
    """Creates a dictionary of type {gene1:[(GO_ID, Evidence, GO_name, category),
    (GO_ID, Evidence, GO_name, category), ...], }

    :param file_g2go: file with relations gene to GO
    :return: dict of type {gene1:[(GO_ID, Evidence, GO_name, category),
             (GO_ID, Evidence, GO_name, category), ...], }
    """

    os.system('gunzip -k ' + file_g2go + '.gz')

    gene2go = open(file_g2go, 'r', encoding = 'utf-8')

    gene2go.readline()  # skip header

    relations_g2go = gene2go.readlines()
    gene2go.close()

    relations_g2go.pop()

    dict_gene_go = {}

    for line in relations_g2go:

        line = line.split('\t')

        gene_id = line[1]
        go = line[2]
        evidence = line[3]
        name = line[5]
        category = line[7][:-1]

        if gene_id not in dict_gene_go:
            dict_gene_go[gene_id] = []
            dict_gene_go[gene_id].append((go, evidence, name, category))

        else:
            dict_gene_go[gene_id].append((go, evidence, name, category))

    os.system('rm ' + file_g2go)

    return dict_gene_go


#### REPLACEMENT OF GENE ANNOTATIONS IN DIVIDED BY SENTENCES ANNOTATIONS FOR THEIR MOST REPRESENTATIVE GO ANNOTATION TO DIVIDED BY SENTENCES GO ANNOTATIONS ####

def go_annotations(annotations_path, file_g2go, destination_path):
    """Generates a file for each abstract with the correspondent phenotype and GO annotations,
       creates a dictionary of type {gene_id:go_id, gene_id2:go_id, } and
       a dictionary of type {gene_name:go_name, gene_name:go_name, }

    :param annotations_path: divided by sentences annotations path
    :param file_g2go: file with relations gene to GO
    :param destination_path: destination path
    :return: file for each abstract with the correspondent phenotype and GO annotations,
             creates a dictionary of type {gene_id:go_id, gene_id2:go_id, } and
             a dictionary of type {gene_name:go_name, gene_name:go_name, }
             annotation file example:

             26 29  negative regulation of cell proliferation   GO_0008285
             279	288	bilateral	HP_0012832
             313	323	unilateral	HP_0012833

    """

    dict_gene_id_go = dict_g2go(file_g2go)
    dict_gene_go_id = {}
    dict_gene_go_name = {}

    for (dir_path, dir_names, file_names) in os.walk(annotations_path):

        for filename in file_names:

            annotation_file = open(annotations_path + filename, 'r', encoding = 'utf-8')
            contents = annotation_file.readlines()

            annotation_file.close()

            annotation_file_go = open(destination_path + filename, 'w', encoding = 'utf-8')

            save = 0

            for line in contents:
                line = line.split('\t')

                start = int(line[0])
                end = int(line[1])

                if save != 0:
                    start = start + save
                    end = end + save

                if line[3].startswith('HP'):
                    annotation_file_go.write(str(start) + '\t' + str(end) + '\t' + line[2] + '\t' + line[3])

                else:
                    value = False
                    for key_g, value_tup in dict_gene_id_go.items():

                        if line[3][:-1] == key_g:

                            list_evidence = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP',
                                             'HDA', 'HMP', 'HGI', 'HEP', 'ISS', 'ISO', 'ISA',
                                             'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA',
                                             'TAS', 'NAS', 'IC', 'ND', 'IEA']  # order criteria

                            bins = collections.defaultdict(list)
                            for pair in value_tup:
                                bins[pair[1]].append(pair)

                            value_tup = [pair for i in list_evidence for pair in bins[i]]

                            for v_tup in value_tup:

                                if v_tup[3] == 'Process':  # Biological Process

                                    value = True
                                    end = int(line[1]) + len(v_tup[2]) - len(line[2]) + save
                                    save = save + len(v_tup[2]) - len(line[2])

                                    dict_gene_go_id[line[3][:-1]] = v_tup[0].replace(':', '_')
                                    dict_gene_go_name[line[2]] = v_tup[2]
                                    annotation_file_go.write(str(start) + '\t' + str(end) + '\t' + \
                                                             v_tup[2] + '\t' + v_tup[0].replace(':', '_') + '\n')

                                    break

                    if not value:  # genes with no associated GO terms or no associated GO terms from Biological Processes

                        end = int(line[1]) + len('biological_process') - len(line[2]) + save
                        save = save + len('biological_process') - len(line[2])

                        dict_gene_go_id[line[3][:-1]] = 'GO_0008150'
                        dict_gene_go_name[line[2]] = 'biological_process'
                        annotation_file_go.write(str(start) + '\t' + str(end) + '\t' + \
                                                 'biological_process' + '\t' + 'GO_0008150' + '\n')

            annotation_file_go.close()

    return dict_gene_go_id, dict_gene_go_name


#### XML FORMAT GENE AND PHENOTYPE CORPUS ####

def pgr_gene(verify_file, destination_path, type = None):
    """Generates a .xml file for each abstract with sentences with relations in corpus with the correspondent phenotype and gene annotations

    :param verify_file: file with sentences with relations verified (for test corpus)
                        or file with sentences with relations not verified (for train corpus)
                        or file with all the relations (for corpus without curator correction)
    :param destination_path: destination path
    :param type: type (optional) if pretended file is a test corpus file
    :return: .xml file for each abstract with sentences with relations in corpus with the correspondent phenotype and gene annotations of type:

    <sentence id="s0" text="In addition, the coexistence of high MACC1 and low NM23-H1 expression and tumor budding
    was associated with short OS (p AAAA 0.001).">
		<entity id="s0.e1" charOffset="51-55"
			type="GENE" text="NM23" ontology_id="4830"/>
		<entity id="s0.e2" charOffset="74-79"
			type="HP" text="tumor" ontology_id="HP_0002664"/>
		<pair id="s0.p1" e1="s0.e1"
		    e2="s0.e2" pgr="true"/>
	</sentence>

    """

    verify = open(verify_file, 'r', encoding = 'utf-8')

    verify.readline()  # skip header

    verify_relations = [line.split('\t') for line in verify]
    verify.close()
    verify_relations.sort(key=lambda x: int(x[0]))  # sort by abstract identifier

    iterator = 1
    sentence_number = 1
    entity_number = 1
    pair_number = 1

    dict_entities = {}
    dict_pairs = {}

    for line in verify_relations:

        abstract = line[0]
        sentence = line[1]
        gene = line[2]
        phenotype = line[3]
        gene_id = line[4]
        phenotype_id = line[5]
        gene_start_position = line[6]
        gene_end_position = line[7]
        phenotype_start_position = line[8]
        phenotype_end_position = line[9]

        if type:
            relation = line[10]

        else:
            relation = line[10][:-1]

        if verify_relations[iterator - 2][0] == abstract:  # same abstract

            if verify_relations[iterator - 2][1] == sentence:  # same sentence

                if int(gene_start_position) < int(phenotype_start_position):

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + gene_start_position + '-' + \
                                 gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
                                 + gene + '" ontology_id="' + gene_id + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                                 phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                                 + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                else:

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                                 phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                                 + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + gene_start_position + '-' + \
                                 gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
                                 + gene + '" ontology_id="' + gene_id + '"/>' + '\n'))

                    entity_number += 1

                dict_pairs[pair_number] = []
                dict_pairs[pair_number].append(((entity_number - 2, entity_number - 1),'\t\t' + '<pair id="s' + str(sentence_number) + '.p' + str(pair_number) + '" e1="s' + str(sentence_number) + \
                             '.e' + str(entity_number - 2) + '"\n\t\t    ' + 'e2="s' + str(sentence_number) + '.e' + str(entity_number - 1) + '" pgr="' + relation.lower() + '"/>' + '\n'))

                pair_number += 1

            else:  # different sentence

                pair_number = 1
                entity_number = 1

                list_entities = sorted(dict_entities.items())

                used_entities_list = []
                used_numbers_list = []
                to_write_entities = []
                right_number = 1
                save_alterations = {}

                for element in range(1, len(list_entities) + 1):

                    if list_entities[element - 1][1][0][0] not in used_entities_list:

                        to_write_entities.append(str(list_entities[element - 1][1][0][1]).replace('e' + str(list_entities[element - 1][0]), 'e' + str(right_number)))
                        used_entities_list.append(list_entities[element - 1][1][0][0])
                        used_numbers_list.append((list_entities[element - 1][1][0][0], element))
                        save_alterations['e' + str(list_entities[element - 1][0])] = 'e' + str(right_number)
                        right_number += 1

                    else:

                        for used_number in used_numbers_list:
                            if used_number[0] == list_entities[element - 1][1][0][0]:
                                save_alterations['e' + str(element)] = 'e' + str(used_number[1])

                organized_writing = []
                for line_to_write in to_write_entities:
                    first_offset = int(line_to_write.split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[0])
                    organized_writing.append((first_offset, line_to_write))

                organized_writing = sorted(organized_writing, key=lambda tup: tup[0])

                new_entity_number = 1
                used_keys = []

                for organized_tuple in organized_writing:
                    original_entity_number = int(organized_tuple[1].split('.e')[1].split('" charOffset="')[0])
                    writer.write(re.sub(r'.e[0-9]+', '.e' + str(new_entity_number), organized_tuple[1]))

                    for key, value in save_alterations.items():
                        if value == 'e' + str(original_entity_number) and key not in used_keys:
                            save_alterations[key] = 'e' + str(new_entity_number)
                            used_keys.append(key)

                    new_entity_number += 1

                dict_entities = {}

                list_pairs = sorted(dict_pairs.items())

                for pair in list_pairs:

                    writer.write(str(pair[1][0][1].replace('.e' + str(pair[1][0][0][0]), '.' + save_alterations['e' + str(pair[1][0][0][0])]).replace('.e' + str(pair[1][0][0][1]), '.' + save_alterations['e' + str(pair[1][0][0][1])])))

                dict_pairs = {}
                writer.write('\t' + '</sentence>' + '\n')
                sentence_number += 1
                sentence = sentence.replace(' <', ' l').replace('(<', '(l').replace('(p<', '(pl').replace(' < ', ' l ').replace('.&quot', '.AAAAA').replace('&gt;', 'AAAA').replace('&quot;', 'AAAAAA').replace('&lt;','AAAA').replace('&amp;', 'AAAAA').split('\n')[0]  # avoid invalid (bad/not well-formed) XML
                writer.write('\t' + '<sentence id="s' + str(sentence_number) + '" text="' + sentence + '">' + '\n')

                if int(gene_start_position) < int(phenotype_start_position):

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene,
                        '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + gene_start_position + '-' + \
                        gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
                        + gene + '" ontology_id="' + gene_id + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype,
                        '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                        phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                        + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                else:

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype,
                        '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                        phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                        + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene,
                        '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + gene_start_position + '-' + \
                        gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
                        + gene + '" ontology_id="' + gene_id + '"/>' + '\n'))

                    entity_number += 1

                dict_pairs[pair_number] = []
                dict_pairs[pair_number].append(((entity_number - 2, entity_number - 1),'\t\t' + '<pair id="s' + str(sentence_number) + '.p' + str(pair_number) + '" e1="s' + str(sentence_number) + \
                             '.e' + str(entity_number - 2) + '"\n\t\t    ' + 'e2="s' + str(sentence_number) + '.e' + str(entity_number - 1) + '" pgr="' + relation.lower() + '"/>' + '\n'))

                pair_number += 1

        else:  # different abstract

            if iterator != 1:  # different for first one
                pair_number = 1
                entity_number = 1
                sentence_number = 1

                list_entities = sorted(dict_entities.items())

                used_entities_list = []
                used_numbers_list = []
                to_write_entities = []
                right_number = 1
                save_alterations = {}

                for element in range(1, len(list_entities) + 1):

                    if list_entities[element - 1][1][0][0] not in used_entities_list:

                        to_write_entities.append(str(list_entities[element - 1][1][0][1]).replace('e' + str(list_entities[element - 1][0]), 'e' + str(right_number)))
                        used_entities_list.append(list_entities[element - 1][1][0][0])
                        used_numbers_list.append((list_entities[element - 1][1][0][0], element))
                        save_alterations['e' + str(list_entities[element - 1][0])] = 'e' + str(right_number)
                        right_number += 1

                    else:

                        for used_number in used_numbers_list:
                            if used_number[0] == list_entities[element - 1][1][0][0]:
                                save_alterations['e' + str(element)] = 'e' + str(used_number[1])

                organized_writing = []
                for line_to_write in to_write_entities:
                    first_offset = int(line_to_write.split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[0])
                    organized_writing.append((first_offset, line_to_write))

                organized_writing = sorted(organized_writing, key=lambda tup: tup[0])

                new_entity_number = 1
                used_keys = []

                for organized_tuple in organized_writing:
                    original_entity_number = int(organized_tuple[1].split('.e')[1].split('" charOffset="')[0])
                    writer.write(re.sub(r'.e[0-9]+', '.e' + str(new_entity_number), organized_tuple[1]))

                    for key, value in save_alterations.items():
                        if value == 'e' + str(original_entity_number) and key not in used_keys:
                            save_alterations[key] = 'e' + str(new_entity_number)
                            used_keys.append(key)

                    new_entity_number += 1

                dict_entities = {}

                list_pairs = sorted(dict_pairs.items())

                for pair in list_pairs:

                    writer.write(str(pair[1][0][1].replace('.e' + str(pair[1][0][0][0]), '.' + save_alterations['e' + str(pair[1][0][0][0])]).replace('.e' + str(pair[1][0][0][1]), '.' + save_alterations['e' + str(pair[1][0][0][1])])))

                dict_pairs = {}
                writer.write('\t' + '</sentence>' + '\n')
                writer.write('</document>' + '\n')
                writer.close()

            writer = open(destination_path + abstract + '.xml', 'w', encoding = 'utf-8')
            writer.write('<?xml version="1.0" encoding="UTF-8"?>' + '\n')
            writer.write('<document id="' + abstract + '">' + '\n')

            sentence = sentence.replace(' <', ' l').replace('(<', '(l').replace('(p<', '(pl').replace(' < ', ' l ').replace('.&quot', '.AAAAA').replace('&gt;', 'AAAA').replace('&quot;', 'AAAAAA').replace('&lt;', 'AAAA').replace('&amp;', 'AAAAA').split('\n')[0]  # avoid invalid (bad/not well-formed) XML
            writer.write('\t' + '<sentence id="s' + str(sentence_number) + '" text="' + sentence + '">' + '\n')

            if int(gene_start_position) < int(phenotype_start_position):

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + gene_start_position + '-' + \
                             gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
                             + gene + '" ontology_id="' + gene_id + '"/>' + '\n'))

                entity_number += 1

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((phenotype,
                    '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                    phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                    + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                entity_number += 1

            else:

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((phenotype,
                    '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                    phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                    + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                entity_number += 1

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + gene_start_position + '-' + \
                             gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
                             + gene + '" ontology_id="' + gene_id + '"/>' + '\n'))

                entity_number += 1

            dict_pairs[pair_number] = []
            dict_pairs[pair_number].append(((entity_number - 2, entity_number - 1),'\t\t' + '<pair id="s' + str(sentence_number) + '.p' + str(pair_number) + '" e1="s' + str(sentence_number) + \
                         '.e' + str(entity_number - 2) + '"\n\t\t    ' + 'e2="s' + str(sentence_number) + '.e' + str(entity_number - 1) + '" pgr="' + relation.lower() + '"/>' + '\n'))

            pair_number += 1

        iterator += 1

    return


#### XML FORMAT GO AND PHENOTYPE CORPUS ####

def pgr_go(file_g2go, go_annotations_path, annotations_path, verify_file, destination_path, type = None):
    """Generates a .xml file for each abstract with sentences with relations in corpus with the correspondent phenotype and GO annotations

    :param file_g2go: file with relations gene to GO
    :param go_annotations_path: divided by sentences go annotations path
    :param annotations_path: final annotations path
    :param verify_file: file with sentences with relations verified (for test corpus)
                        or file with sentences with relations not verified (for train corpus)
                        or file with all the relations (for corpus without curator correction)
    :param destination_path: destination path
    :param type: type (optional) if pretended file is a test corpus file
    :return: .xml file for each abstract with sentences with relations in corpus with the correspondent phenotype and GO annotations of type:

    <sentence id="s0" text="In addition, the coexistence of high MACC1 and low positive regulation of DNA binding-H1
    expression and tumor budding was associated with short OS (p AAAA 0.001).">
		<entity id="s0.e1" charOffset="51-85"
			type="GO" text="positive regulation of DNA binding" ontology_id="GO_0043388"/>
		<entity id="s0.e2" charOffset="104-109"
			type="HP" text="tumor" ontology_id="HP_0002664"/>
		<pair id="s0.p1" e1="s0.e1"
		    e2="s0.e2" pgr="true"/>
	</sentence>

    """

    dict_g2go_id, dict_g2go_name = go_annotations(annotations_path, file_g2go, go_annotations_path)

    verify = open(verify_file, 'r', encoding = 'utf-8')

    verify.readline()  # skip header

    verify_relations = [line.split('\t') for line in verify]
    verify.close()
    verify_relations.sort(key=lambda x: int(x[0]))  # sort by abstract identifier

    iterator = 1
    sentence_number = 1
    entity_number = 1
    pair_number = 1

    dict_entities = {}
    dict_pairs = {}
    save_sentence = ''

    for line in verify_relations:

        abstract = line[0]
        sentence = line[1]
        gene = line[2]
        phenotype = line[3]
        gene_id = line[4]
        phenotype_id = line[5]
        gene_start_position = line[6]
        gene_end_position = line[7]
        phenotype_start_position = line[8]
        phenotype_end_position = line[9]

        if type:
            relation = line[10]

        else:
            relation = line[10][:-1]

        if verify_relations[iterator - 2][0] == abstract:  # same abstract

            if verify_relations[iterator - 2][1] == sentence:  # same sentence

                go_start_position = gene_start_position
                go_end_position = str(int(gene_start_position) + len(dict_g2go_name[gene]))

                if int(gene_start_position) < int(phenotype_start_position):

                    phenotype_start_position = str(int(phenotype_start_position) + len(dict_g2go_name[gene]) - len(gene))
                    phenotype_end_position = str(int(phenotype_end_position) + len(dict_g2go_name[gene]) - len(gene))

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + go_start_position + '-' + \
                             go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
                             + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                                 phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                                 + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                else:

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                                 phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                                 + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + go_start_position + '-' + \
                             go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
                             + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n'))

                    entity_number += 1

                dict_pairs[pair_number] = []
                dict_pairs[pair_number].append(((entity_number - 2, entity_number - 1),'\t\t' + '<pair id="s' + str(sentence_number) + '.p' + str(pair_number) + '" e1="s' + str(sentence_number) + \
                             '.e' + str(entity_number - 2) + '"\n\t\t    ' + 'e2="s' + str(sentence_number) + '.e' + str(entity_number - 1) + '" pgr="' + relation.lower() + '"/>' + '\n'))

                pair_number += 1

            else:  # different sentence

                pair_number = 1
                entity_number = 1

                for entity_n, info_tuple in dict_entities.items():

                    first_offset = info_tuple[0][1].split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[0]
                    second_offset = info_tuple[0][1].split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[1]

                    if 'GO' in info_tuple[0][1]:
                        go_start_position = str(len(save_sentence.split(info_tuple[0][0], 1)[0]))
                        go_end_position = str(int(go_start_position) + len(dict_g2go_name[info_tuple[0][0]]))

                        dict_entities[entity_n] = [(info_tuple[0][0],
                                                    info_tuple[0][1].replace(first_offset, go_start_position).replace(
                                                        second_offset, go_end_position))]

                        save_sentence = save_sentence.replace(info_tuple[0][0], dict_g2go_name[info_tuple[0][0]], 1)

                writer.write('\t' + '<sentence id="s' + str(sentence_number) + '" text="' + save_sentence + '">' + '\n')

                list_entities = sorted(dict_entities.items())

                used_entities_list = []
                used_numbers_list = []
                to_write_entities = []
                right_number = 1
                save_alterations = {}

                for element in range(1, len(list_entities) + 1):

                    if list_entities[element - 1][1][0][0] not in used_entities_list:

                        to_write_entities.append(str(list_entities[element - 1][1][0][1]).replace('e' + str(list_entities[element - 1][0]), 'e' + str(right_number)))
                        used_entities_list.append(list_entities[element - 1][1][0][0])
                        used_numbers_list.append((list_entities[element - 1][1][0][0], element))
                        save_alterations['e' + str(list_entities[element - 1][0])] = 'e' + str(right_number)
                        right_number += 1

                    else:

                        for used_number in used_numbers_list:
                            if used_number[0] == list_entities[element - 1][1][0][0]:
                                save_alterations['e' + str(element)] = 'e' + str(used_number[1])

                organized_writing = []
                for line_to_write in to_write_entities:
                    first_offset = int(line_to_write.split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[0])
                    organized_writing.append((first_offset, line_to_write))

                organized_writing = sorted(organized_writing, key=lambda tup: tup[0])

                new_entity_number = 1
                used_keys = []

                for organized_tuple in organized_writing:
                    original_entity_number = int(organized_tuple[1].split('.e')[1].split('" charOffset="')[0])
                    writer.write(re.sub(r'.e[0-9]+', '.e' + str(new_entity_number), organized_tuple[1]))

                    for key, value in save_alterations.items():
                        if value == 'e' + str(original_entity_number) and key not in used_keys:
                            save_alterations[key] = 'e' + str(new_entity_number)
                            used_keys.append(key)

                    new_entity_number += 1

                dict_entities = {}

                list_pairs = sorted(dict_pairs.items())

                for pair in list_pairs:

                    writer.write(str(pair[1][0][1].replace('.e' + str(pair[1][0][0][0]), '.' + save_alterations['e' + str(pair[1][0][0][0])]).replace('.e' + str(pair[1][0][0][1]), '.' + save_alterations['e' + str(pair[1][0][0][1])])))

                dict_pairs = {}
                save_sentence = ''
                writer.write('\t' + '</sentence>' + '\n')
                sentence_number += 1
                sentence = sentence.replace(' <', ' l').replace('(<', '(l').replace('(p<', '(pl').replace(' < ', ' l ').replace('.&quot', '.AAAAA').replace('&gt;', 'AAAA').replace('&quot;', 'AAAAAA').replace('&lt;','AAAA').replace('&amp;', 'AAAAA').split('\n')[0]  # avoid invalid (bad/not well-formed) XML
                save_sentence = sentence

                go_start_position = gene_start_position
                go_end_position = str(int(gene_start_position) + len(dict_g2go_name[gene]))

                if int(gene_start_position) < int(phenotype_start_position):

                    phenotype_start_position = str(int(phenotype_start_position) + len(dict_g2go_name[gene]) - len(gene))
                    phenotype_end_position = str(int(phenotype_end_position) + len(dict_g2go_name[gene]) - len(gene))

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + go_start_position + '-' + \
                             go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
                             + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype,
                        '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                        phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                        + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                else:

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((phenotype,
                        '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                        phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                        + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                    entity_number += 1

                    dict_entities[entity_number] = []
                    dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + go_start_position + '-' + \
                             go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
                             + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n'))

                    entity_number += 1

                dict_pairs[pair_number] = []
                dict_pairs[pair_number].append(((entity_number - 2, entity_number - 1),'\t\t' + '<pair id="s' + str(sentence_number) + '.p' + str(pair_number) + '" e1="s' + str(sentence_number) + \
                             '.e' + str(entity_number - 2) + '"\n\t\t    ' + 'e2="s' + str(sentence_number) + '.e' + str(entity_number - 1) + '" pgr="' + relation.lower() + '"/>' + '\n'))

                pair_number += 1

        else:  # different abstract

            if iterator != 1:  # different for first one
                pair_number = 1
                entity_number = 1

                for entity_n, info_tuple in dict_entities.items():

                    first_offset = info_tuple[0][1].split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[0]
                    second_offset = info_tuple[0][1].split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[1]

                    if 'GO' in info_tuple[0][1]:

                        go_start_position = str(len(save_sentence.split(info_tuple[0][0], 1)[0]))
                        go_end_position = str(int(go_start_position) + len(dict_g2go_name[info_tuple[0][0]]))

                        dict_entities[entity_n] = [(info_tuple[0][0], info_tuple[0][1].replace(first_offset, go_start_position).replace(second_offset, go_end_position))]

                        save_sentence = save_sentence.replace(info_tuple[0][0], dict_g2go_name[info_tuple[0][0]], 1)

                writer.write('\t' + '<sentence id="s' + str(sentence_number) + '" text="' + save_sentence + '">' + '\n')

                sentence_number = 1
                list_entities = sorted(dict_entities.items())

                used_entities_list = []
                used_numbers_list = []
                to_write_entities = []
                right_number = 1
                save_alterations = {}

                for element in range(1, len(list_entities) + 1):

                    if list_entities[element - 1][1][0][0] not in used_entities_list:

                        to_write_entities.append(str(list_entities[element - 1][1][0][1]).replace('e' + str(list_entities[element - 1][0]), 'e' + str(right_number)))
                        used_entities_list.append(list_entities[element - 1][1][0][0])
                        used_numbers_list.append((list_entities[element - 1][1][0][0], element))
                        save_alterations['e' + str(list_entities[element - 1][0])] = 'e' + str(right_number)
                        right_number += 1

                    else:

                        for used_number in used_numbers_list:
                            if used_number[0] == list_entities[element - 1][1][0][0]:
                                save_alterations['e' + str(element)] = 'e' + str(used_number[1])

                organized_writing = []
                for line_to_write in to_write_entities:
                    first_offset = int(line_to_write.split('charOffset="')[1].split('"\n\t\t\t')[0].split('-')[0])
                    organized_writing.append((first_offset, line_to_write))

                organized_writing = sorted(organized_writing, key=lambda tup: tup[0])

                new_entity_number = 1
                used_keys = []

                for organized_tuple in organized_writing:
                    original_entity_number = int(organized_tuple[1].split('.e')[1].split('" charOffset="')[0])
                    writer.write(re.sub(r'.e[0-9]+', '.e' + str(new_entity_number), organized_tuple[1]))

                    for key, value in save_alterations.items():
                        if value == 'e' + str(original_entity_number) and key not in used_keys:
                            save_alterations[key] = 'e' + str(new_entity_number)
                            used_keys.append(key)

                    new_entity_number += 1

                dict_entities = {}

                list_pairs = sorted(dict_pairs.items())

                for pair in list_pairs:

                    writer.write(str(pair[1][0][1].replace('.e' + str(pair[1][0][0][0]), '.' + save_alterations['e' + str(pair[1][0][0][0])]).replace('.e' + str(pair[1][0][0][1]), '.' + save_alterations['e' + str(pair[1][0][0][1])])))

                dict_pairs = {}
                save_sentence = ''
                writer.write('\t' + '</sentence>' + '\n')
                writer.write('</document>' + '\n')
                writer.close()

            writer = open(destination_path + abstract + '.xml', 'w', encoding = 'utf-8')
            writer.write('<?xml version="1.0" encoding="UTF-8"?>' + '\n')
            writer.write('<document id="' + abstract + '">' + '\n')

            sentence = sentence.replace(' <', ' l').replace('(<', '(l').replace('(p<', '(pl').replace(' < ', ' l ').replace('.&quot', '.AAAAA').replace('&gt;', 'AAAA').replace('&quot;', 'AAAAAA').replace('&lt;', 'AAAA').replace('&amp;', 'AAAAA').split('\n')[0]  # avoid invalid (bad/not well-formed) XML
            save_sentence = sentence

            go_start_position = gene_start_position
            go_end_position = str(int(gene_start_position) + len(dict_g2go_name[gene]))

            if int(gene_start_position) < int(phenotype_start_position):

                phenotype_start_position = str(int(phenotype_start_position) + len(dict_g2go_name[gene]) - len(gene))
                phenotype_end_position =  str(int(phenotype_end_position) + len(dict_g2go_name[gene]) - len(gene))

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + go_start_position + '-' + \
                             go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
                             + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n'))

                entity_number += 1

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((phenotype,
                    '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                    phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                    + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                entity_number += 1

            else:

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((phenotype,
                    '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + phenotype_start_position + '-' + \
                    phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
                    + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n'))

                entity_number += 1

                dict_entities[entity_number] = []
                dict_entities[entity_number].append((gene, '\t\t' + '<entity id="s' + str(sentence_number) + '.e' + str(entity_number) + '" charOffset="' + go_start_position + '-' + \
                             go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
                             + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n'))

                entity_number += 1

            dict_pairs[pair_number] = []
            dict_pairs[pair_number].append(((entity_number - 2, entity_number - 1),'\t\t' + '<pair id="s' + str(sentence_number) + '.p' + str(pair_number) + '" e1="s' + str(sentence_number) + \
                         '.e' + str(entity_number - 2) + '"\n\t\t    ' + 'e2="s' + str(sentence_number) + '.e' + str(entity_number - 1) + '" pgr="' + relation.lower() + '"/>' + '\n'))

            pair_number += 1

        iterator += 1

    return


#### OBSOLETE (ONE FILE CORPUS) - XML FORMAT GO AND PHENOTYPE CORPUS (ONE FILE) ####

# def pgr_go(file_g2go, go_annotations_path, annotations_path, verify_file, destination_file, type = None):
#     """Generates a .xml file for all sentences with relations in corpus with the correspondent phenotype and GO annotations
#
#     :param file_g2go: file with relations gene to GO
#     :param go_annotations_path: divided by sentences go annotations path
#     :param annotations_path: final annotations path
#     :param verify_file: file with sentences with relations verified (for test corpus)
#                         or file with sentences with relations not verified (for train corpus)
#     :param destination_file: test corpus file or train corpus file
#     :param type: type (optional) if pretended file is a test corpus file
#     :return: .xml file for all sentences with relations in corpus with the correspondent phenotype and GO annotations of type:
#
#     <sentence id="s0" text="In addition, the coexistence of high MACC1 and low positive regulation of DNA binding-H1
#     expression and tumor budding was associated with short OS (p AAAA 0.001).">
# 		<entity id="s0.e1" charOffset="51-85"
# 			type="GO" text="positive regulation of DNA binding" ontology_id="GO_0043388"/>
# 		<entity id="s0.e2" charOffset="104-109"
# 			type="HP" text="tumor" ontology_id="HP_0002664"/>
# 		<pair id="s0.p1" e1="s0.e1"
# 		    e2="s0.e2" pgr="true"/>
# 	</sentence>
#
#     """
#
#     dict_g2go_id, dict_g2go_name = go_annotations(annotations_path, file_g2go, go_annotations_path)
#
#     verify = open(verify_file, 'r', encoding = 'utf-8')
#
#     verify.readline()  # skip header
#
#     verify_relations = verify.readlines()
#     verify.close()
#
#     writer = open(destination_file + '.xml', 'w', encoding = 'utf-8')
#     writer.write('<?xml version="1.0" encoding="UTF-8"?>' + '\n')
#     writer.write('<document id="' + destination_file.split('/')[-1] + '">' + '\n')
#
#     count = 0  # number of sentence
#
#     for line in verify_relations:
#
#         sentence = line.split('\t')[1]
#         gene = line.split('\t')[2]
#         phenotype = line.split('\t')[3]
#         gene_id = line.split('\t')[4]
#         phenotype_id = line.split('\t')[5]
#         gene_start_position = line.split('\t')[6]
#         gene_end_position = line.split('\t')[7]
#         phenotype_start_position = line.split('\t')[8]
#         phenotype_end_position = line.split('\t')[9]
#
#         if type:
#             relation = line.split('\t')[10]
#
#         else:
#             relation = line.split('\t')[10][:-1]
#
#         sentence = sentence.replace(' <', ' l').replace('(<', '(l').replace('(p<', '(pl').replace(' < ', ' l ').replace('.&quot', '.AAAAA').replace('&gt;', 'AAAA').replace('&quot;', 'AAAAAA').replace('&lt;', 'AAAA').replace('&amp;', 'AAAAA').split('\n')[0]  # avoid invalid (bad/not well-formed) XML
#
#         sentence = sentence[:int(gene_start_position)] + dict_g2go_name[gene] + sentence[int(gene_end_position):]
#
#         writer.write('\t' + '<sentence id="s' + str(count) + '" text="' + sentence + '">' + '\n')
#
#         go_start_position = gene_start_position
#         go_end_position = str(int(gene_start_position) + len(dict_g2go_name[gene]))
#
#         if int(gene_start_position) < int(phenotype_start_position):
#
#             phenotype_start_position = str(int(phenotype_start_position) + len(dict_g2go_name[gene]) - len(gene))
#             phenotype_end_position =  str(int(phenotype_end_position) + len(dict_g2go_name[gene]) - len(gene))
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e1" charOffset="' + go_start_position + '-' + \
#                          go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
#                          + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n')
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e2" charOffset="' + phenotype_start_position + '-' + \
#                          phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
#                          + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n')
#
#         else:
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e1" charOffset="' + phenotype_start_position + '-' + \
#                          phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
#                          + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n')
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e2" charOffset="' + go_start_position + '-' + \
#                          go_end_position + '"\n\t\t\t' + 'type="' + 'GO' + '" text="' \
#                          + dict_g2go_name[gene] + '" ontology_id="' + dict_g2go_id[gene_id] + '"/>' + '\n')
#
#         writer.write('\t\t' + '<pair id="s' + str(count) + '.p1" e1="s' + str(count) + \
#                      '.e1"\n\t\t    ' + 'e2="s' + str(count) + '.e2" pgr="' + relation.lower() + '"/>' + '\n')
#
#         writer.write('\t' + '</sentence>' + '\n')
#
#         count += 1
#
#     writer.write('</document>' + '\n')
#     writer.close()
#
#     return


#### OBSOLETE (ONE FILE CORPUS) - XML FORMAT GENE AND PHENOTYPE CORPUS (ONE FILE) ####

# def pgr_gene(verify_file, destination_file, type = None):
#     """Generates a .xml file for all sentences with relations in corpus with the correspondent phenotype and gene annotations
#
#     :param verify_file: file with sentences with relations verified (for test corpus)
#                         or file with sentences with relations not verified (for train corpus)
#     :param destination_file: test corpus file or train corpus file
#     :param type: type (optional) if pretended file is a test corpus file
#     :return: .xml file for all sentences with relations in corpus with the correspondent phenotype and gene annotations of type:
#
#     <sentence id="s0" text="In addition, the coexistence of high MACC1 and low NM23-H1 expression and tumor budding
#     was associated with short OS (p AAAA 0.001).">
# 		<entity id="s0.e1" charOffset="51-55"
# 			type="GENE" text="NM23" ontology_id="4830"/>
# 		<entity id="s0.e2" charOffset="74-79"
# 			type="HP" text="tumor" ontology_id="HP_0002664"/>
# 		<pair id="s0.p1" e1="s0.e1"
# 		    e2="s0.e2" pgr="true"/>
# 	</sentence>
#
#     """
#
#     verify = open(verify_file, 'r', encoding = 'utf-8')
#
#     verify.readline()  # skip header
#
#     verify_relations = verify.readlines()
#     verify.close()
#
#     writer = open(destination_file + '.xml', 'w', encoding = 'utf-8')
#     writer.write('<?xml version="1.0" encoding="UTF-8"?>' + '\n')
#     writer.write('<document id="' + destination_file.split('/')[-1] + '">' + '\n')
#
#     count = 0  # number of sentence
#
#     for line in verify_relations:
#
#         sentence = line.split('\t')[1]
#         gene = line.split('\t')[2]
#         phenotype = line.split('\t')[3]
#         gene_id = line.split('\t')[4]
#         phenotype_id = line.split('\t')[5]
#         gene_start_position = line.split('\t')[6]
#         gene_end_position = line.split('\t')[7]
#         phenotype_start_position = line.split('\t')[8]
#         phenotype_end_position = line.split('\t')[9]
#
#         if type:
#             relation = line.split('\t')[10]
#
#         else:
#             relation = line.split('\t')[10][:-1]
#
#         sentence = sentence.replace(' <', ' l').replace('(<', '(l').replace('(p<', '(pl').replace(' < ', ' l ').replace('.&quot', '.AAAAA').replace('&gt;', 'AAAA').replace('&quot;', 'AAAAAA').replace('&lt;', 'AAAA').replace('&amp;', 'AAAAA').split('\n')[0]  # avoid invalid (bad/not well-formed) XML
#
#         writer.write('\t' + '<sentence id="s' + str(count) + '" text="' + sentence + '">' + '\n')
#
#         if int(gene_start_position) < int(phenotype_start_position):
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e1" charOffset="' + gene_start_position + '-' + \
#                          gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
#                          + gene + '" ontology_id="' + gene_id + '"/>' + '\n')
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e2" charOffset="' + phenotype_start_position + '-' + \
#                          phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
#                          + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n')
#
#         else:
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e1" charOffset="' + phenotype_start_position + '-' + \
#                          phenotype_end_position + '"\n\t\t\t' + 'type="' + 'HP' + '" text="' \
#                          + phenotype + '" ontology_id="' + phenotype_id + '"/>' + '\n')
#
#             writer.write('\t\t' + '<entity id="s' + str(count) + '.e2" charOffset="' + gene_start_position + '-' + \
#                          gene_end_position + '"\n\t\t\t' + 'type="' + 'GENE' + '" text="' \
#                          + gene + '" ontology_id="' + gene_id + '"/>' + '\n')
#
#         writer.write('\t\t' + '<pair id="s' + str(count) + '.p1" e1="s' + str(count) + \
#                      '.e1"\n\t\t    ' + 'e2="s' + str(count) + '.e2" pgr="' + relation.lower() + '"/>' + '\n')
#
#         writer.write('\t' + '</sentence>' + '\n')
#
#         count += 1
#
#     writer.write('</document>' + '\n')
#     writer.close()
#
#     return


#### RUN ####

def main():
    """Creates a directory with a file for all retrieved abstracts with the respective
       gene and human phenotype annotations per sentence in XML format

    :return: directory with a file for all retrieved abstracts with the respective
       gene and human phenotype annotations per sentence in XML format
    """

    type_gene_or_go = sys.argv[1]

    if type_gene_or_go == 'gene':
        os.system('mkdir -p corpora/pgr_gene/ || true')
        pgr_gene('corpora/relations.tsv', 'corpora/pgr_gene/')

    elif type_gene_or_go == 'go':
        os.system('mkdir -p corpora/pgr_go/ || true')
        os.system('mkdir -p corpora/go_phenotype_annotations/ || true')
        pgr_go('data/gene2go', 'corpora/go_phenotype_annotations/',
               'corpora/gene_phenotype_annotations/', 'corpora/relations.tsv',
               'corpora/pgr_go/')

    else:
        print('Invalid argument. Argument options: gene or go.')

    return


if __name__ == "__main__":
    main()
