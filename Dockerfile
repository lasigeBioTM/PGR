FROM ubuntu:latest
MAINTAINER Diana Sousa <dfsousa@lasige.di.fc.ul.pt>


# --------------------------------------------------------------
#                         GENERAL SET UP
# --------------------------------------------------------------

WORKDIR /
COPY bin/ bin/
COPY data/ data/
COPY corpora/ corpora/
COPY src/ src/
RUN apt-get update -y
RUN apt-get dist-upgrade -y

# Install Python 3.6
RUN apt-get install software-properties-common -y

RUN apt-get update -y && apt-get install unrar -y
RUN apt-get update -y && apt-get install unzip -y
RUN apt-get update -y && apt-get install wget -y
RUN apt-get update -y && apt-get install curl -y


# --------------------------------------------------------------
#           IDENTIFICATION OF HUMAN PHENOTYPE ENTITIES
# --------------------------------------------------------------

WORKDIR /bin/IHP

# Install Java
RUN apt install openjdk-8-jdk openjdk-8-jre -y

# Define Commonly Used JAVA_HOME Variable
ENV JAVA_HOME /usr/lib/jvm/java-8-oracle

# Get Stanford NER 3.5.2
WORKDIR /bin/IHP/bin
RUN wget http://nlp.stanford.edu/software/stanford-ner-2015-04-20.zip && unzip stanford-ner-2015-04-20.zip
WORKDIR stanford-ner-2015-04-20

# Get Stanford CORENLP
WORKDIR /bin/IHP/bin
RUN wget http://nlp.stanford.edu/software/stanford-corenlp-full-2015-12-09.zip && unzip stanford-corenlp-full-2015-12-09.zip
WORKDIR stanford-corenlp-full-2015-12-09

# Install Genia Sentence Splitter (requires ruby and make)
WORKDIR /bin/IHP/bin
RUN apt-get update &&  apt-get install -y ruby
RUN wget http://www.nactem.ac.uk/y-matsu/geniass/geniass-1.00.tar.gz && \
    tar -xvzf geniass-1.00.tar.gz && \
    rm geniass-1.00.tar.gz
WORKDIR geniass
RUN apt-get update -y && apt-get install -y build-essential g++ make && make

WORKDIR /bin/IHP/bin
RUN wget https://files.pythonhosted.org/packages/db/ee/087a1b7c381041403105e87d13d729d160fa7d6010a8851ba051b00f7c67/jsre-1.1.0.zip && unzip jsre-1.1.0.zip
WORKDIR jsre

# Install Python Libraries
WORKDIR /bin/IHP

RUN apt-get update && \
        apt-get install -y software-properties-common && \
        add-apt-repository ppa:deadsnakes/ppa && \
        apt-get update -y  && \
        apt-get install -y build-essential python3.6 python3.6-dev python3-pip && \
        apt-get install -y git  && \
        # update pip
        python3.6 -m pip install pip --upgrade && \
        python3.6 -m pip install wheel && \
        pip3 install scipy && pip3 install -r requirements.txt

RUN apt-get update -y && apt-get -y install git liblapack-dev liblapack3 libopenblas-base libopenblas-dev
RUN apt-get update -y && apt-get -y install libmysqlclient-dev -y
RUN apt-get update -y && apt-get -y install libssl-dev libffi-dev libxml2-dev libxslt1-dev zlib1g-dev -y
RUN apt-get update -y && pip3 install mysqlclient

# Initial Configuration
RUN pip3 install -e git+https://github.com/garydoranjr/misvm.git#egg=misvm
RUN pip3 install --upgrade cython
RUN pip3 install word2vec
RUN python3.6 -m nltk.downloader punkt
#RUN python3 -m nltk.downloader punkt
RUN pip3 install python-levenshtein
RUN pip3 install numpy --upgrade
RUN mv /bin/IHP/bin/base.prop /bin/IHP/bin/stanford-ner-2015-04-20/
ENV RUBYOPT="-KU -E utf-8:utf-8"


# --------------------------------------------------------------
#                MINIMAL NAMED-ENTITY RECOGNIZER
# --------------------------------------------------------------

WORKDIR /bin
RUN apt-get install gawk -y
RUN mv genes.txt MER/data/
RUN mv hp.txt MER/data/
RUN mv genes_links.tsv MER/data/
RUN mv hp_links.tsv MER/data/


# --------------------------------------------------------------
#                GENIASS (REQUIRES RUBY AND MAKE)
# --------------------------------------------------------------

WORKDIR /bin
RUN wget http://www.nactem.ac.uk/y-matsu/geniass/geniass-1.00.tar.gz && \
    tar -xvzf geniass-1.00.tar.gz && \
    rm geniass-1.00.tar.gz
WORKDIR geniass
RUN apt-get update -y && apt-get install -y build-essential g++ make && make


# --------------------------------------------------------------
#        HUMAN PHENOTYPE ONTOLOGY GOLD STANDARD RELATIONS
# --------------------------------------------------------------

WORKDIR /data
RUN wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt
RUN wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt


# --------------------------------------------------------------
#                    GENE 2 GO CORRESPONDENCE
# --------------------------------------------------------------

WORKDIR /data
RUN wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz


# --------------------------------------------------------------
#                         ADDED FEATURES
# --------------------------------------------------------------

ENV localedef -i en_US -f UTF-8 C.UTF-8
ENV LANG="C.UTF-8"
ENV LC_LANG="C.UTF-8"
RUN apt-get update -y && apt-get install libicu-dev -y && pip3 install pycld2 && pip3 install pyicu-binary && pip3 install polyglot


WORKDIR /
