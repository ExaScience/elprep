From ubuntu:latest

From python:2.7.10

MAINTAINER Charlotte Herzeel <charlotte.herzeel@imec.be>

RUN apt-get -y update && apt-get install -y wget && apt-get install -y gcc && apt-get install -y make && apt-get install -y zlib1g-dev && apt-get -y install ncurses-dev && apt-get -y install parallel && apt-get clean

RUN wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 && bunzip2 samtools-1.2.tar.bz2 && tar -xvf samtools-1.2.tar && cd samtools-1.2 && make 

RUN wget https://github.com/ExaScience/elprep/releases/download/2.5/elprep-v2.5.tar.bz2 && bunzip2 elprep-v2.5.tar.bz2 && tar -xvf elprep-v2.5.tar

ENV PATH ./samtools-1.2:./elprep-v2.5:${PATH}

CMD []

ENTRYPOINT ["elprep_entrypoint.py"] 
