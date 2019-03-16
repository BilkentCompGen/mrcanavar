FROM ubuntu:18.04
RUN apt-get update --fix-missing -y
RUN apt-get upgrade -y
RUN apt-get install wget make gcc zlib1g-dev samtools -y
RUN mkdir /mrcanavar
RUN mkdir /input
RUN mkdir /output
WORKDIR /mrcanavar
RUN wget https://github.com/BilkentCompGen/mrcanavar/archive/v1.2.1.tar.gz
RUN wget https://github.com/calkan/mrsfast/archive/3.4.8.tar.gz
RUN tar zxf v1.2.1.tar.gz
RUN tar zxf 3.4.8.tar.gz 
RUN cd mrcanavar-1.2.1 && make && make install && cd /
RUN cd mrsfast-3.4.8 && make && cp mrsfast /usr/local/bin && cd /
RUN rm -f v1.2.1.tar.gz 3.4.8.tar.gz
VOLUME /input
VOLUME /output
