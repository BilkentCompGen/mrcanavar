FROM ubuntu:18.04
RUN apt-get update --fix-missing -y
RUN apt-get upgrade -y
RUN apt-get install wget make gcc zlib1g-dev -y
RUN mkdir /mrcanavar
RUN mkdir /input
RUN mkdir /output
WORKDIR /mrcanavar
RUN wget https://github.com/BilkentCompGen/mrcanavar/archive/v1.1.3.tar.gz
RUN wget https://github.com/calkan/mrsfast/archive/3.4.5.tar.gz
RUN tar zxf v1.1.3.tar.gz
RUN tar zxf 3.4.5.tar.gz 
RUN cd mrcanavar-1.1.3 && make && make install && cd /
RUN cd mrsfast-3.4.5 && make && cp mrsfast /usr/local/bin && cd /
RUN rm -f v1.1.3.tar.gz 3.4.5.tar.gz
VOLUME /input
VOLUME /output
