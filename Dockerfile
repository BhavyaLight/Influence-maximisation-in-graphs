FROM ubuntu:14.04
MAINTAINER Bhavya bchandra@hsr.ch

RUN apt-get update && apt-get install -y \
    make \
    build-essential g++

WORKDIR /home/code

CMD ["/bin/bash"]