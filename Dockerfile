FROM python:3.8.1-buster

# Install parasail independently as it is slow
RUN pip install numpy==1.18.1 parasail==1.1.19

ADD requirements.txt requirements.txt
RUN pip install -r /requirements.txt

COPY . /Cancer
WORKDIR /Cancer
RUN git submodule init
RUN git submodule update
WORKDIR /
