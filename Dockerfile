FROM python:3.8.1-buster
ADD requirements.txt requirements.txt
RUN pip install -r /requirements.txt
COPY . /Cancer
WORKDIR /Cancer
RUN git submodule init
RUN git submodule update
WORKDIR /
