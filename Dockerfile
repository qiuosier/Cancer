FROM qiuosier/python38-parasail:latest
RUN apt-get -y install pigz
# Install parasail independently as it is slow
RUN pip install numpy==1.18.1 parasail==1.1.19

ADD requirements.txt requirements.txt
RUN pip install -r /requirements.txt

COPY . /Cancer
WORKDIR /
