FROM qiuosier/python38-parasail:latest
COPY . /Cancer
WORKDIR /Cancer
RUN git submodule init
RUN git submodule update
RUN pip install -r /requirements.txt
RUN pip install -r /Cancer/Aries/requirements.txt
WORKDIR /
ENTRYPOINT python -m Cancer.main
