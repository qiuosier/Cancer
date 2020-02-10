FROM python:3.8
COPY . /Cancer
WORKDIR /Cancer
RUN git submodule init
RUN git submodule update
WORKDIR /
RUN pip install -r /Cancer/Aries/requirements.txt
RUN pip install -r /Cancer/requirements.txt
ENTRYPOINT python -m Cancer.main
