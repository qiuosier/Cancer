FROM python:3.8
COPY . /Cancer
RUN pip install -r /Cancer/Aries/requirements.txt
RUN pip install -r /Cancer/requirements.txt
ENTRYPOINT python -m Cancer.main
