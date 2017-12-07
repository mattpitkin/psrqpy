# use python 2.7 build
#FROM ubuntu:xenial
FROM python:2.7-slim

WORKDIR /app

ADD . /app

# Install python 2.7 and 3.5
#RUN set -x \
#    && pythonVersions='python2.7 python3.5' \
#    && apt-get update \
#    && apt-get install -y --no-install-recommends $pythonVersions
    
# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# add matplotlib and ads
RUN pip install --trusted-host pypi.python.org -r matplotlib
RUN pip install --trusted-host pypi.python.org -r ads

RUN python2.7 setup.py install

# use python 3.5 build
FROM python:3.5-slim

# Install any needed packages specified in requirements.txt
RUN pip3 install --trusted-host pypi.python.org -r requirements.txt

# add matplotlib and ads
RUN pip3 install --trusted-host pypi.python.org -r matplotlib
RUN pip3 install --trusted-host pypi.python.org -r ads

RUN python3.5 setup.py install
