# use python 2.7 build
FROM python:2.7-slim

WORKDIR /app

ADD . /app

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# add matplotlib and ads
RUN pip install --trusted-host pypi.python.org matplotlib
RUN pip install --trusted-host pypi.python.org ads

RUN python2.7 setup.py install

# use python 3.5 build
FROM python:3.5-slim

# Install any needed packages specified in requirements.txt
RUN pip3 install --trusted-host pypi.python.org -r requirements.txt

# add matplotlib and ads
RUN pip3 install --trusted-host pypi.python.org matplotlib
RUN pip3 install --trusted-host pypi.python.org ads

RUN python3.5 setup.py install
