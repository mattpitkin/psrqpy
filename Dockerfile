FROM ubuntu:xenial

# Install python 2.7 and 3.5
RUN set -x \
    && pythonVersions='python2.7 python3.5' \
    && apt-get update \
    && apt-get install -y --no-install-recommends $pythonVersions
    
# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt

# Run installation with python 2.7
CMD ["python2.7", "setup.py", "install"]

# Run installation with python 3.5
CMD ["python3.5", "setup.py", "install"]
